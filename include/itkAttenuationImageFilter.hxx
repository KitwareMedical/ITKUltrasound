/*=========================================================================
 *
 *  Copyright NumFOCUS
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkAttenuationImageFilter_hxx
#define itkAttenuationImageFilter_hxx

#include "itk_eigen.h"
#include ITK_EIGEN(Dense)
#include "itkMath.h"
#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkImageSink.h"
#include "itkImageRegionSplitterDirection.h"
#include "itkSimpleDataObjectDecorator.h"

namespace itk
{

template <typename TInputImage, typename TOutputImage, typename TMaskImage>
AttenuationImageFilter<TInputImage, TOutputImage, TMaskImage>::AttenuationImageFilter()
{
  this->SetNumberOfRequiredInputs(1);
  this->DynamicMultiThreadingOn();
}

template <typename TInputImage, typename TOutputImage, typename TMaskImage>
void
AttenuationImageFilter<TInputImage, TOutputImage, TMaskImage>::SetFixedEstimationDepthMM(const float distanceMM)
{
  this->SetFixedEstimationDepth(TransformPhysicalToPixelScanLineDistance(distanceMM));
}

template <typename TInputImage, typename TOutputImage, typename TMaskImage>
float
AttenuationImageFilter<TInputImage, TOutputImage, TMaskImage>::GetFixedEstimationDepthMM() const
{
  return this->TransformPixelToPhysicalScanLineDistance(this->GetFixedEstimationDepth());
}

template <typename TInputImage, typename TOutputImage, typename TMaskImage>
void
AttenuationImageFilter<TInputImage, TOutputImage, TMaskImage>::SetPadUpperBoundsMM(const float distanceMM)
{
  this->SetPadUpperBounds(TransformPhysicalToPixelScanLineDistance(distanceMM));
}

template <typename TInputImage, typename TOutputImage, typename TMaskImage>
float
AttenuationImageFilter<TInputImage, TOutputImage, TMaskImage>::GetPadUpperBoundsMM() const
{
  return this->TransformPixelToPhysicalScanLineDistance(this->GetPadUpperBounds());
}

template <typename TInputImage, typename TOutputImage, typename TMaskImage>
void
AttenuationImageFilter<TInputImage, TOutputImage, TMaskImage>::SetPadLowerBoundsMM(const float distanceMM)
{
  this->SetPadLowerBounds(TransformPhysicalToPixelScanLineDistance(distanceMM));
}

template <typename TInputImage, typename TOutputImage, typename TMaskImage>
float
AttenuationImageFilter<TInputImage, TOutputImage, TMaskImage>::GetPadLowerBoundsMM() const
{
  return this->TransformPixelToPhysicalScanLineDistance(this->GetPadLowerBounds());
}

template <typename TInputImage, typename TOutputImage, typename TMaskImage>
const ImageRegionSplitterBase *
AttenuationImageFilter<TInputImage, TOutputImage, TMaskImage>::GetImageRegionSplitter() const
{
  m_RegionSplitter->SetDirection(m_Direction);
  return m_RegionSplitter.GetPointer();
}

template <typename TInputImage, typename TOutputImage, typename TMaskImage>
void
AttenuationImageFilter<TInputImage, TOutputImage, TMaskImage>::BeforeThreadedGenerateData()
{
  // Verify inputs
  if (this->GetInputMaskImage() == nullptr)
  {
    itkExceptionMacro("Filter requires a mask image for inclusion estimates!");
  }
  else
  {
    m_ThreadedInputMaskImage = this->GetInputMaskImage();
  }

  if (this->GetDirection() >= ImageDimension)
  {
    itkExceptionMacro("Scan line direction must be a valid image dimension!");
  }

  if (this->GetSamplingFrequencyMHz() < itk::Math::eps)
  {
    itkExceptionMacro("RF sampling frequency was not set!");
  }

  Superclass::BeforeThreadedGenerateData();

  // Initialize metric image
  this->GetOutput()->Allocate();
  this->GetOutput()->FillBuffer(0.0f);

  // Initialize output mask image
  const MaskImageType * inputMaskImage = this->GetInputMaskImage();
  m_OutputMaskImage->CopyInformation(inputMaskImage);
  m_OutputMaskImage->SetRegions(inputMaskImage->GetLargestPossibleRegion());
  m_OutputMaskImage->Allocate();
  m_OutputMaskImage->FillBuffer(0.0f);
}

template <typename TInputImage, typename TOutputImage, typename TMaskImage>
void
AttenuationImageFilter<TInputImage, TOutputImage, TMaskImage>::DynamicThreadedGenerateData(
  const OutputRegionType & regionForThread)
{
  if (regionForThread.GetNumberOfPixels() == 0)
  {
    return;
  }

  const InputImageType * input = this->GetInput();
  OutputImageType *      output = this->GetOutput();
  const MaskImageType *  inputMaskImage = this->GetInputMaskImage();

  ImageLinearConstIteratorWithIndex<TInputImage> it(input, regionForThread);
  it.SetDirection(m_Direction);
  it.GoToBegin();

  unsigned int   inclusionLength;
  InputIndexType start, end;

  const float scanStepMM = input->GetSpacing()[m_Direction];

  // do the work
  while (!it.IsAtEnd())
  {
    inclusionLength = 0;
    while (!it.IsAtEndOfLine())
    {
      // Advance until an inclusion is found
      InputIndexType index = it.GetIndex();
      if (ThreadedIsIncluded(index) && inclusionLength == 0)
      {
        // Step into inclusion
        start = it.GetIndex();
      }
      else if (!ThreadedIsIncluded(index) && inclusionLength > 0)
      {
        // Stay at last pixel in the inclusion
        end = it.GetIndex();
        end[m_Direction] -= 1;

        // Adjust for pixel padding
        start[m_Direction] += m_PadLowerBounds;
        end[m_Direction] -= m_PadUpperBounds;

        // Estimate attenuation for each inclusion pixel
        // with respect to the last pixel in the inclusion
        while (start[m_Direction] < end[m_Direction])
        {
          // If no fixed estimation depth is set or too few pixels remain for fixed-depth estimation
          // then take the attenuation between the given pixel and the last pixel in the inclusion
          InputIndexType target = end;

          if (m_FixedEstimationDepth != 0 && end[m_Direction] - start[m_Direction] > m_FixedEstimationDepth)
          {
            target[m_Direction] = start[m_Direction] + m_FixedEstimationDepth;
          }

          float estimatedAttenuation = ComputeAttenuation(target, start);

          // Update the corresponding pixel in the metric image
          if (estimatedAttenuation > 0.0 || m_ConsiderNegativeAttenuations)
          {
            // Assignment is thread-safe for internal operations because
            // each thread writes only within its allotted output region.
            output->SetPixel(start, estimatedAttenuation);
          }

          // Dynamically generate the output mask with values corresponding to input
          m_OutputMaskImage->SetPixel(start, inputMaskImage->GetPixel(start));

          ++start[m_Direction];
        }

        inclusionLength = 0;
      }

      if (ThreadedIsIncluded(index))
      {
        ++inclusionLength;
      }

      ++it;
    }
    it.NextLine();
  }
};

template <typename TInputImage, typename TOutputImage, typename TMaskImage>
typename AttenuationImageFilter<TInputImage, TOutputImage, TMaskImage>::OutputPixelType
AttenuationImageFilter<TInputImage, TOutputImage, TMaskImage>::ComputeAttenuation(const InputIndexType & end,
                                                                                  const InputIndexType & start) const
{
  const float nyquistFrequency = m_SamplingFrequencyMHz / 2;

  // Number and width of RF spectra frequency bins over the range (0, nyquist_frequency]
  const unsigned int numComponents = this->GetInput()->GetNumberOfComponentsPerPixel();
  const float        frequencyDelta = nyquistFrequency / numComponents;

  // Frequency band to consider for attenuation
  const unsigned int startComponent = m_FrequencyBandStartMHz / frequencyDelta;
  const unsigned int endComponent = m_FrequencyBandEndMHz / frequencyDelta;
  const unsigned int consideredComponents = endComponent - startComponent + 1;

  // Get RF spectra frequency bins at start and end pixel positions
  auto           input = this->GetInput();
  InputPixelType endSample = input->GetPixel(end);
  InputPixelType startSample = input->GetPixel(start);

  // Get distance between start and end pixel positions (assume mm units)
  const float        scanStepMM = input->GetSpacing()[m_Direction];
  const unsigned int pixelDistance = end[m_Direction] - start[m_Direction];
  float              distanceMM = pixelDistance * scanStepMM;

  Eigen::Matrix<float, Eigen::Dynamic, 2> A(consideredComponents, 2);
  Eigen::Matrix<float, Eigen::Dynamic, 1> b(consideredComponents);
  for (unsigned i = 0; i < consideredComponents; i++)
  {
    A(i, 0) = 1;
    A(i, 1) = (1 + i + startComponent) * frequencyDelta;                                       // x_i = frequency
    b(i) = endSample[i + startComponent] / (startSample[i + startComponent] + itk::Math::eps); // y_i = ratio
  }

  // from https://eigen.tuxfamily.org/dox/group__LeastSquares.html
  Eigen::Matrix<float, 1, 2> lineFit = A.householderQr().solve(b);
  float                      frequencySlope = -lineFit(1); // we expect attenuation to increase with frequency

  // https://www.electronics-notes.com/articles/basic_concepts/decibel/neper-to-db-conversion.php
  // Neper to dB conversion: 1Np = 20 log10e dB, approximately 1Np = 8.6858896 dB
  float neper = 20 * itk::Math::log10e;

  return 10 * neper * frequencySlope / distanceMM; // 10 converts mm into cm
}

template <typename TInputImage, typename TOutputImage, typename TMaskImage>
float
AttenuationImageFilter<TInputImage, TOutputImage, TMaskImage>::TransformPhysicalToPixelScanLineDistance(
  float distanceMM) const
{
  if (distanceMM < 0)
  {
    itkExceptionMacro("Expected nonnegative spatial distance!");
  }

  auto input = this->GetInput();
  if (input == nullptr)
  {
    itkExceptionMacro("Tried to translate spatial distance to pixel distance without reference input image!");
  }

  const float scanStepMM = input->GetSpacing()[m_Direction];
  float       distanceInPixels = distanceMM / scanStepMM;
  return static_cast<unsigned int>(std::round(distanceInPixels));
}

template <typename TInputImage, typename TOutputImage, typename TMaskImage>
float
AttenuationImageFilter<TInputImage, TOutputImage, TMaskImage>::TransformPixelToPhysicalScanLineDistance(
  unsigned int distance) const
{
  auto input = this->GetInput();
  if (input == nullptr)
  {
    itkExceptionMacro("Tried to translate spatial distance to pixel distance without reference input image!");
  }

  const float scanStepMM = input->GetSpacing()[m_Direction];
  return scanStepMM * distance;
}

template <typename TInputImage, typename TOutputImage, typename TMaskImage>
bool
AttenuationImageFilter<TInputImage, TOutputImage, TMaskImage>::ThreadedIsIncluded(InputIndexType index) const
{
  auto maskValue = m_ThreadedInputMaskImage->GetPixel(index);
  return (m_LabelValue == 0 && maskValue > 0) || (m_LabelValue > 0 && maskValue == m_LabelValue);
}

template <typename TInputImage, typename TOutputImage, typename TMaskImage>
void
AttenuationImageFilter<TInputImage, TOutputImage, TMaskImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Image axis representing RF scanline: " << this->GetDirection() << std::endl;
  os << indent << "Label value: " << static_cast<unsigned int>(this->GetLabelValue()) << std::endl;
  os << indent << "Sampling frequency (MHz): " << this->GetSamplingFrequencyMHz() << std::endl;
  os << indent << "Frequency band: [" << this->GetFrequencyBandStartMHz() << "," << this->GetFrequencyBandEndMHz()
     << "]" << std::endl;
  os << indent << "Consider negative attenuations: " << (this->GetConsiderNegativeAttenuations() ? "Yes" : "No")
     << std::endl;
  os << indent << "Fixed estimation distance: " << this->GetFixedEstimationDepthMM()
     << "mm == " << this->GetFixedEstimationDepth() << "px" << std::endl;
  os << indent << "Inclusion padding on scanline: Lower: " << this->GetPadLowerBoundsMM()
     << "mm == " << this->GetPadLowerBounds() << "px ; Upper: " << this->GetPadUpperBoundsMM()
     << "mm == " << this->GetPadUpperBounds() << " px" << std::endl;
}
} // end namespace itk
#endif // itkAttenuationImageFilter_hxx
