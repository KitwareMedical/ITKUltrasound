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

#include <algorithm>
#include <cmath>

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
  this->DynamicMultiThreadingOff();
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
AttenuationImageFilter<TInputImage, TOutputImage, TMaskImage>::VerifyPreconditions() const
{
  Superclass::VerifyPreconditions();

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
}

template <typename TInputImage, typename TOutputImage, typename TMaskImage>
void
AttenuationImageFilter<TInputImage, TOutputImage, TMaskImage>::BeforeThreadedGenerateData()
{
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

  // Initialize distance weights
  unsigned fourSigma = std::max(m_FixedEstimationDepth, 16u); // An arbitrary default.
  m_DistanceWeights.resize(fourSigma);
  float twoSigmaSquared = fourSigma * fourSigma / 8.0f;
  for (unsigned i = 0; i < fourSigma; ++i)
  {
    m_DistanceWeights[i] = 1.0f - std::exp(i * i / -twoSigmaSquared);
  }

  // Initialize iVars used in ComputeAttenuation()
  float nyquistFrequency = m_SamplingFrequencyMHz / 2;
  float numComponents = this->GetInput()->GetNumberOfComponentsPerPixel();
  m_FrequencyDelta = nyquistFrequency / numComponents;
  m_StartComponent = m_FrequencyBandStartMHz / m_FrequencyDelta;
  m_EndComponent = m_FrequencyBandEndMHz / m_FrequencyDelta;
  if (m_EndComponent == 0) // If m_FrequencyBandEndMHz is not set
  {
    m_EndComponent = numComponents - 1; // Use all components
  }
  m_ConsideredComponents = m_EndComponent - m_StartComponent + 1;
  m_ScanStepMM = this->GetInput()->GetSpacing()[m_Direction];
}

template <typename TInputImage, typename TOutputImage, typename TMaskImage>
void
AttenuationImageFilter<TInputImage, TOutputImage, TMaskImage>::ThreadedGenerateData(
  const OutputRegionType & regionForThread,
  ThreadIdType)
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

  thread_local std::vector<float> accumulatedWeight;
  accumulatedWeight.resize(regionForThread.GetSize(m_Direction));

  // make sure this is not recomputed in the inner loop
  const unsigned distanceWeightsSize = m_DistanceWeights.size();

  // do the work
  while (!it.IsAtEnd())
  {
    start = it.GetIndex();
    inclusionLength = 0;
    while (!it.IsAtEndOfLine())
    {
      // Advance until an inclusion is found
      InputIndexType index = it.GetIndex();
      bool           inside = ThreadedIsIncluded(index);
      if (inside)
      {
        if (inclusionLength == 0)
        {
          start = it.GetIndex(); // Mark the start
        }
        ++inclusionLength;
      }
      else if (inclusionLength > 0) // End of a segment
      {
        inclusionLength = 0; // Prepare for the next one

        // Stay at last pixel in the inclusion
        end = it.GetIndex();
        end[m_Direction] -= 1;

        // Adjust for pixel padding
        start[m_Direction] += m_PadLowerBounds;
        end[m_Direction] -= m_PadUpperBounds;

        if (start[m_Direction] < end[m_Direction]) // We need at least a pair of pixels to estimate attenuation
        {
          // Estimate attenuation for each inclusion pixel
          // by weighted average of pair-wise attenuations for all pairs
          while (start[m_Direction] <= end[m_Direction])
          {
            for (IndexValueType k = start[m_Direction] + 1; k <= end[m_Direction]; ++k)
            {
              unsigned pixelDistance = k - start[m_Direction];

              InputIndexType target = start;
              target[m_Direction] = k;
              float estimatedAttenuation = ComputeAttenuation(target, start);
              float weight = 1.0;                      // Weight for this pair's attenuation. 1 for large distances.
              if (pixelDistance < distanceWeightsSize) // If pixels are close, weight is lower than 1.
              {
                weight = m_DistanceWeights[pixelDistance];
              }

              // Update this pixel
              accumulatedWeight[start[m_Direction]] += weight;
              output->SetPixel(start, estimatedAttenuation * weight + output->GetPixel(start));

              // Update distant pair
              accumulatedWeight[k] += weight;
              output->SetPixel(target, estimatedAttenuation * weight + output->GetPixel(target));
            } // for k

            // Normalize output by accumulated weight
            output->SetPixel(start, output->GetPixel(start) / accumulatedWeight[start[m_Direction]]);
            accumulatedWeight[start[m_Direction]] = 0.0f; // reset for next next inclusion segment

            // Possibly eliminate negative attenuations
            if (!m_ConsiderNegativeAttenuations && output->GetPixel(start) < 0.0)
            {
              output->SetPixel(start, 0.0);
            }

            // Dynamically generate the output mask with values corresponding to input
            m_OutputMaskImage->SetPixel(start, inputMaskImage->GetPixel(start));

            ++start[m_Direction];
          } // while start<=end
        } // if start<end
      } // else !inside

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
  // Get RF spectra frequency bins at start and end pixel positions
  auto           input = this->GetInput();
  InputPixelType endSample = input->GetPixel(end);
  InputPixelType startSample = input->GetPixel(start);

  // Get distance between start and end pixel positions (assume mm units)
  const unsigned int pixelDistance = end[m_Direction] - start[m_Direction];
  float              distanceMM = pixelDistance * m_ScanStepMM;

  Eigen::Matrix<float, Eigen::Dynamic, 2> A(m_ConsideredComponents, 2);
  Eigen::Matrix<float, Eigen::Dynamic, 1> b(m_ConsideredComponents);
  for (unsigned i = 0; i < m_ConsideredComponents; i++)
  {
    A(i, 0) = 1;
    A(i, 1) = (1 + i + m_StartComponent) * m_FrequencyDelta;                                       // x_i = frequency
    b(i) = endSample[i + m_StartComponent] / (startSample[i + m_StartComponent] + itk::Math::eps); // y_i = ratio
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
