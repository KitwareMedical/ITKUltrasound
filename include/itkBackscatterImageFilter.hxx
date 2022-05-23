/*=========================================================================
 *
 *  Copyright NumFOCUS
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         https://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkBackscatterImageFilter_hxx
#define itkBackscatterImageFilter_hxx

#include <algorithm>
#include <cmath>

#include "itk_eigen.h"
#include ITK_EIGEN(Dense)
#include "itkMath.h"
#include "itkImageScanlineConstIterator.h"

namespace itk
{

template <typename TInputImage, typename TOutputImage>
BackscatterImageFilter<TInputImage, TOutputImage>::BackscatterImageFilter()
{
  this->SetNumberOfRequiredInputs(1);
  this->SetNumberOfRequiredOutputs(3);
  this->SetNthOutput(0, this->MakeOutput(0));
  this->SetNthOutput(1, this->MakeOutput(1));
  this->SetNthOutput(2, this->MakeOutput(2));
}

template <typename TInputImage, typename TOutputImage>
void
BackscatterImageFilter<TInputImage, TOutputImage>::VerifyPreconditions() const
{
  Superclass::VerifyPreconditions();

  if (this->GetSamplingFrequencyMHz() < itk::Math::eps)
  {
    itkExceptionMacro("RF sampling frequency was not set!");
  }

  if (this->GetFrequencyBandEndMHz() <= this->GetFrequencyBandStartMHz())
  {
    itkExceptionMacro(<< "FrequencyBandStart must be less than FrequencyBandEnd!"
                      << "\n  FrequencyBandStartMHz: " << this->GetFrequencyBandStartMHz()
                      << "\n  FrequencyBandEndMHz: " << this->GetFrequencyBandEndMHz());
  }
}

template <typename TInputImage, typename TOutputImage>
void
BackscatterImageFilter<TInputImage, TOutputImage>::BeforeThreadedGenerateData()
{
  Superclass::BeforeThreadedGenerateData();

  // Initialize metric image
  this->GetOutput(0)->Allocate();
  this->GetOutput(1)->Allocate();
  this->GetOutput(2)->Allocate();

  // Initialize iVars used in ComputeBackscatter()
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
}

template <typename TInputImage, typename TOutputImage>
void
BackscatterImageFilter<TInputImage, TOutputImage>::DynamicThreadedGenerateData(const OutputRegionType & regionForThread)
{
  if (regionForThread.GetNumberOfPixels() == 0)
  {
    return;
  }

  const InputImageType * input = this->GetInput();
  OutputImageType *      averageOutput = this->GetOutput(0);
  OutputImageType *      slopeOutput = this->GetOutput(1);
  OutputImageType *      interceptOutput = this->GetOutput(2);

  ImageScanlineConstIterator<TInputImage> it(input, regionForThread);
  it.GoToBegin();

  // do the work
  while (!it.IsAtEnd())
  {
    while (!it.IsAtEndOfLine())
    {
      InputIndexType index = it.GetIndex();

      // Pass output pixel references directly
      ComputeBackscatter(
        index, averageOutput->GetPixel(index), slopeOutput->GetPixel(index), interceptOutput->GetPixel(index));

      ++it;
    }
    it.NextLine();
  }
};

template <typename TInputImage, typename TOutputImage>
void
BackscatterImageFilter<TInputImage, TOutputImage>::ComputeBackscatter(const InputIndexType & index,
                                                                      OutputPixelType &      average,
                                                                      OutputPixelType &      slope,
                                                                      OutputPixelType &      intercept) const
{
  using ScalarType = typename NumericTraits<OutputPixelType>::ValueType;

  // Get RF spectra frequency bins at start and end pixel positions
  auto           input = this->GetInput();
  InputPixelType sample = input->GetPixel(index);
  ScalarType     sum = 0;

  Eigen::Matrix<float, Eigen::Dynamic, 2> A(m_ConsideredComponents, 2);
  Eigen::Matrix<float, Eigen::Dynamic, 1> b(m_ConsideredComponents);
  for (unsigned i = 0; i < m_ConsideredComponents; i++)
  {
    A(i, 0) = 1;
    A(i, 1) = (1 + i + m_StartComponent) * m_FrequencyDelta; // x_i = frequency
    b(i) = sample[i + m_StartComponent];                     // y_i = intensity
    sum += sample[i + m_StartComponent];
  }

  average = sum / m_ConsideredComponents;

  // from https://eigen.tuxfamily.org/dox/group__LeastSquares.html
  Eigen::Matrix<float, 1, 2> lineFit = A.householderQr().solve(b);
  slope = -lineFit(1);
  intercept = lineFit(0);
}

template <typename TInputImage, typename TOutputImage>
void
BackscatterImageFilter<TInputImage, TOutputImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Sampling frequency (MHz): " << this->GetSamplingFrequencyMHz() << std::endl;
  os << indent << "Analysis frequency band: [" << this->GetFrequencyBandStartMHz() << ","
     << this->GetFrequencyBandEndMHz() << "]" << std::endl;
}
} // end namespace itk
#endif // itkBackscatterImageFilter_hxx
