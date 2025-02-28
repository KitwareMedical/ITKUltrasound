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
#ifndef itkAnalyticSignalImageFilter_h
#define itkAnalyticSignalImageFilter_h

#include <complex>

#include "itkComplexToComplex1DFFTImageFilter.h"
#include "itkForward1DFFTImageFilter.h"
#include "itkFrequencyDomain1DImageFilter.h"

namespace itk
{
/** \class AnalyticSignalImageFilter
 * \brief Generates the analytic signal from one direction of an image.
 *
 * This filter generates the complex valued analytic signal along one direction
 * of an image.  This input is a real valued image, and the output is a complex
 * image.
 *
 * The analytic signal is given by
 *
 * f_a(x) = f(x) - i f_H(x)
 *
 * Where i is the square root of one and f_H(x) is the Hibert transform of f(x).
 *
 * Since the Hilbert transform in the Fourier domain is
 *
 * F_H(k) = F(k) i sign(k),
 *
 * f_a(x) is calculated by
 *
 * f_a(x) = F^{-1}( F(k) 2 U(k) )
 *
 * where U(k) is the unit step function.
 *
 * \ingroup FourierTransform
 * \ingroup Ultrasound
 */
template <typename TInputImage, typename TOutputImage>
class ITK_TEMPLATE_EXPORT AnalyticSignalImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(AnalyticSignalImageFilter);

  /** Standard class type alias. */
  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;
  using OutputImageRegionType = typename OutputImageType::RegionType;

  static constexpr unsigned int ImageDimension = InputImageType::ImageDimension;

  using Self = AnalyticSignalImageFilter;
  using Superclass = ImageToImageFilter<InputImageType, OutputImageType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  itkOverrideGetNameOfClassMacro(AnalyticSignalImageFilter);
  itkNewMacro(Self);

  using FrequencyFilterType = FrequencyDomain1DImageFilter<OutputImageType, OutputImageType>;

  /** Get the direction in which the filter is to be applied. */
  virtual unsigned int
  GetDirection() const
  {
    return this->m_FFTRealToComplexFilter->GetDirection();
  }

  /** Set the direction in which the filter is to be applied. */
  virtual void
  SetDirection(const unsigned int direction)
  {
    if (this->m_FFTRealToComplexFilter->GetDirection() != direction)
    {
      this->m_FFTRealToComplexFilter->SetDirection(direction);
      this->m_FFTComplexToComplexFilter->SetDirection(direction);
      if (this->m_FrequencyFilter.IsNotNull())
      {
        this->m_FrequencyFilter->SetDirection(direction);
      }
      this->Modified();
    }
  }

  virtual void
  SetFrequencyFilter(FrequencyFilterType * filter)
  {
    if (filter != this->m_FrequencyFilter.GetPointer())
    {
      this->m_FrequencyFilter = filter;
      this->m_FrequencyFilter->SetDirection(this->GetDirection());
      this->Modified();
    }
  }

protected:
  AnalyticSignalImageFilter();
  ~AnalyticSignalImageFilter() override = default;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  // These behave like their analogs in Forward1DFFTImageFilter.
  void
  GenerateInputRequestedRegion() override;
  void
  EnlargeOutputRequestedRegion(DataObject * output) override;

  void
  GenerateData() override;

  using FFTRealToComplexType = Forward1DFFTImageFilter<InputImageType, OutputImageType>;
  typename FFTRealToComplexType::Pointer m_FFTRealToComplexFilter;

  using FFTComplexToComplexType = ComplexToComplex1DFFTImageFilter<OutputImageType, OutputImageType>;
  typename FFTComplexToComplexType::Pointer m_FFTComplexToComplexFilter;

private:
  typename FrequencyFilterType::Pointer m_FrequencyFilter;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkAnalyticSignalImageFilter.hxx"
#endif

#endif // itkAnalyticSignalImageFilter_h
