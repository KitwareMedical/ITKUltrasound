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
#ifndef itkBModeImageFilter_h
#define itkBModeImageFilter_h

#include "itkAddImageFilter.h"
#include "itkComplexToModulusImageFilter.h"
#include "itkConstantPadImageFilter.h"
#include "itkRegionFromReferenceImageFilter.h"
#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkLog10ImageFilter.h"

#include "itkAnalyticSignalImageFilter.h"

namespace itk
{

/**
 * \class BModeImageFilter
 *
 * \brief Create an ultrasound B-Mode (Brightness-Mode) image from raw
 * "RF" data.  The RF's envelope is calculated from the analytic signal and
 * logarithmic intensity transform is applied.
 *
 * Use SetDirection() to define the axis of propagation.
 *
 * Use SetFrequencyFilter() to add a filtering step before the analytic
 * signal computation.
 *
 * \sa AnalyticSignalImageFilter
 *
 * \ingroup Ultrasound
 */
template <typename TInputImage,
          typename TOutputImage = TInputImage,
          typename TComplexImage = Image<std::complex<typename TInputImage::PixelType>, TInputImage::ImageDimension>>
class ITK_TEMPLATE_EXPORT BModeImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class type alias.   */
  using Self = BModeImageFilter;
  using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** The type of input image.   */
  using InputImageType = TInputImage;

  /** Dimension of the input and output images. */
  itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension);

  /** Typedef support for the input image scalar value type. */
  using InputPixelType = typename InputImageType::PixelType;

  /** The type of output image.   */
  using OutputImageType = TOutputImage;

  /** Typedef support for the output image scalar value type. */
  using OutputPixelType = typename OutputImageType::PixelType;

  /** Typedef of the image used for internal computations that has
   * std::complex pixels. */
  using ComplexImageType = TComplexImage;

  /** Other convenient type alias   */
  using InputRegionType = typename InputImageType::RegionType;
  using InputSizeType = typename InputImageType::SizeType;
  using InputIndexType = typename InputImageType::IndexType;

  /** Run-time type information (and related methods) */
  itkOverrideGetNameOfClassMacro(BModeImageFilter);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  using FrequencyFilterType = FrequencyDomain1DImageFilter<ComplexImageType, ComplexImageType>;

  /** Set the direction in which the envelope is to be calculated. */
  virtual void
  SetDirection(unsigned int direction)
  {
    this->m_AnalyticFilter->SetDirection(direction);
    this->Modified();
  }

  /** Get the direction in which the envelope is to be calculated. */
  virtual unsigned int
  GetDirection() const
  {
    return m_AnalyticFilter->GetDirection();
  }

  void
  SetFrequencyFilter(FrequencyFilterType * filter)
  {
    m_AnalyticFilter->SetFrequencyFilter(filter);
  }

protected:
  BModeImageFilter();
  ~BModeImageFilter() override = default;

  virtual void
  PrintSelf(std::ostream & os, Indent indent) const override;

  virtual void
  GenerateData() override;

  // These behave like their analogs in Forward1DFFTImageFilter.
  virtual void
  GenerateInputRequestedRegion() override;
  virtual void
  EnlargeOutputRequestedRegion(DataObject * output) override;

  /** Component filters. */
  using AnalyticType = AnalyticSignalImageFilter<InputImageType, ComplexImageType>;
  using ComplexToModulusType = ComplexToModulusImageFilter<typename AnalyticType::OutputImageType, OutputImageType>;
  using PadType = ConstantPadImageFilter<InputImageType, InputImageType>;
  using AddConstantType = AddImageFilter<InputImageType, InputImageType>;
  using LogType = Log10ImageFilter<InputImageType, OutputImageType>;
  using ROIType = RegionFromReferenceImageFilter<OutputImageType, OutputImageType>;

private:
  BModeImageFilter(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented

  typename AnalyticType::Pointer         m_AnalyticFilter;
  typename ComplexToModulusType::Pointer m_ComplexToModulusFilter;
  typename PadType::Pointer              m_PadFilter;
  typename AddConstantType::Pointer      m_AddConstantFilter;
  typename LogType::Pointer              m_LogFilter;
  typename ROIType::Pointer              m_ROIFilter;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkBModeImageFilter.hxx"
#endif

#endif // itkBModeImageFilter_h
