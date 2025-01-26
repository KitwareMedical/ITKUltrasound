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
#ifndef itkSpectra1DNormalizeImageFilter_h
#define itkSpectra1DNormalizeImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkFrequencyDomain1DFilterFunction.h"

namespace itk
{
/** \class Spectra1DNormalizeImageFilter
 * \brief Normalize (divide) a spectral image by reference line spectra.
 *
 * Single transducer line is the fastest-varying dimension.
 *
 * \ingroup Ultrasound
 */
template <typename TInputImage, typename TReferenceImage>
class ITK_TEMPLATE_EXPORT Spectra1DNormalizeImageFilter : public ImageToImageFilter<TInputImage, TInputImage>
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(Spectra1DNormalizeImageFilter);

  /** Standard class type alias. */
  using InputImageType = TInputImage;
  using OutputImageType = TInputImage;
  using OutputImageRegionType = typename OutputImageType::RegionType;
  using ReferenceImageType = TReferenceImage;

  itkStaticConstMacro(ImageDimension, unsigned int, InputImageType::ImageDimension);

  using Self = Spectra1DNormalizeImageFilter;
  using Superclass = ImageToImageFilter<InputImageType, OutputImageType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  itkOverrideGetNameOfClassMacro(Spectra1DNormalizeImageFilter);
  itkNewMacro(Self);

  void
  SetReferenceImage(ReferenceImageType * referenceImage)
  {
    this->SetInput("ReferenceImage", referenceImage);
  }

protected:
  Spectra1DNormalizeImageFilter() { this->AddRequiredInputName("ReferenceImage", 1); }
  ~Spectra1DNormalizeImageFilter() override = default;

  void
  GenerateInputRequestedRegion() override;

  void
  DynamicThreadedGenerateData(const OutputImageRegionType & outputRegionForThread) override;
};
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkSpectra1DNormalizeImageFilter.hxx"
#endif

#endif // itkSpectra1DNormalizeImageFilter_h
