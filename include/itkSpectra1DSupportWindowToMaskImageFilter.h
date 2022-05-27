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
#ifndef itkSpectra1DSupportWindowToMaskImageFilter_h
#define itkSpectra1DSupportWindowToMaskImageFilter_h

#include "itkImageToImageFilter.h"

namespace itk
{

/** \class Spectra1DSupportWindowToMaskImageFilter
 * \brief Generate a mask image from the support window at a given index.
 *
 * \ingroup Ultrasound
 */
template <typename TInputImage, typename TOutputImage>
class ITK_TEMPLATE_EXPORT Spectra1DSupportWindowToMaskImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension);

  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;

  using IndexType = typename InputImageType::IndexType;
  using OutputPixelType = typename OutputImageType::PixelType;

  /** Standard class type alias. */
  using Self = Spectra1DSupportWindowToMaskImageFilter;
  using Superclass = ImageToImageFilter<InputImageType, OutputImageType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  itkTypeMacro(Spectra1DSupportWindowToMaskImageFilter, ImageToImageFilter);
  itkNewMacro(Self);

  /** Set/Get the index of the support window to create the mask for. */
  itkGetConstReferenceMacro(MaskIndex, IndexType);
  itkSetMacro(MaskIndex, IndexType);

  /** Set/Get the value to consider as "background". Defaults to zero. */
  itkSetMacro(BackgroundValue, OutputPixelType);
  itkGetConstMacro(BackgroundValue, OutputPixelType);

  /** Set/Get the value in the image to consider as "foreground". Defaults to
   * maximum value of the OutputPixelType. */
  itkSetMacro(ForegroundValue, OutputPixelType);
  itkGetConstMacro(ForegroundValue, OutputPixelType);


protected:
  Spectra1DSupportWindowToMaskImageFilter();
  ~Spectra1DSupportWindowToMaskImageFilter() override = default;

  virtual void
  GenerateData() override;

private:
  Spectra1DSupportWindowToMaskImageFilter(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented

  IndexType m_MaskIndex;

  OutputPixelType m_BackgroundValue;
  OutputPixelType m_ForegroundValue;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkSpectra1DSupportWindowToMaskImageFilter.hxx"
#endif

#endif // itkSpectra1DSupportWindowToMaskImageFilter_h
