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
#ifndef itkSpectra1DSupportWindowImageFilter_h
#define itkSpectra1DSupportWindowImageFilter_h

#include <list>

#include "itkImageToImageFilter.h"

namespace itk
{

/** \class Spectra1DSupportWindowImageFilter
 * \brief Generate an image of local spectra computation support windows.
 *
 * The information from the input image is used to determine the output image
 * information. The pixel value of the input image is used to specify the
 * nominal number of lines on either side of the central FFT line to add to
 * the window. The nominal size of the 1D FFT is specified with SetFFTSize()
 *
 * The overlap between windows is specified with SetStep(). By default, the
 * Step is only one sample.
 *
 * \ingroup Ultrasound
 */
template <typename TInputImage>
class ITK_TEMPLATE_EXPORT Spectra1DSupportWindowImageFilter
  : public ImageToImageFilter<TInputImage,
                              Image<std::list<typename TInputImage::IndexType>, TInputImage::ImageDimension>>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(Spectra1DSupportWindowImageFilter);

  itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension);

  using InputImageType = TInputImage;
  using IndexType = typename InputImageType::IndexType;

  using OutputPixelType = std::list<IndexType>;
  using OutputImageType = Image<OutputPixelType, ImageDimension>;

  using FFT1DSizeType = unsigned int;

  /** Standard class type alias. */
  using Self = Spectra1DSupportWindowImageFilter;
  using Superclass = ImageToImageFilter<InputImageType, OutputImageType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  itkTypeMacro(Spectra1DSupportWindowImageFilter, ImageToImageFilter);
  itkNewMacro(Self);

  /** Set/Get the nominal size of the FFT.  This will be truncated at the
   * boundary of image. */
  itkGetConstMacro(FFT1DSize, FFT1DSizeType);
  itkSetMacro(FFT1DSize, FFT1DSizeType);

  /** Set/Get the number of samples between windows -- defaults to 1. */
  itkGetConstMacro(Step, SizeValueType);
  itkSetMacro(Step, SizeValueType);

protected:
  Spectra1DSupportWindowImageFilter();
  virtual ~Spectra1DSupportWindowImageFilter(){};

  using OutputImageRegionType = typename OutputImageType::RegionType;

  void
  GenerateOutputInformation() override;

  void
  DynamicThreadedGenerateData(const OutputImageRegionType & outputRegionForThread) override;
  void
  AfterThreadedGenerateData() override;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

private:
  FFT1DSizeType m_FFT1DSize;
  SizeValueType m_Step;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkSpectra1DSupportWindowImageFilter.hxx"
#endif

#endif // itkSpectra1DSupportWindowImageFilter_h
