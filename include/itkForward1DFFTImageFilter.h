/*=========================================================================
 *
 *  Copyright Insight Software Consortium
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
#ifndef itkForward1DFFTImageFilter_h
#define itkForward1DFFTImageFilter_h

#include <complex>

#include "itkImageToImageFilter.h"
#include "itkImageRegionSplitterDirection.h"

namespace itk
{
/** \class Forward1DFFTImageFilter
 * \brief Perform the Fast Fourier Transform, in the forward direction, with
 * real inputs, but only along one dimension.
 *
 * \ingroup FourierTransform
 * \ingroup Ultrasound
 */
template< typename TInputImage, typename TOutputImage=Image< std::complex< typename TInputImage::PixelType >, TInputImage::ImageDimension > >
class ITK_TEMPLATE_EXPORT Forward1DFFTImageFilter:
  public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(Forward1DFFTImageFilter);

  /** Standard class typedefs. */
  typedef TInputImage                                           InputImageType;
  typedef TOutputImage                                          OutputImageType;
  typedef typename OutputImageType::RegionType                  OutputImageRegionType;

  typedef Forward1DFFTImageFilter                               Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType > Superclass;
  typedef SmartPointer< Self >                                  Pointer;
  typedef SmartPointer< const Self >                            ConstPointer;

  itkStaticConstMacro( ImageDimension, unsigned int, InputImageType::ImageDimension );

  itkTypeMacro( Forward1DFFTImageFilter, ImageToImageFilter );

  /** Customized object creation methods that support configuration-based
    * selection of FFT implementation.
    *
    * Default implementation is VnlFFT1D.
    */
  static Pointer New();

  /** Get the direction in which the filter is to be applied. */
  itkGetMacro(Direction, unsigned int);

  /** Set the direction in which the filter is to be applied. */
  itkSetClampMacro(Direction, unsigned int, 0, ImageDimension - 1);

protected:
  Forward1DFFTImageFilter();
  virtual ~Forward1DFFTImageFilter() {}

  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;

  void GenerateInputRequestedRegion() ITK_OVERRIDE;
  void EnlargeOutputRequestedRegion(DataObject *output) ITK_OVERRIDE;

  void BeforeThreadedGenerateData() ITK_OVERRIDE;

  /** Override to return a splitter that does not split along the direction we
   * are performing the transform. */
  virtual const ImageRegionSplitterBase* GetImageRegionSplitter() const ITK_OVERRIDE;

private:
  ImageRegionSplitterDirection::Pointer m_ImageRegionSplitter;

  /** Direction in which the filter is to be applied
   * this should be in the range [0,ImageDimension-1]. */
  unsigned int m_Direction;

};
}

#ifndef ITK_MANUAL_INSTANTIATION
#ifndef itkVnlForward1DFFTImageFilter_h
#ifndef itkVnlForward1DFFTImageFilter_hxx
#ifndef itkFFTWForward1DFFTImageFilter_h
#ifndef itkFFTWForward1DFFTImageFilter_hxx
#include "itkForward1DFFTImageFilter.hxx"
#endif
#endif
#endif
#endif
#endif

#endif // itkForward1DFFTImageFilter_h
