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
#ifndef itkInverse1DFFTImageFilter_h
#define itkInverse1DFFTImageFilter_h

#include <complex>

#include "itkImageToImageFilter.h"
#include "itkImageRegionSplitterDirection.h"

namespace itk
{
/** \class Inverse1DFFTImageFilter
 * \brief Perform the Fast Fourier Transform, in the reverse direction, with
 * real output, but only along one dimension.
 *
 * \ingroup FourierTransform
 * \ingroup Ultrasound
 */
template< typename TInputImage, typename TOutputImage=Image< typename NumericTraits< typename TInputImage::PixelType >::ValueType, TInputImage::ImageDimension > >
class ITK_TEMPLATE_EXPORT Inverse1DFFTImageFilter:
  public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(Inverse1DFFTImageFilter);

  /** Standard class typedefs. */
  typedef TInputImage                          InputImageType;
  typedef TOutputImage                         OutputImageType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  typedef Inverse1DFFTImageFilter                               Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType > Superclass;
  typedef SmartPointer< Self >                                  Pointer;
  typedef SmartPointer< const Self >                            ConstPointer;

  itkStaticConstMacro(ImageDimension, unsigned int, InputImageType::ImageDimension );

  itkTypeMacro( Inverse1DFFTImageFilter, ImageToImageFilter );

  /** Customized object creation methods that support configuration-based
    * selection of FFT implementation.
    *
    * Default implementation is VnlFFT1D.
    */
  static Pointer New(void);

  /** Get the direction in which the filter is to be applied. */
  itkGetMacro(Direction, unsigned int);

  /** Set the direction in which the filter is to be applied. */
  itkSetClampMacro(Direction, unsigned int, 0, ImageDimension - 1);

protected:
  Inverse1DFFTImageFilter();
  virtual ~Inverse1DFFTImageFilter() {}

  void PrintSelf(std::ostream& os, Indent indent) const ITK_OVERRIDE;

  void GenerateInputRequestedRegion() ITK_OVERRIDE;
  void EnlargeOutputRequestedRegion(DataObject *output) ITK_OVERRIDE;

  void BeforeThreadedGenerateData() ITK_OVERRIDE;

  /** Override to return a splitter that does not split along the direction we
   * are performing the transform. */
  const ImageRegionSplitterBase* GetImageRegionSplitter() const ITK_OVERRIDE;

  /** Direction in which the filter is to be applied
   * this should be in the range [0,ImageDimension-1]. */
  unsigned int m_Direction;

private:
  ImageRegionSplitterDirection::Pointer m_ImageRegionSplitter;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#ifndef itkVnlInverse1DFFTImageFilter_h
#ifndef itkVnlInverse1DFFTImageFilter_hxx
#ifndef itkFFTWInverse1DFFTImageFilter_h
#ifndef itkFFTWInverse1DFFTImageFilter_hxx
#include "itkInverse1DFFTImageFilter.hxx"
#endif
#endif
#endif
#endif
#endif

#endif // itkInverse1DFFTImageFilter_h
