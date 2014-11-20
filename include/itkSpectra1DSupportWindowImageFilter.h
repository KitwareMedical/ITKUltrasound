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
#ifndef __itkSpectra1DSupportWindowImageFilter_h
#define __itkSpectra1DSupportWindowImageFilter_h

#include <list>

#include "itkImageToImageFilter.h"

namespace itk
{

/** \class Spectra1DSupportWindowImageFilter
 * \brief Generate an image of local spectra computation support windows.
 *
 * \ingroup FFT1D
 */
template< class TInputImage >
class Spectra1DSupportWindowImageFilter:
  public ImageToImageFilter< TInputImage, TInputImage >
{
public:
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension );

  typedef TInputImage                        InputImageType;

  typedef typename InputImageType::IndexType IndexType;
  struct FFT1DRegionType
    {
    IndexType     Index;
    SizeValueType Size;
    };
  typedef std::list< FFT1DRegionType >             OutputPixelType;
  typedef Image< OutputPixelType, ImageDimension > OutputImageType;

  /** Standard class typedefs. */
  typedef Spectra1DSupportWindowImageFilter                     Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType > Superclass;
  typedef SmartPointer< Self >                                  Pointer;
  typedef SmartPointer< const Self >                            ConstPointer;

  itkTypeMacro( Spectra1DSupportWindowImageFilter, ImageToImageFilter );
  itkNewMacro( Self );


protected:
  Spectra1DSupportWindowImageFilter();
  virtual ~Spectra1DSupportWindowImageFilter() {};

private:
  Spectra1DSupportWindowImageFilter( const Self & ); // purposely not implemented
  void operator=( const Self & ); // purposely not implemented
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSpectra1DSupportWindowImageFilter.hxx"
#endif

#endif // __itkSpectra1DSupportWindowImageFilter_h
