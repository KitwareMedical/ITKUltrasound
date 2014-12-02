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
#ifndef __itkSpectra1DImageFilter_h
#define __itkSpectra1DImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkDefaultConvertPixelTraits.h"

#include "vnl/algo/vnl_fft_base.h"
#include "vnl/algo/vnl_fft_1d.h"

#include <utility>

#include "itkSpectra1DSupportWindowImageFilter.h"

namespace itk
{

/** \class Spectra1DImageFilter
 * \brief Generate an image of local spectra.
 *
 * This image takes in the input image and image that has indexes of the local
 * lines used to compute the local spectra.
 *
 * \ingroup Ultrasound
 *
 * \sa Spectra1DSupportWindowImageFilter
 */
template< typename TInputImage, typename TSupportWindowImage, typename TOutputImage >
class Spectra1DImageFilter:
  public ImageToImageFilter< TInputImage,
                             TOutputImage >
{
public:
  itkStaticConstMacro( ImageDimension, unsigned int, TInputImage::ImageDimension );

  typedef TInputImage                              InputImageType;
  typedef TSupportWindowImage                      SupportWindowImageType;
  typedef TOutputImage                             OutputImageType;

  typedef typename DefaultConvertPixelTraits< typename OutputImageType::PixelType >::ComponentType
    ScalarType;

  /** Standard class typedefs. */
  typedef Spectra1DImageFilter                                  Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType > Superclass;
  typedef SmartPointer< Self >                                  Pointer;
  typedef SmartPointer< const Self >                            ConstPointer;

  itkTypeMacro( Spectra1DImageFilter, ImageToImageFilter );
  itkNewMacro( Self );

  /** Set/get the input image containning the support window for local spectra
   * computation. */
  itkSetInputMacro( SupportWindowImage, SupportWindowImageType );
  itkGetInputMacro( SupportWindowImage, SupportWindowImageType );

protected:
  Spectra1DImageFilter();
  virtual ~Spectra1DImageFilter() {};

  typedef typename OutputImageType::RegionType OutputImageRegionType;

  virtual void BeforeThreadedGenerateData() ITK_OVERRIDE;
  virtual void ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId ) ITK_OVERRIDE;

private:
  Spectra1DImageFilter( const Self & ); // purposely not implemented
  void operator=( const Self & ); // purposely not implemented

  typedef vcl_complex< ScalarType >                  ComplexType;
  typedef vnl_vector< ComplexType >                  ComplexVectorType;
  typedef vnl_vector< ScalarType >                   SpectraVectorType;
  typedef typename InputImageType::IndexType         IndexType;
  typedef std::pair< IndexType, SpectraVectorType >  SpectraLineType;
  typedef std::deque< SpectraLineType >              SpectraLinesContainerType;
  typedef typename SupportWindowImageType::PixelType SupportWindowType;
  typedef ImageRegionConstIterator< InputImageType > InputImageIteratorType;
  typedef vnl_fft_1d< ScalarType >                   FFT1DType;

  typedef Spectra1DSupportWindowImageFilter< OutputImageType >     Spectra1DSupportWindowFilterType;
  typedef typename Spectra1DSupportWindowFilterType::FFT1DSizeType FFT1DSizeType;

  struct PerThreadData
    {
    ComplexVectorType                 ComplexVector;
    SpectraVectorType                 SpectraVector;
    typename InputImageType::SizeType LineImageRegionSize;
    };
  typedef std::vector< PerThreadData > PerThreadDataContainerType;
  PerThreadDataContainerType m_PerThreadDataContainer;

  SpectraLineType ComputeSpectra( const IndexType & lineIndex, ThreadIdType threadId );
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSpectra1DImageFilter.hxx"
#endif

#endif // __itkSpectra1DImageFilter_h