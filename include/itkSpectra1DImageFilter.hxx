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
#ifndef __itkSpectra1DImageFilter_hxx
#define __itkSpectra1DImageFilter_hxx

#include "itkSpectra1DImageFilter.h"

#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkMetaDataObject.h"

#include "vnl/algo/vnl_fft_base.h"
#include "vnl/algo/vnl_fft_1d.h"

#include "itkSpectra1DSupportWindowImageFilter.h"

#include <utility>

namespace itk
{

template< typename TInputImage, typename TSupportWindowImage, typename TOutputImage >
Spectra1DImageFilter< TInputImage, TSupportWindowImage, TOutputImage >
::Spectra1DImageFilter()
{
  this->AddRequiredInputName( "SupportWindowImage" );
}


template< typename TInputImage, typename TSupportWindowImage, typename TOutputImage >
void
Spectra1DImageFilter< TInputImage, TSupportWindowImage, TOutputImage >
::ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread, ThreadIdType itkNotUsed( threadId ) )
{
  OutputImageType * output = this->GetOutput();
  const InputImageType * input = this->GetInput();
  const SupportWindowImageType * supportWindowImage = this->GetSupportWindowImage();

  typedef ImageLinearIteratorWithIndex< OutputImageType > OutputIteratorType;
  OutputIteratorType outputIt( output, outputRegionForThread );
  outputIt.SetDirection( 1 );

  typedef Spectra1DSupportWindowImageFilter< OutputImageType > Spectra1DSupportWindowFilterType;
  typedef typename Spectra1DSupportWindowFilterType::FFT1DSizeType FFT1DSizeType;
  const MetaDataDictionary & dict = supportWindowImage->GetMetaDataDictionary();
  FFT1DSizeType fft1DSize = 32;
  ExposeMetaData< FFT1DSizeType >( dict, "FFT1DSize", fft1DSize );

  typedef vcl_complex< ScalarType >                  ComplexType;
  typedef vnl_vector< ComplexType >                  ComplexVectorType;
  typedef vnl_vector< ScalarType >                   SpectraVectorType;
  typedef typename InputImageType::IndexType         IndexType;
  typedef std::pair< IndexType, SpectraVectorType >  SpectraLineType;
  typedef std::deque< SpectraLineType >              SpectraLinesContainerType;
  typedef typename SupportWindowImageType::PixelType SupportWindowType;
  typedef ImageRegionConstIterator< InputImageType > InputImageIteratorType;

  ComplexVectorType complexVector( fft1DSize );
  SpectraVectorType spectraVector( fft1DSize );
  SpectraLinesContainerType spectraLines;
  typename InputImageType::SizeType lineImageRegionSize;
  lineImageRegionSize.Fill( 1 );
  lineImageRegionSize[0] = fft1DSize;
  vnl_fft_1d< ScalarType > fft1D( fft1DSize );

  typedef ImageLinearConstIteratorWithIndex< SupportWindowImageType > SupportWindowIteratorType;
  SupportWindowIteratorType supportWindowIt( supportWindowImage, outputRegionForThread );
  supportWindowIt.SetDirection( 1 );


  for( outputIt.GoToBegin(), supportWindowIt.GoToBegin();
       !outputIt.IsAtEnd();
       outputIt.NextLine(), supportWindowIt.NextLine() )
    {
    spectraLines.clear();
    while( ! outputIt.IsAtEndOfLine() )
      {
      const SupportWindowType & supportWindow = supportWindowIt.Value();
      if( spectraLines.size() == 0 ) // first window in this lateral direction
        {
        const typename SupportWindowType::const_iterator windowLineEnd = supportWindow.end();
        for( typename SupportWindowType::const_iterator windowLine = supportWindow.begin();
             windowLine != windowLineEnd;
             ++windowLine )
          {
          const IndexType & lineIndex = *windowLine;
          const typename InputImageType::RegionType lineRegion( lineIndex, lineImageRegionSize );
          InputImageIteratorType inputIt( input, lineRegion );
          inputIt.GoToBegin();
          complexVector.fill( 0 );
          typename ComplexVectorType::iterator complexVectorIt = complexVector.begin();
          while( !inputIt.IsAtEnd() )
            {
            *complexVectorIt = inputIt.Value();
            ++inputIt;
            ++complexVectorIt;
            }
          fft1D.bwd_transform( complexVector );
          typename ComplexVectorType::const_iterator complexVectorConstIt = complexVector.begin();
          typename SpectraVectorType::iterator spectraVectorIt = spectraVector.begin();
          const typename SpectraVectorType::iterator spectraVectorItEnd = spectraVector.end();
          while( spectraVectorIt != spectraVectorItEnd )
            {
            *spectraVectorIt = std::real(*complexVectorConstIt * std::conj(*complexVectorConstIt));
            ++spectraVectorIt;
            ++complexVectorConstIt;
            }
          const SpectraLineType spectraLine = std::make_pair( lineIndex, spectraVector );
          spectraLines.push_back( spectraLine );
          }
        }
      else
        {
        // todo
        }
      ++outputIt;
      ++supportWindowIt;
      }
    }
}


} // end namespace itk

#endif
