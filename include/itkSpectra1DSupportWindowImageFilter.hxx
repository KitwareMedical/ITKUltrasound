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
#ifndef __itkSpectra1DSupportWindowImageFilter_hxx
#define __itkSpectra1DSupportWindowImageFilter_hxx

#include "itkSpectra1DSupportWindowImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkMetaDataObject.h"

namespace itk
{

template< typename TInputImage >
Spectra1DSupportWindowImageFilter< TInputImage >
::Spectra1DSupportWindowImageFilter():
  m_FFT1DSize( 32 )
{
}


template< typename TInputImage >
void
Spectra1DSupportWindowImageFilter< TInputImage >
::ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread, ThreadIdType itkNotUsed( threadId ) )
{
  OutputImageType * output = this->GetOutput();
  const InputImageType * input = this->GetInput();

  const OutputImageRegionType outputLargestRegion = output->GetLargestPossibleRegion();
  typedef typename OutputImageType::IndexType IndexType;
  const IndexType largestIndexStart = outputLargestRegion.GetIndex();
  IndexType largestIndexStop = largestIndexStart + outputLargestRegion.GetSize();
  for( unsigned int dim = 0; dim < ImageDimension; ++dim )
    {
    largestIndexStop[dim] -= 1;
    }

  typedef ImageRegionConstIteratorWithIndex< InputImageType > InputIteratorType;
  InputIteratorType inputIt( input, outputRegionForThread );
  typedef ImageRegionIterator< OutputImageType > OutputIteratorType;
  OutputIteratorType outputIt( output, outputRegionForThread );
  const FFT1DSizeType fftSize = this->GetFFT1DSize();
  if( outputLargestRegion.GetSize()[0] < fftSize )
    {
    itkExceptionMacro( "Insufficient size in the FFT direction." );
    }
  for( inputIt.GoToBegin(), outputIt.GoToBegin(); !outputIt.IsAtEnd(); ++inputIt, ++outputIt )
    {
    OutputPixelType & supportWindow = outputIt.Value();
    supportWindow.clear();

    const IndexType inputIndex = inputIt.GetIndex();

    IndexType lineIndex;
    lineIndex[0] = inputIndex[0] - fftSize / 2;
    if( lineIndex[0] < largestIndexStart[0] )
      {
      lineIndex[0] = largestIndexStart[0];
      }

    if( lineIndex[0] + fftSize > largestIndexStop[0] )
      {
      lineIndex[0] = largestIndexStop[0] - fftSize;
      }

    const IndexValueType sideLines = static_cast< IndexValueType >( inputIt.Get() );
    for( IndexValueType line = inputIndex[1] - sideLines;
         line < inputIndex[1] + sideLines;
         ++line )
      {
      if( line < largestIndexStart[1] || line > largestIndexStop[1] )
        {
        continue;
        }
      lineIndex[1] = line;
      supportWindow.push_back( lineIndex );
      }
    }
}


template< typename TInputImage >
void
Spectra1DSupportWindowImageFilter< TInputImage >
::AfterThreadedGenerateData()
{
  OutputImageType * output = this->GetOutput();
  MetaDataDictionary & dict = output->GetMetaDataDictionary();
  EncapsulateMetaData< FFT1DSizeType >( dict, "FFT1DSize", this->GetFFT1DSize() );
}

} // end namespace itk

#endif
