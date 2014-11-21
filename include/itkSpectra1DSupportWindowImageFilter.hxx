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

namespace itk
{

template< typename TInputImage >
Spectra1DSupportWindowImageFilter< TInputImage >
::Spectra1DSupportWindowImageFilter():
  m_FFTSize( 32 )
{
}


template< typename TInputImage >
void
Spectra1DSupportWindowImageFilter< TInputImage >
::ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread, ThreadIdType itkNotUsed( threadId ) )
{
  OutputImageType * output = this->GetOutput();

  const OutputImageRegionType outputLargestRegion = output->GetLargestPossibleRegion();

  typedef ImageRegionIterator< OutputImageType > ImageIteratorType;
  ImageIteratorType imageIt( output, outputRegionForThread );
  for( imageIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt )
    {
    OutputPixelType & supportWindow = imageIt.Value();
    supportWindow.clear();
    }
}

} // end namespace itk

#endif
