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
#ifndef itkTimeGainCompensationImageFilter_hxx
#define itkTimeGainCompensationImageFilter_hxx

#include "itkTimeGainCompensationImageFilter.h"

#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"

namespace itk
{

template< typename TInputImage, typename TOutputImage >
TimeGainCompensationImageFilter< TInputImage, TOutputImage >
::TimeGainCompensationImageFilter()
{
}

template< typename TInputImage, typename TOutputImage >
void
TimeGainCompensationImageFilter< TInputImage, TOutputImage >
::ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread, ThreadIdType itkNotUsed( threadId ) )
{
  const InputImageType * inputImage = this->GetInput();
  OutputImageType * outputImage = this->GetOutput();

  typedef ImageLinearConstIteratorWithIndex< InputImageType > InputIteratorType;
  InputIteratorType inputIt( inputImage, outputRegionForThread );
  inputIt.SetDirection( 0 );
  inputIt.GoToBegin();

  typedef ImageLinearIteratorWithIndex< OutputImageType > OutputIteratorType;
  OutputIteratorType outputIt( outputImage, outputRegionForThread );
  outputIt.SetDirection( 0 );
  outputIt.GoToBegin();

  for( inputIt.GoToBegin(), outputIt.GoToBegin();
       !outputIt.IsAtEnd();
       inputIt.NextLine(), outputIt.NextLine() )
    {
    inputIt.GoToBeginOfLine();
    outputIt.GoToBeginOfLine();
    while( ! outputIt.IsAtEndOfLine() )
      {
      outputIt.Set( inputIt.Value() );
      ++inputIt;
      ++outputIt;
      }
    }
}

} // end namespace itk

#endif // itkTimeGainCompensationImageFilter_hxx
