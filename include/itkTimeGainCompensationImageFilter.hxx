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
::TimeGainCompensationImageFilter():
  m_Gain( 2, 2 )
{
  m_Gain(0, 0) = NumericTraits< double >::min();
  m_Gain(0, 1) = NumericTraits< double >::OneValue();
  m_Gain(1, 0) = NumericTraits< double >::max();
  m_Gain(1, 1) = NumericTraits< double >::OneValue();
}


template< typename TInputImage, typename TOutputImage >
void
TimeGainCompensationImageFilter< TInputImage, TOutputImage >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Gain:" << std::endl;
  for( unsigned int ii = 0; ii < m_Gain.rows(); ++ii )
    {
    os << indent.GetNextIndent() << "[" << m_Gain( ii, 0 ) << ", " << m_Gain( ii, 1 ) << "]" << std::endl;
    }
}


template< typename TInputImage, typename TOutputImage >
void
TimeGainCompensationImageFilter< TInputImage, TOutputImage >
::BeforeThreadedGenerateData()
{
  const GainType & gain = this->GetGain();
  if( gain.cols() != 2 )
    {
    itkExceptionMacro( "Gain should have two columns." );
    }
  if( gain.rows() < 2 )
    {
    itkExceptionMacro( "Insufficient depths specified in Gain." );
    }
  double depth = gain( 0, 0 );
  for( unsigned int ii = 1; ii < gain.rows(); ++ii )
    {
    if( gain( ii, 0 ) <= depth )
      {
      itkExceptionMacro( "Gain depths must be strictly increasing." );
      }
    depth = gain( ii, 0 );
    }
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
