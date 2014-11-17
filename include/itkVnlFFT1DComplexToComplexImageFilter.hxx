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
#ifndef __itkVnlFFT1DComplexToComplexImageFilter_hxx
#define __itkVnlFFT1DComplexToComplexImageFilter_hxx

#include "itkVnlFFT1DComplexToComplexImageFilter.h"

#include "itkFFT1DComplexToComplexImageFilter.hxx"
#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkIndent.h"
#include "itkMetaDataObject.h"
#include "itkExceptionObject.h"
#include "vnl/algo/vnl_fft_base.h"
#include "vnl/algo/vnl_fft_1d.h"

namespace itk
{

template <class TPixel, unsigned int VDimension>
void
VnlFFT1DComplexToComplexImageFilter<TPixel,VDimension>
::ThreadedGenerateData( const OutputImageRegionType& outputRegion, ThreadIdType itkNotUsed( threadID ) )
{
  // get pointers to the input and output
  typename Superclass::InputImageType::ConstPointer  inputPtr  = this->GetInput();
  typename Superclass::OutputImageType::Pointer      outputPtr = this->GetOutput();

  if ( !inputPtr || !outputPtr )
    {
    return;
    }

  const typename Superclass::InputImageType::SizeType&   inputSize
    = inputPtr->GetRequestedRegion().GetSize();

  unsigned int vec_size = inputSize[this->m_Direction];

  typedef itk::ImageLinearConstIteratorWithIndex< InputImageType >  InputIteratorType;
  typedef itk::ImageLinearIteratorWithIndex< OutputImageType >      OutputIteratorType;
  InputIteratorType inputIt( inputPtr, outputRegion );
  OutputIteratorType outputIt( outputPtr, outputRegion );

  inputIt.SetDirection(this->m_Direction);
  outputIt.SetDirection(this->m_Direction);

  vnl_vector< vcl_complex<TPixel> > inputBuffer( vec_size );
  typename vnl_vector< vcl_complex< TPixel > >::iterator inputBufferIt  = inputBuffer.begin();
    // fft is done in-place
  typename vnl_vector< vcl_complex< TPixel > >::iterator outputBufferIt = inputBuffer.begin();
  vnl_fft_1d<TPixel> v1d(vec_size);

  // for every fft line
  for( inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd();
    outputIt.NextLine(), inputIt.NextLine() )
    {
    // copy the input line into our buffer
    inputIt.GoToBeginOfLine();
    inputBufferIt = inputBuffer.begin();
    while( !inputIt.IsAtEndOfLine() )
      {
      *inputBufferIt = inputIt.Get();
      ++inputIt;
      ++inputBufferIt;
      }

    // do the transform
    if( this->m_TransformDirection == Superclass::DIRECT )
      {
      v1d.bwd_transform(inputBuffer);
      // copy the output from the buffer into our line
      outputBufferIt = inputBuffer.begin();
      outputIt.GoToBeginOfLine();
      while( !outputIt.IsAtEndOfLine() )
	{
	outputIt.Set( *outputBufferIt );
	++outputIt;
	++outputBufferIt;
	}
      }
    else // m_TransformDirection == INVERSE
      {
      v1d.fwd_transform(inputBuffer);
      // copy the output from the buffer into our line
      outputBufferIt = inputBuffer.begin();
      outputIt.GoToBeginOfLine();
      while( !outputIt.IsAtEndOfLine() )
	{
	outputIt.Set( (*outputBufferIt) / static_cast< TPixel >( vec_size ));
	++outputIt;
	++outputBufferIt;
	}
      }
    }
}

} // end namespace itk

#endif
