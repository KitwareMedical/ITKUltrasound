/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: itkVnlFFT1DRealToComplexConjugateImageFilter.txx,v $
Language:  C++
Date:      $Date: 2009-01-27 19:30:16 $
Version:   $Revision: 1.12 $

Copyright (c) 2002 Insight Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVnlFFT1DRealToComplexConjugateImageFilter_txx
#define __itkVnlFFT1DRealToComplexConjugateImageFilter_txx

#include "itkVnlFFT1DRealToComplexConjugateImageFilter.h"

#include "itkFFT1DRealToComplexConjugateImageFilter.txx"
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
bool VnlFFT1DRealToComplexConjugateImageFilter<TPixel,VDimension>::
Legaldim(int n)
{
  int ifac = 2;
  for (int l = 1; l <= 3; l++) 
    {
    for(; n % ifac == 0;)
      {
      n /= ifac;
      }
    ifac += l;
    }
  return (n == 1); // return false if decomposition failed
}


template <class TPixel, unsigned int VDimension>
void
VnlFFT1DRealToComplexConjugateImageFilter<TPixel,VDimension>::
ThreadedGenerateData( const OutputImageRegionType& outputRegion, int threadID )
{
  // get pointers to the input and output
  typename Superclass::TInputImageType::ConstPointer  inputPtr  = this->GetInput();
  typename Superclass::TOutputImageType::Pointer      outputPtr = this->GetOutput();

  if ( !inputPtr || !outputPtr )
    {
    return;
    }
  
  const typename Superclass::TInputImageType::SizeType&   inputSize
    = inputPtr->GetRequestedRegion().GetSize();

  unsigned int vec_size = inputSize[this->m_Direction];
  if( !this->Legaldim(vec_size) )
    {
    ExceptionObject exception(__FILE__, __LINE__);
    exception.SetDescription("Illegal Array DIM for FFT");
    exception.SetLocation(ITK_LOCATION);
    throw exception;
    }


  typedef itk::ImageLinearConstIteratorWithIndex< TInputImageType >  InputIteratorType;
  typedef itk::ImageLinearIteratorWithIndex< TOutputImageType >      OutputIteratorType;
  InputIteratorType inputIt( inputPtr, outputRegion );
  OutputIteratorType outputIt( outputPtr, outputRegion );

  inputIt.SetDirection(this->m_Direction);
  outputIt.SetDirection(this->m_Direction);

  vnl_vector< vcl_complex<TPixel> > inputBuffer( vec_size );
  typename vnl_vector< vcl_complex< TPixel > >::iterator inputBufferIt = inputBuffer.begin();
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
}

}


#endif
