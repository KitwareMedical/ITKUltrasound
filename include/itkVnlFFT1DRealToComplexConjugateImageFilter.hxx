/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: itkVnlFFT1DRealToComplexConjugateImageFilter.hxx,v $
Language:  C++
Date:      $Date: 2009-01-27 19:30:16 $
Version:   $Revision: 1.12 $

Copyright (c) 2002 Insight Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVnlFFT1DRealToComplexConjugateImageFilter_hxx
#define __itkVnlFFT1DRealToComplexConjugateImageFilter_hxx

#include "itkVnlFFT1DRealToComplexConjugateImageFilter.h"

#include "itkFFT1DRealToComplexConjugateImageFilter.hxx"
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
bool VnlFFT1DRealToComplexConjugateImageFilter<TPixel,VDimension>
::Legaldim(int n)
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
VnlFFT1DRealToComplexConjugateImageFilter<TPixel,VDimension>
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

  std::cout << "***direction: " << this->m_Direction << std::endl;
  unsigned int vecSize = inputSize[this->m_Direction];
  std::cout << "***vecSize: " << vecSize << std::endl;
  if( !this->Legaldim(vecSize) )
    {
    ExceptionObject exception(__FILE__, __LINE__);
    exception.SetDescription("Illegal Array DIM for FFT");
    exception.SetLocation(ITK_LOCATION);
    throw exception;
    }


  typedef itk::ImageLinearConstIteratorWithIndex< InputImageType >  InputIteratorType;
  typedef itk::ImageLinearIteratorWithIndex< OutputImageType >      OutputIteratorType;
  InputIteratorType inputIt( inputPtr, outputRegion );
  OutputIteratorType outputIt( outputPtr, outputRegion );

  inputIt.SetDirection(this->m_Direction);
  outputIt.SetDirection(this->m_Direction);

  vnl_vector< vcl_complex<TPixel> > inputBuffer( vecSize );
  typename vnl_vector< vcl_complex< TPixel > >::iterator inputBufferIt = inputBuffer.begin();
    // fft is done in-place
  typename vnl_vector< vcl_complex< TPixel > >::iterator outputBufferIt = inputBuffer.begin();
  vnl_fft_1d<TPixel> v1d(vecSize);

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
