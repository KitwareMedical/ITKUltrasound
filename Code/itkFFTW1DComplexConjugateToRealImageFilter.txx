#ifndef __itkFFTW1DComplexConjugateToRealImageFilter_txx
#define __itkFFTW1DComplexConjugateToRealImageFilter_txx

#include "itkFFT1DComplexConjugateToRealImageFilter.txx"
#include "itkFFTW1DComplexConjugateToRealImageFilter.h"

#include "itkFFTWCommon.h"
#include "itkIndent.h"
#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkMetaDataObject.h"

namespace itk
{

template <typename TPixel, unsigned int Dimension>
void
FFTW1DComplexConjugateToRealImageFilter<TPixel,Dimension>::
GenerateData()
{
  // get pointers to the input and output
  typename TInputImageType::ConstPointer  inputPtr  = this->GetInput();
  typename TOutputImageType::Pointer      outputPtr = this->GetOutput();

  if ( !inputPtr || !outputPtr )
    {
    return;
    }

  // allocate output buffer memory
  outputPtr->SetBufferedRegion( outputPtr->GetRequestedRegion() );
  outputPtr->Allocate();

  const typename TInputImageType::SizeType&   inputSize
    = inputPtr->GetRequestedRegion().GetSize();
  const typename TOutputImageType::SizeType&   outputSize
    = outputPtr->GetRequestedRegion().GetSize();

  // figure out sizes
  // size of input and output aren't the same which is handled in the superclass,
  // sort of.
  // the input size and output size only differ in the fastest moving dimension
  unsigned int total_inputSize = inputSize[this->m_Direction];
  unsigned int total_outputSize = outputSize[this->m_Direction];

  if(this->m_PlanComputed)            // if we've already computed a plan
    {
    // if the image sizes aren't the same,
    // we have to compute the plan again
    if(this->m_LastImageSize != total_inputSize)
      {
      delete [] this->m_InputBuffer;
      delete [] this->m_OutputBuffer;
      FFTW1DProxyType::DestroyPlan(this->m_Plan);
      this->m_PlanComputed = false;
      }
    }
  if(!this->m_PlanComputed)
    {
    try
      {
      this->m_InputBuffer =
        new typename FFTW1DProxyType::ComplexType[total_inputSize];
      this->m_OutputBuffer =
        new typename FFTW1DProxyType::ComplexType[total_inputSize];
      }
    catch( std::bad_alloc & )
      {
      itkExceptionMacro("Problem allocating memory for internal computations");
      }
    this->m_LastImageSize = total_inputSize;
    this->m_Plan = FFTW1DProxyType::Plan_dft_1d(outputSize[this->m_Direction],
                                         this->m_InputBuffer,
                                         this->m_OutputBuffer,
					 FFTW_BACKWARD,
                                         FFTW_ESTIMATE);
    this->m_PlanComputed = true;
    }

  typedef itk::ImageLinearConstIteratorWithIndex< TInputImageType >  InputIteratorType;
  typedef itk::ImageLinearIteratorWithIndex< TOutputImageType >      OutputIteratorType;
  InputIteratorType inputIt( inputPtr, inputPtr->GetRequestedRegion() );
  // the output region should be the same as the input region in the non-fft directions
  OutputIteratorType outputIt( outputPtr, outputPtr->GetRequestedRegion() );

  inputIt.SetDirection(this->m_Direction);
  outputIt.SetDirection(this->m_Direction);

  typename InputIteratorType::PixelType* inputBufferIt;
  typename FFTW1DProxyType::ComplexType* outputBufferIt;

  // for every fft line
  for( inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd();
    outputIt.NextLine(), inputIt.NextLine() )
    {
    // copy the input line into our buffer
    inputIt.GoToBeginOfLine();
    inputBufferIt = reinterpret_cast< typename InputIteratorType::PixelType* >( this->m_InputBuffer );
    while( !inputIt.IsAtEndOfLine() )
      {
      *inputBufferIt = inputIt.Get();
      ++inputIt;
      ++inputBufferIt;
      }

    // do the transform
    FFTW1DProxyType::Execute(this->m_Plan);

    // copy the output from the buffer into our line
    outputBufferIt = this->m_OutputBuffer;
    outputIt.GoToBeginOfLine();
    while( !outputIt.IsAtEndOfLine() )
      {
      outputIt.Set( (*outputBufferIt)[0] / total_outputSize);
      ++outputIt;
      ++outputBufferIt;
      }
    }
}

} // namespace itk

#endif //_itkFFTW1DComplexConjugateToRealImageFilter_txx
