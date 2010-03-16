#ifndef __itkFFTW1DComplexToComplexImageFilter_txx
#define __itkFFTW1DComplexToComplexImageFilter_txx

#include "itkFFT1DComplexToComplexImageFilter.txx"
#include "itkFFTW1DComplexToComplexImageFilter.h"

#include "itkFFTWCommonExtended.h"
#include "itkIndent.h"
#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkMetaDataObject.h"

namespace itk
{

template <typename TPixel, unsigned int Dimension>
void
FFTW1DComplexToComplexImageFilter<TPixel,Dimension>::
GenerateData()
{
  // get pointers to the input and output
  typename InputImageType::ConstPointer  inputPtr  = this->GetInput();
  typename OutputImageType::Pointer      outputPtr = this->GetOutput();

  if ( !inputPtr || !outputPtr )
    {
    return;
    }

  // allocate output buffer memory
  outputPtr->SetBufferedRegion( outputPtr->GetRequestedRegion() );
  outputPtr->Allocate();

  const typename InputImageType::SizeType&   inputSize
    = inputPtr->GetRequestedRegion().GetSize();
  const typename OutputImageType::SizeType&   outputSize
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
        new typename FFTW1DProxyType::ComplexType[total_outputSize];
      }
    catch( std::bad_alloc & )
      {
      itkExceptionMacro("Problem allocating memory for internal computations");
      }
    this->m_LastImageSize = total_inputSize;
    int sign;
    if( this->m_TransformDirection == Superclass::DIRECT )
      sign = FFTW_FORWARD;
    else
      sign = FFTW_BACKWARD;
    this->m_Plan = FFTW1DProxyType::Plan_dft_1d(outputSize[this->m_Direction],
                                         this->m_InputBuffer,
                                         this->m_OutputBuffer,
					 sign,
                                         FFTW_ESTIMATE);
    this->m_PlanComputed = true;
    }

  typedef itk::ImageLinearConstIteratorWithIndex< InputImageType >  InputIteratorType;
  typedef itk::ImageLinearIteratorWithIndex< OutputImageType >      OutputIteratorType;
  InputIteratorType inputIt( inputPtr, inputPtr->GetRequestedRegion() );
  // the output region should be the same as the input region in the non-fft directions
  OutputIteratorType outputIt( outputPtr, outputPtr->GetRequestedRegion() );

  inputIt.SetDirection(this->m_Direction);
  outputIt.SetDirection(this->m_Direction);

  typename InputIteratorType::PixelType* inputBufferIt;
  typename OutputIteratorType::PixelType* outputBufferIt;

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

    if( this->m_TransformDirection == Superclass::DIRECT )
      {
      // copy the output from the buffer into our line
      outputBufferIt = reinterpret_cast< typename OutputIteratorType::PixelType* >( this->m_OutputBuffer );
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
      // copy the output from the buffer into our line
      outputBufferIt = reinterpret_cast< typename OutputIteratorType::PixelType* >( this->m_OutputBuffer );
      outputIt.GoToBeginOfLine();
      while( !outputIt.IsAtEndOfLine() )
	{
	outputIt.Set( *outputBufferIt / static_cast< TPixel >( total_outputSize ));
	++outputIt;
	++outputBufferIt;
	}
      }
    }
}

} // namespace itk

#endif //_itkFFTW1DComplexToComplexImageFilter_txx
