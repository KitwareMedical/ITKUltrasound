#ifndef __itkFFTW1DRealToComplexConjugateImageFilter_txx
#define __itkFFTW1DRealToComplexConjugateImageFilter_txx

#include "itkFFT1DRealToComplexConjugateImageFilter.txx"
#include "itkFFTW1DRealToComplexConjugateImageFilter.h"

#include "itkFFTWCommonExtended.h"
#include "itkIndent.h"
#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkMetaDataObject.h"

namespace itk
{

template <typename TPixel, unsigned int Dimension>
void
FFTW1DRealToComplexConjugateImageFilter<TPixel,Dimension>::
DestroyPlans()
{
  for( unsigned int i = 0; i < m_PlanArray.size(); i++ )
    {
    FFTW1DProxyType::DestroyPlan( m_PlanArray[i]  );
    delete[] m_InputBufferArray[i];
    delete[] m_OutputBufferArray[i];
    this->m_PlanComputed = false;
    }
}

template <typename TPixel, unsigned int Dimension>
void
FFTW1DRealToComplexConjugateImageFilter<TPixel,Dimension>::
BeforeThreadedGenerateData()
{
  typename OutputImageType::Pointer      outputPtr = this->GetOutput();
  if ( !outputPtr )
    {
    return;
    }
  
  const typename OutputImageType::SizeType& outputSize = outputPtr->GetRequestedRegion().GetSize();
  const unsigned int lineSize = outputSize[this->m_Direction];

  if( this->m_PlanComputed )
    {
    // if the image sizes aren't the same,
    // we have to comput the plan again
    if( this->m_LastImageSize != lineSize )
      {
      this->DestroyPlans();
      }
    }
  if( ! this->m_PlanComputed )
    {
    const int threads = this->GetNumberOfThreads();
    m_PlanArray.resize( threads );
    m_InputBufferArray.resize( threads );
    m_OutputBufferArray.resize( threads );
    for( int i = 0; i < threads; i++ )
      {
      try
	{
	m_InputBufferArray[i]  = new typename FFTW1DProxyType::ComplexType[lineSize];
	m_OutputBufferArray[i] = new typename FFTW1DProxyType::ComplexType[lineSize];
	}
      catch( std::bad_alloc & )
	{
	itkExceptionMacro("Problem allocating memory for internal computations");
	}
      m_PlanArray[i] = FFTW1DProxyType::Plan_dft_1d( lineSize,
				       m_InputBufferArray[i],
				       m_OutputBufferArray[i],
				       FFTW_FORWARD,
				       FFTW_ESTIMATE,
				       threads );
      }
    this->m_LastImageSize = lineSize;
    this->m_PlanComputed = true;
    }
}

template <typename TPixel, unsigned int Dimension>
void
FFTW1DRealToComplexConjugateImageFilter<TPixel,Dimension>::
ThreadedGenerateData( const OutputImageRegionType& outputRegion, int threadID )
{
  // get pointers to the input and output
  typename InputImageType::ConstPointer  inputPtr  = this->GetInput();
  typename OutputImageType::Pointer      outputPtr = this->GetOutput();

  if ( !inputPtr || !outputPtr )
    {
    return;
    }

  typedef itk::ImageLinearConstIteratorWithIndex< InputImageType >  InputIteratorType;
  typedef itk::ImageLinearIteratorWithIndex< OutputImageType >      OutputIteratorType;
  InputIteratorType inputIt( inputPtr, outputRegion );
  OutputIteratorType outputIt( outputPtr, outputRegion );

  inputIt.SetDirection(this->m_Direction);
  outputIt.SetDirection(this->m_Direction);

  typename FFTW1DProxyType::ComplexType* inputBufferIt;
  typename FFTW1DProxyType::ComplexType* outputBufferIt;

  // for every fft line
  for( inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd();
    outputIt.NextLine(), inputIt.NextLine() )
    {
    // copy the input line into our buffer
    inputIt.GoToBeginOfLine();
    inputBufferIt = m_InputBufferArray[threadID];
    while( !inputIt.IsAtEndOfLine() )
      {
      (*inputBufferIt)[0] = inputIt.Get();
      (*inputBufferIt)[1] = 0.;
      ++inputIt;
      ++inputBufferIt;
      }

    // do the transform
    FFTW1DProxyType::Execute( m_PlanArray[threadID] );

    // copy the output from the buffer into our line
    outputBufferIt = m_OutputBufferArray[threadID];
    outputIt.GoToBeginOfLine();
    while( !outputIt.IsAtEndOfLine() )
      {
      outputIt.Set( *(reinterpret_cast<typename OutputIteratorType::PixelType*>(outputBufferIt)) );
      ++outputIt;
      ++outputBufferIt;
      }
    }
}

} // namespace itk

#endif //_itkFFTW1DRealToComplexConjugateImageFilter_txx
