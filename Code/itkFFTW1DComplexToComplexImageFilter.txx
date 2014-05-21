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
FFTW1DComplexToComplexImageFilter<TPixel,Dimension>::
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
      if( this->m_TransformDirection == Superclass::DIRECT )
	{
	m_PlanArray[i] = FFTW1DProxyType::Plan_dft_1d( lineSize,
					 m_InputBufferArray[i],
					 m_OutputBufferArray[i],
					 FFTW_FORWARD,
					 FFTW_ESTIMATE,
					 threads );
	}
      else // m_TransformDirection == INVERSE
	{
	m_PlanArray[i] = FFTW1DProxyType::Plan_dft_1d( lineSize,
					 m_InputBufferArray[i],
					 m_OutputBufferArray[i],
					 FFTW_BACKWARD,
					 FFTW_ESTIMATE,
					 threads );
	}
      }
    this->m_LastImageSize = lineSize;
    this->m_PlanComputed = true;
    }
}

template <typename TPixel, unsigned int Dimension>
void
FFTW1DComplexToComplexImageFilter<TPixel,Dimension>::
ThreadedGenerateData( const OutputImageRegionType& outputRegion, ThreadIdType threadID )
{
  // get pointers to the input and output
  typename InputImageType::ConstPointer  inputPtr  = this->GetInput();
  typename OutputImageType::Pointer      outputPtr = this->GetOutput();

  if ( !inputPtr || !outputPtr )
    {
    return;
    }

  const typename OutputImageType::SizeType& outputSize = outputPtr->GetRequestedRegion().GetSize();
  const unsigned int lineSize = outputSize[this->m_Direction];

  typedef itk::ImageLinearConstIteratorWithIndex< InputImageType >  InputIteratorType;
  typedef itk::ImageLinearIteratorWithIndex< OutputImageType >      OutputIteratorType;
  InputIteratorType inputIt( inputPtr, outputRegion );
  // the output region should be the same as the input region in the non-fft directions
  OutputIteratorType outputIt( outputPtr, outputRegion );

  inputIt.SetDirection(this->m_Direction);
  outputIt.SetDirection(this->m_Direction);

  typename InputIteratorType::PixelType*  inputBufferIt;
  typename OutputIteratorType::PixelType* outputBufferIt;

  // for every fft line
  for( inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd();
    outputIt.NextLine(), inputIt.NextLine() )
    {
    // copy the input line into our buffer
    inputIt.GoToBeginOfLine();
    inputBufferIt = reinterpret_cast< typename InputIteratorType::PixelType* >( m_InputBufferArray[threadID] );
    while( !inputIt.IsAtEndOfLine() )
      {
      *inputBufferIt = inputIt.Get();
      ++inputIt;
      ++inputBufferIt;
      }

    // do the transform
    FFTW1DProxyType::Execute( m_PlanArray[threadID] );

    if( this->m_TransformDirection == Superclass::DIRECT )
      {
      // copy the output from the buffer into our line
      outputBufferIt = reinterpret_cast< typename OutputIteratorType::PixelType* >( m_OutputBufferArray[threadID] );
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
      outputBufferIt = reinterpret_cast< typename OutputIteratorType::PixelType* >( m_OutputBufferArray[threadID] );
      outputIt.GoToBeginOfLine();
      while( !outputIt.IsAtEndOfLine() )
	{
	outputIt.Set( *outputBufferIt / static_cast< TPixel >( lineSize ));
	++outputIt;
	++outputBufferIt;
	}
      }
    }
}

} // namespace itk

#endif //_itkFFTW1DComplexToComplexImageFilter_txx
