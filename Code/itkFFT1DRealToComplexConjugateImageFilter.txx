#ifndef __itkFFT1DRealToComplexConjugateImageFilter_txx
#define __itkFFT1DRealToComplexConjugateImageFilter_txx

#include "itkFFT1DRealToComplexConjugateImageFilter.h"

#include "itkVnlFFT1DRealToComplexConjugateImageFilter.h"

#if defined(USE_OPENCL_FFT)
#include "itkOpenCL1DRealToComplexConjugateImageFilter.h"
#endif
#if defined(USE_FFTWD) || defined(USE_FFTWF)
#include "itkFFTW1DRealToComplexConjugateImageFilter.h"
#endif

#include "itkMetaDataObject.h"

namespace itk
{

template < class TPixel, unsigned int VDimension >
typename FFT1DRealToComplexConjugateImageFilter< TPixel, VDimension >::Pointer
FFT1DRealToComplexConjugateImageFilter< TPixel, VDimension >
::New()
{
  Pointer smartPtr = ::itk::ObjectFactory< Self >::Create();

#ifdef USE_OPENCL_FFT
  if( smartPtr.IsNull() )
    {
    if( typeid( TPixel ) == typeid( float ) )
      {
      smartPtr = dynamic_cast<Self *>(
	OpenCL1DRealToComplexConjugateImageFilter< float, VDimension >
	::New().GetPointer() );
      }
    }
#endif
#ifdef USE_FFTWD
  if( smartPtr.IsNull() )
    {
    if( typeid( TPixel ) == typeid( double ) )
      {
      smartPtr = dynamic_cast< Self* >(
	FFTW1DRealToComplexConjugateImageFilter< double, VDimension >
	::New().GetPointer() );
      }
    }
#endif
#ifdef USE_FFTWF
  if( smartPtr.IsNull() )
    {
    if( typeid( TPixel ) == typeid( float ) )
      {
      smartPtr = dynamic_cast<Self *>(
	FFTW1DRealToComplexConjugateImageFilter< float, VDimension >
	::New().GetPointer() );
      }
    }
#endif

  if( smartPtr.IsNull() )
    {
    smartPtr = VnlFFT1DRealToComplexConjugateImageFilter< TPixel, VDimension >
      ::New().GetPointer();
    }

  return smartPtr;
}

template < class TPixel , unsigned int Dimension >
void
FFT1DRealToComplexConjugateImageFilter < TPixel , Dimension >
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // get pointers to the inputs
  typename InputImageType::Pointer inputPtr  =
    const_cast<InputImageType *> (this->GetInput());
  typename OutputImageType::Pointer outputPtr = this->GetOutput();

  if ( !inputPtr || !outputPtr )
    {
    return;
    }

  // we need to compute the input requested region (size and start index)
  typedef const typename OutputImageType::SizeType& OutputSizeType;
  OutputSizeType outputRequestedRegionSize =
    outputPtr->GetRequestedRegion().GetSize();
  typedef const typename OutputImageType::IndexType& OutputIndexType;
  OutputIndexType outputRequestedRegionStartIndex =
    outputPtr->GetRequestedRegion().GetIndex();

  //// the regions other than the fft direction are fine
  typename InputImageType::SizeType  inputRequestedRegionSize = outputRequestedRegionSize;
  typename InputImageType::IndexType inputRequestedRegionStartIndex = outputRequestedRegionStartIndex;

  // we but need all of the input in the fft direction
  const unsigned int direction = this->m_Direction;
  const typename InputImageType::SizeType& inputLargeSize =
    inputPtr->GetLargestPossibleRegion().GetSize();
  inputRequestedRegionSize[direction] = inputLargeSize[direction];
  const typename InputImageType::IndexType& inputLargeIndex =
    inputPtr->GetLargestPossibleRegion().GetIndex();
  inputRequestedRegionStartIndex[direction] = inputLargeIndex[direction];

  typename InputImageType::RegionType inputRequestedRegion;
  inputRequestedRegion.SetSize( inputRequestedRegionSize );
  inputRequestedRegion.SetIndex( inputRequestedRegionStartIndex );

  inputPtr->SetRequestedRegion( inputRequestedRegion );
}


template < class TPixel , unsigned int Dimension >
void
FFT1DRealToComplexConjugateImageFilter < TPixel , Dimension >
::EnlargeOutputRequestedRegion(DataObject *output)
{
  OutputImageType* outputPtr = dynamic_cast<OutputImageType*>( output );

  // we need to enlarge the region in the fft direction to the
  // largest possible in that direction
  typedef const typename OutputImageType::SizeType& ConstOutputSizeType;
  ConstOutputSizeType requestedSize =
    outputPtr->GetRequestedRegion().GetSize();
  ConstOutputSizeType outputLargeSize =
    outputPtr->GetLargestPossibleRegion().GetSize();
  typedef const typename OutputImageType::IndexType& ConstOutputIndexType;
  ConstOutputIndexType requestedIndex =
    outputPtr->GetRequestedRegion().GetIndex();
  ConstOutputIndexType outputLargeIndex =
    outputPtr->GetLargestPossibleRegion().GetIndex();

  typename OutputImageType::SizeType enlargedSize   = requestedSize;
  typename OutputImageType::IndexType enlargedIndex = requestedIndex;
  enlargedSize[this->m_Direction]  = outputLargeSize[this->m_Direction];
  enlargedIndex[this->m_Direction] = outputLargeIndex[this->m_Direction];

  typename OutputImageType::RegionType enlargedRegion;
  enlargedRegion.SetSize( enlargedSize );
  enlargedRegion.SetIndex( enlargedIndex );
  outputPtr->SetRequestedRegion( enlargedRegion );
}


template < class TPixel , unsigned int Dimension >
void
FFT1DRealToComplexConjugateImageFilter < TPixel , Dimension >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Direction: " << m_Direction << std::endl;
}


template < class TPixel , unsigned int Dimension >
int 
FFT1DRealToComplexConjugateImageFilter < TPixel , Dimension >
::SplitRequestedRegion(int i, int num, OutputImageRegionType& splitRegion)
{
  // Get the output pointer
  OutputImageType * outputPtr = this->GetOutput();
  const typename OutputImageType::SizeType& requestedRegionSize 
    = outputPtr->GetRequestedRegion().GetSize();

  int splitAxis;
  typename OutputImageType::IndexType splitIndex;
  typename OutputImageType::SizeType splitSize;

  // Initialize the splitRegion to the output requested region
  splitRegion = outputPtr->GetRequestedRegion();
  splitIndex = splitRegion.GetIndex();
  splitSize = splitRegion.GetSize();

  // split on the outermost dimension available
  splitAxis = outputPtr->GetImageDimension() - 1;
  while (requestedRegionSize[splitAxis] == 1)
    {
    --splitAxis;
    if (splitAxis < 0)
      { // cannot split
      itkDebugMacro("  Cannot Split");
      return 1;
      }
    }
  if( splitAxis == static_cast< int >( this->m_Direction ) )
    {
    --splitAxis;
    if (splitAxis < 0)
      { // cannot split
      itkDebugMacro("  Cannot Split");
      return 1;
      }
    }

  // determine the actual number of pieces that will be generated
  typename OutputImageType::SizeType::SizeValueType range = requestedRegionSize[splitAxis];
  int valuesPerThread = (int)::vcl_ceil(range/(double)num);
  int maxThreadIdUsed = (int)::vcl_ceil(range/(double)valuesPerThread) - 1;

  // Split the region
  if (i < maxThreadIdUsed)
    {
    splitIndex[splitAxis] += i*valuesPerThread;
    splitSize[splitAxis] = valuesPerThread;
    }
  if (i == maxThreadIdUsed)
    {
    splitIndex[splitAxis] += i*valuesPerThread;
    // last thread needs to process the "rest" dimension being split
    splitSize[splitAxis] = splitSize[splitAxis] - i*valuesPerThread;
    }
  
  // set the split region ivars
  splitRegion.SetIndex( splitIndex );
  splitRegion.SetSize( splitSize );

  itkDebugMacro("  Split Piece: " << splitRegion );

  return maxThreadIdUsed + 1;
}


} // end namespace itk

#endif // __itkFFT1DRealToComplexConjugateImageFilter_txx
