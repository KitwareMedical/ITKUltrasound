#ifndef __itkFFT1DComplexConjugateToRealImageFilter_txx
#define __itkFFT1DComplexConjugateToRealImageFilter_txx

#include "itkFFT1DComplexConjugateToRealImageFilter.h"

#include "itkVnlFFT1DComplexConjugateToRealImageFilter.h"

//#if defined(USE_FFTWD) || defined(USE_FFTWF)
//#include "itkFFTW1DComplexConjugateToRealImageFilter.h"
//#endif

#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"

namespace itk
{

template < class TPixel, unsigned int VDimension >
typename FFT1DComplexConjugateToRealImageFilter< TPixel, VDimension >::Pointer
FFT1DComplexConjugateToRealImageFilter< TPixel, VDimension >
::New()
{
  Pointer smartPtr = ::itk::ObjectFactory< Self >::Create();

//#ifdef USE_FFTWD
  //if( smartPtr.IsNull() )
    //{
    //if( typeid( TPixel ) == typeid( double ) )
      //{
      //smartPtr = dynamic_cast< Self* >(
	//FFTW1DComplexConjugateToRealImageFilter< double, VDimension >
	//::New().GetPointer() );
      //}
    //}
//#endif
//#ifdef USE_FFTWF
  //if( smartPtr.IsNull() )
    //{
    //if( typeid( TPixel ) == typeid( float ) )
      //{
      //smartPtr = dynamic_cast<Self *>(
	//FFT1DComplexConjugateToRealImageFilter< float, VDimension >
	//::New().GetPointer() );
      //}
    //}
//#endif

  if( smartPtr.IsNull() )
    {
    smartPtr = VnlFFT1DComplexConjugateToRealImageFilter< TPixel, VDimension >
      ::New().GetPointer();
    }

  return smartPtr;
}

template < class TPixel, unsigned int VDimension >
void
FFT1DComplexConjugateToRealImageFilter< TPixel, VDimension >
::GenerateOutputInformation()
{
  // call the superclass' implementation of this method
  Superclass::GenerateOutputInformation();
  //
  // If this implementation returns a full result
  // instead of a 'half-complex' matrix, then none of this
  // is necessary
  if(this->FullMatrix())
    return;

   // get pointers to the input and output
  typename TInputImageType::ConstPointer  inputPtr  = this->GetInput();
  typename TOutputImageType::Pointer      outputPtr = this->GetOutput();

  if ( !inputPtr || !outputPtr )
    {
    return;
    }
  
  // 
  // This is all based on the same function in itk::ShrinkImageFilter
  // ShrinkImageFilter also modifies the image spacing, but spacing
  // has no meaning in the result of an FFT. For an IFFT, since the
  // spacing is propagated to the complex result, we can use the spacing
  // from the input to propagate back to the output.
  const typename TInputImageType::SizeType&   inputSize
    = inputPtr->GetLargestPossibleRegion().GetSize();
  const typename TInputImageType::IndexType&  inputStartIndex
    = inputPtr->GetLargestPossibleRegion().GetIndex();
  
  typename TOutputImageType::SizeType     outputSize = inputSize;
  typename TOutputImageType::IndexType    outputStartIndex = inputStartIndex;
  
  for (unsigned int i = 0; i < TOutputImageType::ImageDimension; i++)
    {
    outputSize[i] = inputSize[i];
    outputStartIndex[i] = inputStartIndex[i];
    }

  //
  // in 4.3.4 of the FFT documentation, they indicate the size of
  // of a real-to-complex FFT is N * N ... + (N /2+1)
  //                              1   2        d
  // complex numbers.
  // going from complex to real, you know the output is at least
  // twice the size in the last dimension as the input, but it might
  // be 2*size+1.  Consequently, the output of the FFT:R2C operation
  // 
  MetaDataDictionary &InputDic = 
    const_cast<MetaDataDictionary &>(inputPtr->GetMetaDataDictionary());
  unsigned int direction = this->m_Direction;
  outputSize[direction] = static_cast<unsigned int>(inputSize[direction])/2 + 1;
  if(this->GetActualDimensionIsOdd())
    {
    outputSize[direction]++;
    }

  // backwards compatible/deprecated version
  typedef typename TInputImageType::SizeType::SizeValueType SizeScalarType;
  SizeScalarType originalSize = 0;

  if(ExposeMetaData < SizeScalarType > 
          (InputDic,std::string("FFT_Actual_RealImage_Size"),originalSize))
    {
    outputSize[direction] = originalSize;
    }

  typename TOutputImageType::RegionType outputLargestPossibleRegion;
  outputLargestPossibleRegion.SetSize( outputSize );
  outputLargestPossibleRegion.SetIndex( outputStartIndex );
  
  outputPtr->SetLargestPossibleRegion( outputLargestPossibleRegion );
}


template < class TPixel , unsigned int Dimension >
void
FFT1DComplexConjugateToRealImageFilter < TPixel , Dimension >
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // get pointers to the inputs
  typename TInputImageType::Pointer inputPtr  =
    const_cast<TInputImageType *> (this->GetInput());
  typename TOutputImageType::Pointer outputPtr = this->GetOutput();

  if ( !inputPtr || !outputPtr )
    {
    return;
    }

  // we need to compute the input requested region (size and start index)
  typedef const typename TOutputImageType::SizeType& OutputSizeType;
  OutputSizeType outputRequestedRegionSize =
    outputPtr->GetRequestedRegion().GetSize();
  typedef const typename TOutputImageType::IndexType& OutputIndexType;
  OutputIndexType outputRequestedRegionStartIndex =
    outputPtr->GetRequestedRegion().GetIndex();

  //// the regions other than the fft direction are fine
  typename TInputImageType::SizeType  inputRequestedRegionSize = outputRequestedRegionSize;
  typename TInputImageType::IndexType inputRequestedRegionStartIndex = outputRequestedRegionStartIndex;

  // we but need all of the input in the fft direction
  const unsigned int direction = this->m_Direction;
  const typename TInputImageType::SizeType& inputLargeSize =
    inputPtr->GetLargestPossibleRegion().GetSize();
  inputRequestedRegionSize[direction] = inputLargeSize[direction];
  const typename TInputImageType::IndexType& inputLargeIndex =
    inputPtr->GetLargestPossibleRegion().GetIndex();
  inputRequestedRegionStartIndex[direction] = inputLargeIndex[direction];

  typename TInputImageType::RegionType inputRequestedRegion;
  inputRequestedRegion.SetSize( inputRequestedRegionSize );
  inputRequestedRegion.SetIndex( inputRequestedRegionStartIndex );

  inputPtr->SetRequestedRegion( inputRequestedRegion );
}


template < class TPixel , unsigned int Dimension >
void
FFT1DComplexConjugateToRealImageFilter < TPixel , Dimension >
::EnlargeOutputRequestedRegion(DataObject *output)
{
  TOutputImageType* outputPtr = dynamic_cast<TOutputImageType*>( output );

  // we need to enlarge the region in the fft direction to the
  // largest possible in that direction
  typedef const typename TOutputImageType::SizeType& ConstOutputSizeType;
  ConstOutputSizeType requestedSize =
    outputPtr->GetRequestedRegion().GetSize();
  ConstOutputSizeType outputLargeSize =
    outputPtr->GetLargestPossibleRegion().GetSize();
  typedef const typename TOutputImageType::IndexType& ConstOutputIndexType;
  ConstOutputIndexType requestedIndex =
    outputPtr->GetRequestedRegion().GetIndex();
  ConstOutputIndexType outputLargeIndex =
    outputPtr->GetLargestPossibleRegion().GetIndex();

  typename TOutputImageType::SizeType enlargedSize   = requestedSize;
  typename TOutputImageType::IndexType enlargedIndex = requestedIndex;
  enlargedSize[this->m_Direction]  = outputLargeSize[this->m_Direction];
  enlargedIndex[this->m_Direction] = outputLargeIndex[this->m_Direction];

  typename TOutputImageType::RegionType enlargedRegion;
  enlargedRegion.SetSize( enlargedSize );
  enlargedRegion.SetIndex( enlargedIndex );
  outputPtr->SetRequestedRegion( enlargedRegion );
}


template < class TPixel , unsigned int Dimension >
void
FFT1DComplexConjugateToRealImageFilter < TPixel , Dimension >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Direction: " << m_Direction << std::endl;
}


template < class TPixel , unsigned int Dimension >
int 
FFT1DComplexConjugateToRealImageFilter < TPixel , Dimension >
::SplitRequestedRegion(int i, int num, OutputImageRegionType& splitRegion)
{
  // Get the output pointer
  TOutputImageType * outputPtr = this->GetOutput();
  const typename TOutputImageType::SizeType& requestedRegionSize 
    = outputPtr->GetRequestedRegion().GetSize();

  int splitAxis;
  typename TOutputImageType::IndexType splitIndex;
  typename TOutputImageType::SizeType splitSize;

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
  typename TOutputImageType::SizeType::SizeValueType range = requestedRegionSize[splitAxis];
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

#endif // __itkFFT1DComplexConjugateToRealImageFilter_txx
