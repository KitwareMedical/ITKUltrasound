#ifndef __itkFFT1DRealToComplexConjugateImageFilter_txx
#define __itkFFT1DRealToComplexConjugateImageFilter_txx

#include "itkFFT1DRealToComplexConjugateImageFilter.h"

#include "itkVnlFFT1DRealToComplexConjugateImageFilter.h"

#if defined(USE_FFTWD) || defined(USE_FFTWF)
#include "itkFFT1DRealToComplexConjugateImageFilter.h"
#endif

namespace itk
{

template < class TPixel, unsigned int VDimension >
typename FFT1DRealToComplexConjugateImageFilter< TPixel, VDimension >::Pointer
FFT1DRealToComplexConjugateImageFilter< TPixel, VDimension >
::New()
{
  Pointer smartPtr = ::itk::ObjectFactory< Self >::Create();

//#ifdef USE_FFTWD
  //if( smartPtr.IsNull() )
    //{
    //if( typeid( TPixel ) == typeid( double ) )
      //{
      //smartPtr = dynamic_cast< Self* >(
	//FFTW1DRealToComplexConjugateImageFilter< double, VDimension >
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
	//FFT1DRealToComplexConjugateImageFilter< float, VDimension >
	//::New().GetPointer() );
      //}
    //}
//#endif

  if( smartPtr.IsNull() )
    {
    smartPtr = VnlFFT1DRealToComplexConjugateImageFilter< TPixel, VDimension >
      ::New().GetPointer();
    }

  return smartPtr;
}

template < class TPixel, unsigned int VDimension >
void
FFT1DRealToComplexConjugateImageFilter< TPixel, VDimension >
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
  // has no meaning in the result of an FFT.
  const typename TInputImageType::SizeType&   inputSize
    = inputPtr->GetLargestPossibleRegion().GetSize();
  const typename TInputImageType::IndexType&  inputStartIndex
    = inputPtr->GetLargestPossibleRegion().GetIndex();
  
  typename TOutputImageType::SizeType     outputSize = inputSize;
  typename TOutputImageType::IndexType    outputStartIndex = inputStartIndex;
  
  //
  // in 4.3.4 of the FFTW documentation, they indicate the size of
  // of a real-to-complex FFT is N * N ... + (N /2+1)
  //                              1   2        d
  // complex numbers.
  // static_cast prob. not necessary but want to make sure integer
  // division is used.
  unsigned int direction = this->m_Direction;
  outputSize[direction] = static_cast<unsigned int>(inputSize[direction])/2 + 1;

  //
  // the halving of the input size hides the actual size of the input.
  // to get the same size image out of the IFFT, need to send it as 
  // Metadata.
  typedef typename TOutputImageType::SizeType::SizeValueType SizeScalarType;
  itk::MetaDataDictionary &OutputDic=outputPtr->GetMetaDataDictionary();
  itk::EncapsulateMetaData<SizeScalarType>(OutputDic,
                                       std::string("FFT_Actual_RealImage_Size"),
                                                     inputSize[direction]);
  typename TOutputImageType::RegionType outputLargestPossibleRegion;
  outputLargestPossibleRegion.SetSize( outputSize );
  outputLargestPossibleRegion.SetIndex( outputStartIndex );
  
  outputPtr->SetLargestPossibleRegion( outputLargestPossibleRegion );
}


template < class TPixel , unsigned int Dimension >
void
FFT1DRealToComplexConjugateImageFilter < TPixel , Dimension >
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
FFT1DRealToComplexConjugateImageFilter < TPixel , Dimension >
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
FFT1DRealToComplexConjugateImageFilter < TPixel , Dimension >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Direction: " << m_Direction << std::endl;
}

} // end namespace itk

#endif // __itkFFT1DRealToComplexConjugateImageFilter_txx
