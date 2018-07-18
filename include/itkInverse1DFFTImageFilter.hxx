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
#ifndef itkInverse1DFFTImageFilter_hxx
#define itkInverse1DFFTImageFilter_hxx

#include "itkInverse1DFFTImageFilter.h"

#include "itkVnlInverse1DFFTImageFilter.h"

#if defined(ITK_USE_FFTWD) || defined(ITK_USE_FFTWF)
#include "itkFFTWInverse1DFFTImageFilter.h"
#endif

#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"

namespace itk
{

template< typename TSelfPointer, typename TInputImage, typename TOutputImage, typename TPixel >
struct Dispatch_1DComplexConjugateToReal_New
{
  static TSelfPointer Apply()
    {
    return VnlInverse1DFFTImageFilter< TInputImage, TOutputImage >
      ::New().GetPointer();
    }
};


#ifdef ITK_USE_FFTWD
template< typename TSelfPointer, typename TInputImage, typename TOutputImage >
struct Dispatch_1DComplexConjugateToReal_New< TSelfPointer, TInputImage, TOutputImage, double >
{
  static TSelfPointer Apply()
    {
    return FFTWInverse1DFFTImageFilter< TInputImage, TOutputImage >
      ::New().GetPointer();
    }
};
#endif


#ifdef ITK_USE_FFTWF
template< typename TSelfPointer, typename TInputImage, typename TOutputImage >
struct Dispatch_1DComplexConjugateToReal_New< TSelfPointer, TInputImage, TOutputImage, float >
{
  static TSelfPointer Apply()
    {
    return FFTWInverse1DFFTImageFilter< TInputImage, TOutputImage >
      ::New().GetPointer();
    }
};
#endif


template < typename TInputImage, typename TOutputImage >
typename Inverse1DFFTImageFilter< TInputImage, TOutputImage >::Pointer
Inverse1DFFTImageFilter< TInputImage, TOutputImage >
::New()
{
  Pointer smartPtr = ObjectFactory< Self >::Create();

  if( smartPtr.IsNull() )
    {
    smartPtr = Dispatch_1DComplexConjugateToReal_New< Pointer, TInputImage, TOutputImage, typename NumericTraits< typename TOutputImage::PixelType >::ValueType >::Apply();
    }

  return smartPtr;
}


template < typename TInputImage, typename TOutputImage >
Inverse1DFFTImageFilter< TInputImage, TOutputImage >
::Inverse1DFFTImageFilter():
  m_Direction( 0 )
{
  this->m_ImageRegionSplitter = ImageRegionSplitterDirection::New();
}


template <typename TInputImage, typename TOutputImage>
const ImageRegionSplitterBase*
Inverse1DFFTImageFilter< TInputImage, TOutputImage >
::GetImageRegionSplitter() const
{
  return this->m_ImageRegionSplitter.GetPointer();
}


template<typename TInputImage, typename TOutputImage >
void
Inverse1DFFTImageFilter< TInputImage, TOutputImage >
::BeforeThreadedGenerateData()
{
  this->m_ImageRegionSplitter->SetDirection( this->GetDirection() );
  this->DynamicMultiThreadingOff();
}


template< typename TInputImage, typename TOutputImage >
void
Inverse1DFFTImageFilter< TInputImage, TOutputImage >
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


template< typename TInputImage, typename TOutputImage >
void
Inverse1DFFTImageFilter< TInputImage, TOutputImage >
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


template< typename TInputImage, typename TOutputImage >
void
Inverse1DFFTImageFilter< TInputImage, TOutputImage >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Direction: " << m_Direction << std::endl;
}

} // end namespace itk

#endif // itkInverse1DFFTImageFilter_hxx
