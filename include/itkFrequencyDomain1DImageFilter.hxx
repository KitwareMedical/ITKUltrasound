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
#ifndef itkFrequencyDomain1DImageFilter_hxx
#define itkFrequencyDomain1DImageFilter_hxx

#include "itkFrequencyDomain1DImageFilter.h"

#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkMetaDataObject.h"
#include "itkMath.h"

namespace itk
{

template< typename TInputImage, typename TOutputImage >
FrequencyDomain1DImageFilter< TInputImage, TOutputImage >
::FrequencyDomain1DImageFilter()
{

  this->SetDirection( 0 );
  this->m_FilterFunction = FrequencyDomain1DFilterFunction::New();
  this->m_ImageRegionSplitter = ImageRegionSplitterDirection::New();
  this->DynamicMultiThreadingOff();
}


template< typename TInputImage, typename TOutputImage >
void
FrequencyDomain1DImageFilter< TInputImage, TOutputImage >
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
  const unsigned int direction = this->GetDirection();
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
FrequencyDomain1DImageFilter< TInputImage, TOutputImage >
::EnlargeOutputRequestedRegion(DataObject *output)
{
  OutputImageType* outputPtr = dynamic_cast< OutputImageType* >( output );

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
  const unsigned int direction = this->GetDirection ();
  enlargedSize[direction]  = outputLargeSize[direction];
  enlargedIndex[direction] = outputLargeIndex[direction];

  typename OutputImageType::RegionType enlargedRegion;
  enlargedRegion.SetSize( enlargedSize );
  enlargedRegion.SetIndex( enlargedIndex );
  outputPtr->SetRequestedRegion( enlargedRegion );
}


template< typename TInputImage, typename TOutputImage >
void
FrequencyDomain1DImageFilter< TInputImage, TOutputImage >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "FilterFunction: " << m_FilterFunction << std::endl;
}


template< typename TInputImage, typename TOutputImage >
const ImageRegionSplitterBase *
FrequencyDomain1DImageFilter< TInputImage, TOutputImage >
::GetImageRegionSplitter() const
{
  return this->m_ImageRegionSplitter.GetPointer();
}


template< typename TInputImage, typename TOutputImage >
void
FrequencyDomain1DImageFilter< TInputImage, TOutputImage >
::BeforeThreadedGenerateData()
{
  this->m_ImageRegionSplitter->SetDirection( this->GetDirection() );

  const InputImageType * inputPtr = this->GetInput();
  const typename OutputImageType::SizeType &inputSize = inputPtr->GetRequestedRegion().GetSize();
  const unsigned int direction = this->GetDirection ();
  const SizeValueType size = inputSize[direction];

  this->m_FilterFunction->SetSignalSize( size );
}


template< typename TInputImage, typename TOutputImage >
void
FrequencyDomain1DImageFilter< TInputImage, TOutputImage >
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, ThreadIdType itkNotUsed( threadId ) )
{
  // get pointers to the input and output
  const InputImageType * inputPtr = this->GetInput();
  OutputImageType * outputPtr = this->GetOutput();

  const unsigned int direction = this->GetDirection ();

  typedef ImageLinearConstIteratorWithIndex< OutputImageType > InputIteratorType;
  typedef ImageLinearIteratorWithIndex< OutputImageType >      OutputIteratorType;
  InputIteratorType inputIt( inputPtr, outputRegionForThread );
  OutputIteratorType outputIt( outputPtr, outputRegionForThread );
  inputIt.SetDirection( direction );
  outputIt.SetDirection( direction );

  // for every fft line
  for( inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd();
    outputIt.NextLine(), inputIt.NextLine() )
    {
    SizeValueType i = 0;
    inputIt.GoToBeginOfLine();
    outputIt.GoToBeginOfLine();
    while( !outputIt.IsAtEndOfLine() )
      {
      outputIt.Set( inputIt.Get() *
            static_cast< typename TInputImage::PixelType>(
                 m_FilterFunction->EvaluateIndex( i ) ) );
      ++outputIt;
      ++inputIt;
      ++i;
      }
    }
}


template< typename TInputImage, typename TOutputImage >
void
FrequencyDomain1DImageFilter< TInputImage, TOutputImage >
::AfterThreadedGenerateData()
{
  this->GraftOutput( this->GetOutput() );
}

} // end namespace itk

#endif // itkFrequencyDomain1DImageFilter_hxx
