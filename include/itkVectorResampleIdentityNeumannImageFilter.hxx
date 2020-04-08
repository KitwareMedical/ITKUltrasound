/*=========================================================================
 *
 *  Copyright NumFOCUS
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
#ifndef itkVectorResampleIdentityNeumannImageFilter_hxx
#define itkVectorResampleIdentityNeumannImageFilter_hxx

#include "itkMath.h"
#include "itkVectorResampleIdentityNeumannImageFilter.h"
#include "itkObjectFactory.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkImageRegionIteratorWithIndex.h"

namespace itk
{

template <typename TInputImage, typename TOutputImage, typename TInterpolatorPrecisionType>
VectorResampleIdentityNeumannImageFilter<TInputImage, TOutputImage, TInterpolatorPrecisionType>
::VectorResampleIdentityNeumannImageFilter()
{
  m_OutputSpacing.Fill(1.0);
  m_OutputOrigin.Fill(0.0);
  m_OutputDirection.SetIdentity();
  m_Size.Fill( 0 );
  m_OutputStartIndex.Fill( 0 );

  m_Interpolator = VectorLinearInterpolateImageFunction<InputImageType, TInterpolatorPrecisionType>::New();
}


template <typename TInputImage, typename TOutputImage, typename TInterpolatorPrecisionType>
void
VectorResampleIdentityNeumannImageFilter<TInputImage, TOutputImage,TInterpolatorPrecisionType>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "Size: " << m_Size << std::endl;
  os << indent << "OutputStartIndex: " << m_OutputStartIndex << std::endl;
  os << indent << "OutputSpacing: " << m_OutputSpacing << std::endl;
  os << indent << "OutputOrigin: " << m_OutputOrigin << std::endl;
  os << indent << "OutputDirection: " << m_OutputDirection << std::endl;
  os << indent << "Interpolator: " << m_Interpolator.GetPointer() << std::endl;

  return;
}


template <typename TInputImage, typename TOutputImage, typename TInterpolatorPrecisionType>
void
VectorResampleIdentityNeumannImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::SetOutputSpacing(const double* spacing)
{
  SpacingType s(spacing);
  this->SetOutputSpacing( s );
}


template <typename TInputImage, typename TOutputImage, typename TInterpolatorPrecisionType>
void
VectorResampleIdentityNeumannImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::SetOutputOrigin(const double* origin)
{
  PointType p(origin);
  this->SetOutputOrigin( p );
}

/**
 * Set up state of filter before multi-threading.
 * InterpolatorType::SetInputImage is not thread-safe and hence
 * has to be set up before ThreadedGenerateData
 */
template <typename TInputImage, typename TOutputImage, typename TInterpolatorPrecisionType>
void
VectorResampleIdentityNeumannImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::BeforeThreadedGenerateData()
{

  if( !m_Interpolator )
    {
    itkExceptionMacro(<< "Interpolator not set");
    }

  // Connect input image to interpolator
  m_Interpolator->SetInputImage( this->GetInput() );

}

/**
 * Set up state of filter after multi-threading.
 */
template <typename TInputImage, typename TOutputImage, typename TInterpolatorPrecisionType>
void
VectorResampleIdentityNeumannImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::AfterThreadedGenerateData()
{
  // Disconnect input image from the interpolator
  m_Interpolator->SetInputImage( nullptr );

}


template <typename TInputImage, typename TOutputImage, typename TInterpolatorPrecisionType>
void
VectorResampleIdentityNeumannImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::DynamicThreadedGenerateData( const OutputImageRegionType& outputRegionForThread )
{
  // Get the output pointers
  OutputImagePointer      outputPtr = this->GetOutput();

  // Get ths input pointers
  InputImageConstPointer inputPtr=this->GetInput();

  // Create an iterator that will walk the output region for this thread.
  using OutputIterator = ImageRegionIteratorWithIndex<TOutputImage>;

  OutputIterator outIt(outputPtr, outputRegionForThread);

  // Define a few indices that will be used to translate from an input pixel
  // to an output pixel
  PointType outputPoint;         // Coordinates of current output pixel
  PointType inputPoint;          // Coordinates of current input pixel

  using ContinuousIndexType = ContinuousIndex<TInterpolatorPrecisionType, ImageDimension>;
  ContinuousIndexType inputIndex;

  const unsigned int numberOfComponents = PixelType::GetNumberOfComponents();

  using OutputType = typename InterpolatorType::OutputType;

  // Cached for determining the closest index when outside the input.
  IndexType startIndex = inputPtr->GetBufferedRegion().GetIndex();
  IndexType endIndex;
  IndexType closestIndex;
  for( unsigned int i = 0; i < ImageDimension; ++i )
    {
    endIndex[i] = startIndex[i] + inputPtr->GetBufferedRegion().GetSize()[i] - 1;
    }

  // Walk the output region
  outIt.GoToBegin();

  while ( !outIt.IsAtEnd() )
    {
    // Determine the index of the current output pixel
    outputPtr->TransformIndexToPhysicalPoint( outIt.GetIndex(), outputPoint );

    // Compute corresponding input pixel position
    inputPoint = outputPoint;
    inputPtr->TransformPhysicalPointToContinuousIndex(inputPoint, inputIndex);

    // Evaluate input at right position and copy to the output
    if( m_Interpolator->IsInsideBuffer(inputIndex) )
      {
      PixelType   pixval;
      const OutputType  value
        = m_Interpolator->EvaluateAtContinuousIndex( inputIndex );
      for( unsigned int i=0; i< numberOfComponents; i++ )
        {
        pixval[i] = static_cast<PixelComponentType>( value[i] );
        }
      outIt.Set( pixval );
      }
    else
      {
      for( unsigned int i = 0; i < ImageDimension; ++i )
        closestIndex[i] = itk::Math::RoundHalfIntegerToEven< typename IndexType::IndexValueType, double >( inputIndex[i] );
      for( unsigned int i = 0; i < ImageDimension; ++i )
        {
        if( inputIndex[i] < startIndex[i] )
          {
          closestIndex[i] = startIndex[i];
          }
        if( inputIndex[i] > endIndex[i] )
          {
          closestIndex[i] = endIndex[i];
          }
        }
      outIt.Set( inputPtr->GetPixel( closestIndex )); // default background value
      }

    ++outIt;
    }
  return;
}

/**
 * Inform pipeline of necessary input image region
 *
 * Determining the actual input region is non-trivial, especially
 * when we cannot assume anything about the transform being used.
 * So we do the easy thing and request the entire input image.
 */
template <typename TInputImage, typename TOutputImage, typename TInterpolatorPrecisionType>
void
VectorResampleIdentityNeumannImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::GenerateInputRequestedRegion()
{
  // call the superclass's implementation of this method
  Superclass::GenerateInputRequestedRegion();

  if ( !this->GetInput() )
    {
    return;
    }

  // get pointers to the input and output
  InputImagePointer  inputPtr  =
    const_cast< TInputImage *>( this->GetInput() );

  // Request the entire input image
  InputImageRegionType inputRegion;
  inputRegion = inputPtr->GetLargestPossibleRegion();
  inputPtr->SetRequestedRegion(inputRegion);

  return;
}


template <typename TInputImage, typename TOutputImage, typename TInterpolatorPrecisionType>
void
VectorResampleIdentityNeumannImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::GenerateOutputInformation()
{
  // call the superclass' implementation of this method
  Superclass::GenerateOutputInformation();

  // get pointers to the input and output
  OutputImagePointer outputPtr = this->GetOutput();
  if ( !outputPtr )
    {
    return;
    }

  // Set the size of the output region
  typename TOutputImage::RegionType outputLargestPossibleRegion;
  outputLargestPossibleRegion.SetSize( m_Size );
  outputLargestPossibleRegion.SetIndex( m_OutputStartIndex );
  outputPtr->SetLargestPossibleRegion( outputLargestPossibleRegion );

  // Set spacing and origin
  outputPtr->SetSpacing( m_OutputSpacing );
  outputPtr->SetOrigin( m_OutputOrigin );
  outputPtr->SetDirection( m_OutputDirection );

  return;
}

/**
 * Verify if any of the components has been modified.
 */
template <typename TInputImage, typename TOutputImage, typename TInterpolatorPrecisionType>
ModifiedTimeType
VectorResampleIdentityNeumannImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::GetMTime() const
{
  unsigned long latestTime = Object::GetMTime();

  if( m_Interpolator )
    {
    if( latestTime < m_Interpolator->GetMTime() )
      {
      latestTime = m_Interpolator->GetMTime();
      }
    }

  return latestTime;
}

} // end namespace itk

#endif
