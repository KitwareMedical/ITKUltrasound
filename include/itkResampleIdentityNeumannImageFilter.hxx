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
#ifndef itkResampleIdentityNeumannImageFilter_hxx
#define itkResampleIdentityNeumannImageFilter_hxx

#include "itkResampleIdentityNeumannImageFilter.h"
#include "itkObjectFactory.h"
#include "itkIdentityTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkProgressReporter.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkSpecialCoordinatesImage.h"

namespace itk
{

template <class TInputImage,
          class TOutputImage,
          class TInterpolatorPrecisionType>
ResampleIdentityNeumannImageFilter<TInputImage, TOutputImage,TInterpolatorPrecisionType>
::ResampleIdentityNeumannImageFilter()
{
  m_OutputOrigin.Fill(0.0);
  m_OutputSpacing.Fill(1.0);
  m_OutputDirection.SetIdentity();

  m_UseReferenceImage = false;

  m_Size.Fill( 0 );
  m_OutputStartIndex.Fill( 0 );

  m_Interpolator = LinearInterpolateImageFunction<InputImageType,
                                  TInterpolatorPrecisionType>::New();
}


template <class TInputImage,
          class TOutputImage,
          class TInterpolatorPrecisionType>
void
ResampleIdentityNeumannImageFilter<TInputImage, TOutputImage,TInterpolatorPrecisionType>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "Size: " << m_Size << std::endl;
  os << indent << "OutputStartIndex: " << m_OutputStartIndex << std::endl;
  os << indent << "OutputSpacing: " << m_OutputSpacing << std::endl;
  os << indent << "OutputOrigin: " << m_OutputOrigin << std::endl;
  os << indent << "OutputDirection: " << m_OutputDirection << std::endl;
  os << indent << "Interpolator: " << m_Interpolator.GetPointer() << std::endl;
  os << indent << "UseReferenceImage: " << (m_UseReferenceImage ? "On" : "Off")
               << std::endl;
  return;
}


template <class TInputImage,
          class TOutputImage,
          class TInterpolatorPrecisionType>
void
ResampleIdentityNeumannImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::SetOutputSpacing( const double* spacing )
{
  SpacingType s(spacing);
  this->SetOutputSpacing( s );
}


template <class TInputImage,
          class TOutputImage,
          class TInterpolatorPrecisionType>
void
ResampleIdentityNeumannImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::SetOutputOrigin( const double* origin )
{
  OriginPointType p(origin);
  this->SetOutputOrigin( p );
}

/** Helper method to set the output parameters based on this image */
template <class TInputImage,
          class TOutputImage,
          class TInterpolatorPrecisionType>
void
ResampleIdentityNeumannImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::SetOutputParametersFromImage ( const ImageBaseType * image )
{
  this->SetOutputOrigin ( image->GetOrigin() );
  this->SetOutputSpacing ( image->GetSpacing() );
  this->SetOutputDirection ( image->GetDirection() );
  this->SetOutputStartIndex ( image->GetLargestPossibleRegion().GetIndex() );
  this->SetSize ( image->GetLargestPossibleRegion().GetSize() );
}


template <class TInputImage,
          class TOutputImage,
          class TInterpolatorPrecisionType>
void
ResampleIdentityNeumannImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::BeforeThreadedGenerateData()
{
  if( !m_Interpolator )
    {
    itkExceptionMacro(<< "Interpolator not set");
    }

  // Connect input image to interpolator
  m_Interpolator->SetInputImage( this->GetInput() );
}


template <class TInputImage,
          class TOutputImage,
          class TInterpolatorPrecisionType>
void
ResampleIdentityNeumannImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::AfterThreadedGenerateData()
{
  // Disconnect input image from the interpolator
  m_Interpolator->SetInputImage( ITK_NULLPTR );
}


template <class TInputImage,
          class TOutputImage,
          class TInterpolatorPrecisionType>
void
ResampleIdentityNeumannImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread,
                        ThreadIdType threadId)
{
  // Get the output pointers
  OutputImagePointer      outputPtr = this->GetOutput();

  // Get ths input pointers
  InputImageConstPointer inputPtr=this->GetInput();

  // Create an iterator that will walk the output region for this thread.
  typedef ImageLinearIteratorWithIndex<TOutputImage> OutputIterator;

  OutputIterator outIt(outputPtr, outputRegionForThread);
  outIt.SetDirection( 0 );

  // Define a few indices that will be used to translate from an input pixel
  // to an output pixel
  PointType outputPoint;         // Coordinates of current output pixel
  PointType inputPoint;          // Coordinates of current input pixel
  PointType tmpOutputPoint;
  PointType tmpInputPoint;

  typedef ContinuousIndex<TInterpolatorPrecisionType, ImageDimension>
                                                            ContinuousIndexType;
  ContinuousIndexType inputIndex;
  ContinuousIndexType tmpInputIndex;

  typedef typename PointType::VectorType VectorType;
  VectorType delta;          // delta in input continuous index coordinate frame
  IndexType index;

  // Support for progress methods/callbacks
  ProgressReporter progress(this,
                            threadId,
                            outputRegionForThread.GetNumberOfPixels());

  typedef typename InterpolatorType::OutputType OutputType;

  // Cached for determining the closest index when outside the input.
  IndexType startIndex = inputPtr->GetBufferedRegion().GetIndex();
  IndexType endIndex;
  IndexType closestIndex;
  for( unsigned int i = 0; i < ImageDimension; ++i )
    {
    endIndex[i] = startIndex[i] + inputPtr->GetBufferedRegion().GetSize()[i] - 1;
    }

  // Min/max values of the output pixel type AND these values
  // represented as the output type of the interpolator
  const PixelType minValue =  NumericTraits<PixelType >::NonpositiveMin();
  const PixelType maxValue =  NumericTraits<PixelType >::max();

  const OutputType minOutputValue = static_cast<OutputType>(minValue);
  const OutputType maxOutputValue = static_cast<OutputType>(maxValue);


  // Determine the position of the first pixel in the scanline
  index = outIt.GetIndex();
  outputPtr->TransformIndexToPhysicalPoint( index, outputPoint );


  // Compute corresponding input pixel position
  inputPoint = outputPoint;
  inputPtr->TransformPhysicalPointToContinuousIndex(inputPoint, inputIndex);

  // As we walk across a scan line in the output image, we trace
  // an oriented/scaled/translated line in the input image.  Cache
  // the delta along this line in continuous index space of the input
  // image. This allows us to use vector addition to model the
  // transformation.
  //
  // To find delta, we take two pixels adjacent in a scanline
  // and determine the continuous indices of these pixels when
  // mapped to the input coordinate frame. We use the difference
  // between these two continuous indices as the delta to apply
  // to an index to trace line in the input image as we move
  // across the scanline of the output image.
  //
  // We determine delta in this manner so that Images
  // are both handled properly (taking into account the direction cosines).
  //
  ++index[0];
  outputPtr->TransformIndexToPhysicalPoint( index, tmpOutputPoint );
  tmpInputPoint = tmpOutputPoint;
  inputPtr->TransformPhysicalPointToContinuousIndex(tmpInputPoint,
                                                    tmpInputIndex);
  delta = tmpInputIndex - inputIndex;


  // This fix works for images up to approximately 2^25 pixels in
  // any dimension.  If the image is larger than this, this constant
  // needs to be made lower.
  double precisionConstant = 1<<(NumericTraits<double>::digits>>1);

  // Delta is precise to many decimal points, but this precision
  // involves some error in the last bits.  This error can accumulate
  // as the delta values are added.
  // Sometimes, when the accumulated delta should be inside of the
  // image, it will be slightly
  // greater than the largest index in the image, like 255.00000000002
  // for a image of size 256.  This can cause an empty column to show up
  // at the right side of the image. If we instead
  // truncate this delta value to some precision, this solves the problem.
  // Therefore, the following routine uses a
  // precisionConstant that specifies the number of relevant bits,
  // and the value is truncated to this precision.
  for( unsigned int i = 0; i < ImageDimension; ++i)
    {
    long roundedInputIndex = (long)(inputIndex[i]);
    if(inputIndex[i]<0.0 && inputIndex[i]!=(double)roundedInputIndex)
      {
      --roundedInputIndex;
      }
    double inputIndexFrac = inputIndex[i] - roundedInputIndex;
    double newInputIndexFrac = (long)(precisionConstant
                                          * inputIndexFrac)
                               / precisionConstant;
    inputIndex[i] = roundedInputIndex + newInputIndexFrac;
    }


  while ( !outIt.IsAtEnd() )
    {
    // Determine the continuous index of the first pixel of output
    // scanline when mapped to the input coordinate frame.
    //

    // First get the position of the pixel in the output coordinate frame
    index = outIt.GetIndex();
    outputPtr->TransformIndexToPhysicalPoint( index, outputPoint );

    // Compute corresponding input pixel continuous index, this index
    // will incremented in the scanline loop
    inputPoint = outputPoint;
    inputPtr->TransformPhysicalPointToContinuousIndex(inputPoint, inputIndex);

    // The inputIndex is precise to many decimal points, but this precision
    // involves some error in the last bits.
    // Sometimes, when an index should be inside of the image, the
    // index will be slightly
    // greater than the largest index in the image, like 255.00000000002
    // for a image of size 256.  This can cause an empty row to show up
    // at the bottom of the image.
    // Therefore, the following routine uses a
    // precisionConstant that specifies the number of relevant bits,
    // and the value is truncated to this precision.
    for( unsigned int i = 0; i < ImageDimension; ++i)
      {
      long roundedInputIndex = (long)(inputIndex[i]);
      if(inputIndex[i]<0.0 && inputIndex[i]!=(double)roundedInputIndex)
        {
        --roundedInputIndex;
        }
      double inputIndexFrac = inputIndex[i] - roundedInputIndex;
      double newInputIndexFrac = (long)(precisionConstant
                                            * inputIndexFrac)
                                 / precisionConstant;
      inputIndex[i] = roundedInputIndex + newInputIndexFrac;
      }

    if(  !outIt.IsAtEnndOfLine() &&
         !m_Interpolator->IsInsideBuffer(inputIndex) )
      {
      for( unsigned int i = 0; i < ImageDimension; ++i )
        {
        closestIndex[i] = Math::Round( inputIndex[i] );
        }
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
      }
    while( !outIt.IsAtEndOfLine() &&
           !m_Interpolator->IsInsideBuffer(inputIndex) )
      {
      outIt.Set(closestIndex);
      progress.CompletedPixel();
      ++outIt;
      inputIndex += delta;
      }

    if( !outIt.IsAtEndOfLine() && m_Interpolator->IsInsideBuffer(inputIndex) )
      {
      PixelType pixval;
      OutputType value;
      value = m_Interpolator->EvaluateAtContinuousIndex(inputIndex);
      if( value <  minOutputValue )
        {
        pixval = minValue;
        }
      else if( value > maxOutputValue )
        {
        pixval = maxValue;
        }
      else
        {
        pixval = static_cast<PixelType>( value );
        }
      outIt.Set( pixval );

      progress.CompletedPixel();
      ++outIt;
      inputIndex += delta;
      }

    while( !outIt.IsAtEndOfLine() )
      {
      // Evaluate input at right position and copy to the output
      if( m_Interpolator->IsInsideBuffer(inputIndex) )
        {
        PixelType pixval;
        OutputType value;
        value = m_Interpolator->EvaluateAtContinuousIndex(inputIndex);
        if( value <  minOutputValue )
          {
          pixval = minValue;
          }
        else if( value > maxOutputValue )
          {
          pixval = maxValue;
          }
        else
          {
          pixval = static_cast<PixelType>( value );
          }
        outIt.Set( pixval );
        }
      else
        {
        for( unsigned int i = 0; i < ImageDimension; ++i )
          closestIndex[i] = Math::Round( inputIndex[i] );
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
        break;
        }

      progress.CompletedPixel();
      ++outIt;
      inputIndex += delta;
      }

    while( !outIt.IsAtEndOfLine() )
      {
      outIt.Set(closestIndex); // default background value
      progress.CompletedPixel();
      ++outIt;
      inputIndex += delta;
      }

    outIt.NextLine();
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
template <class TInputImage,
          class TOutputImage,
          class TInterpolatorPrecisionType>
void
ResampleIdentityNeumannImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
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
  inputPtr->SetRequestedRegionToLargestPossibleRegion();

  return;
}


/**
 * Set the smart pointer to the reference image that will provide
 * the grid parameters for the output image.
 */
template <class TInputImage,
          class TOutputImage,
          class TInterpolatorPrecisionType>
const typename ResampleIdentityNeumannImageFilter<TInputImage,
                                      TOutputImage,
                                      TInterpolatorPrecisionType>
               ::OutputImageType *
ResampleIdentityNeumannImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::GetReferenceImage() const
{
  Self * surrogate = const_cast< Self * >( this );
  const OutputImageType * referenceImage =
    static_cast<const OutputImageType *>(surrogate->ProcessObject::GetInput(1));
  return referenceImage;
}


/**
 * Set the smart pointer to the reference image that will provide
 * the grid parameters for the output image.
 */
template <class TInputImage,
          class TOutputImage,
          class TInterpolatorPrecisionType>
void
ResampleIdentityNeumannImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::SetReferenceImage( const TOutputImage *image )
{
  itkDebugMacro("setting input ReferenceImage to " << image);
  if( image != static_cast<const TOutputImage *>(this->ProcessObject::GetInput( 1 )) )
    {
    this->ProcessObject::SetNthInput(1, const_cast< TOutputImage *>( image ) );
    this->Modified();
    }
}

/**
 * Inform pipeline of required output region
 */
template <class TInputImage,
          class TOutputImage,
          class TInterpolatorPrecisionType>
void
ResampleIdentityNeumannImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
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

  const OutputImageType * referenceImage = this->GetReferenceImage();

  // Set the size of the output region
  if( m_UseReferenceImage && referenceImage )
    {
    outputPtr->SetLargestPossibleRegion(
                                  referenceImage->GetLargestPossibleRegion() );
    }
  else
    {
    typename TOutputImage::RegionType outputLargestPossibleRegion;
    outputLargestPossibleRegion.SetSize( m_Size );
    outputLargestPossibleRegion.SetIndex( m_OutputStartIndex );
    outputPtr->SetLargestPossibleRegion( outputLargestPossibleRegion );
    }

  // Set spacing and origin
  if (m_UseReferenceImage && referenceImage)
    {
    outputPtr->SetSpacing( referenceImage->GetSpacing() );
    outputPtr->SetOrigin( referenceImage->GetOrigin() );
    outputPtr->SetDirection( referenceImage->GetDirection() );
    }
  else
    {
    outputPtr->SetSpacing( m_OutputSpacing );
    outputPtr->SetOrigin( m_OutputOrigin );
    outputPtr->SetDirection( m_OutputDirection );
    }

  return;
}


/**
 * Verify if any of the components has been modified.
 */
template <class TInputImage,
          class TOutputImage,
          class TInterpolatorPrecisionType>
unsigned long
ResampleIdentityNeumannImageFilter<TInputImage,TOutputImage,TInterpolatorPrecisionType>
::GetMTime( void ) const
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
