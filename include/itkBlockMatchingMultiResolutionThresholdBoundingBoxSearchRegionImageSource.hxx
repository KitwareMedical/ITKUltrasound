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
#ifndef itkBlockMatchingMultiResolutionThresholdBoundingBoxSearchRegionImageSource_hxx
#define itkBlockMatchingMultiResolutionThresholdBoundingBoxSearchRegionImageSource_hxx

#include "itkBlockMatchingMultiResolutionThresholdBoundingBoxSearchRegionImageSource.h"

#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkResampleImageFilter.h"

namespace itk
{
namespace BlockMatching
{

template < class TFixedImage, class TMovingImage, class TMetricImage,
         class TDisplacementImage >
MultiResolutionThresholdBoundingBoxSearchRegionImageSource<
TFixedImage, TMovingImage, TMetricImage, TDisplacementImage >
::MultiResolutionThresholdBoundingBoxSearchRegionImageSource():
  m_ThresholdScheduleSpecified( false ),
  m_TopLevelRadiusSpecified( false )
{
  this->m_SearchRegionRadiusImage = nullptr;

  m_DisplacementResampler = DisplacementResamplerType::New();
  m_SearchRegionRadiusResampler = SearchRegionRadiusResamplerType::New();

  // Sane default.
  m_MinimumSearchRegionRadius.Fill( 1 );
}

template < class TFixedImage, class TMovingImage, class TMetricImage,
         class TDisplacementImage >
void
MultiResolutionThresholdBoundingBoxSearchRegionImageSource<
TFixedImage, TMovingImage, TMetricImage, TDisplacementImage >
::SetDisplacementImage( DisplacementImageType * image )
{
  DisplacementCalculatorSuperclass::SetDisplacementImage( image );

  if( m_SearchRegionRadiusImage.GetPointer() == nullptr )
    m_SearchRegionRadiusImage = SearchRegionRadiusImageType::New();
  m_SearchRegionRadiusImage->CopyInformation( image );
  m_SearchRegionRadiusImage->SetRegions( image->GetLargestPossibleRegion() );
  m_SearchRegionRadiusImage->Allocate();
}

//! @todo remove me
template < class TFixedImage, class TMovingImage, class TMetricImage,
         class TDisplacementImage >
void
MultiResolutionThresholdBoundingBoxSearchRegionImageSource<
TFixedImage, TMovingImage, TMetricImage, TDisplacementImage >
::SetThresholdSchedule( const ThresholdType threshold )
{
  // Check to make sure the PyramidSchedule has been set.
  if( this->m_PyramidSchedule.size() == 0 )
    {
    itkExceptionMacro(<<"The PyramidSchedule must be set before calling this method.");
    }

  m_ThresholdSchedule.set_size( this->m_PyramidSchedule.rows() );
  for( unsigned int i = 0; i < this->m_PyramidSchedule.rows(); ++i )
    {
    m_ThresholdSchedule[i] = threshold;
    }
  this->Modified();
  m_ThresholdScheduleSpecified = true;
}

template < class TFixedImage, class TMovingImage, class TMetricImage,
         class TDisplacementImage >
void
MultiResolutionThresholdBoundingBoxSearchRegionImageSource<
TFixedImage, TMovingImage, TMetricImage, TDisplacementImage >
::SetMetricImagePixel( const PointType& point, const IndexType& index,
  MetricImageType* metricImage )
{
  DisplacementCalculatorSuperclass::SetMetricImagePixel( point, index, metricImage );

  if( ! m_ThresholdScheduleSpecified )
    {
    itkExceptionMacro( << "Threshold has not been set." );
    }

  typedef typename MetricImageType::RegionType MetricRegionType;
  typedef typename MetricRegionType::IndexType MetricIndexType;

  // Lower and upper indices of the bounding box above the threshold.
  MetricIndexType upperIndex = metricImage->GetBufferedRegion().GetIndex();
  MetricIndexType lowerIndex;
  unsigned int i;
  for( i = 0; i < ImageDimension; ++i )
    lowerIndex[i] = upperIndex[i] + metricImage->GetBufferedRegion().GetSize()[i];

  typedef itk::ImageLinearConstIteratorWithIndex< MetricImageType > IteratorType;
  IteratorType it( metricImage, metricImage->GetBufferedRegion() );
  for( i = 0; i < ImageDimension; ++i )
    {
    it.SetDirection( i );
    it.GoToBegin();

    while( !it.IsAtEnd() )
      {
      it.GoToBeginOfLine();
      while( ! it.IsAtEndOfLine() )
        {
        if( lowerIndex[i] == it.GetIndex()[i] )
          break;
        if( it.Get() > m_ThresholdSchedule[this->m_CurrentLevel] )
          {
          lowerIndex[i] = it.GetIndex()[i];
          break;
          }
        ++it;
        }

      it.GoToEndOfLine();
      while( ! it.IsAtReverseEndOfLine() )
        {
        if( upperIndex[i] == it.GetIndex()[i] )
          break;
        if( it.Get() > m_ThresholdSchedule[this->m_CurrentLevel] )
          {
          upperIndex[i] = it.GetIndex()[i];
          break;
          }
        --it;
        }

      it.NextLine();
      }
    }

  if( upperIndex[0] < lowerIndex[0] )
    //{
    //itkExceptionMacro( << "None of the metric image values were below the threshold." );
    //}
  // This is a degenerate case that seems to occur when matching signal-less
  // areas to each other.
    {
    // Set to center value.
    lowerIndex = metricImage->GetBufferedRegion().GetIndex();
    for( i = 0; i < ImageDimension; ++i )
      {
      upperIndex[i] = lowerIndex[i] + metricImage->GetBufferedRegion().GetSize()[i];
      lowerIndex[i] = ( upperIndex[i] - lowerIndex[i] ) / 2;
      upperIndex[i] = lowerIndex[i] + 1;
      }
    }

  typename SearchRegionRadiusImageType::PixelType radius;
  MetricIndexType centerIndex;
  typename MetricImageType::SpacingType metricSpacing = metricImage->GetSpacing();
  typename MovingImageType::SpacingType movingSpacing = this->m_MovingImage->GetSpacing();
  for( i = 0; i < ImageDimension; ++i )
    {
    radius[i] = std::ceil(( static_cast< float >( upperIndex[i] ) - static_cast< float >( lowerIndex[i] ) - 1.f ) / 2.f );
    centerIndex[i] = lowerIndex[i] + radius[i];
    // The radius should be specified in terms of the moving image.
    radius[i] = std::ceil( radius[i] * metricSpacing[i] / movingSpacing[i] );
    }

  m_SearchRegionRadiusImage->SetPixel( index, radius );
  // The displacement is taken to be the center of the search region.
  PointType centerPoint;
  metricImage->TransformIndexToPhysicalPoint( centerIndex, centerPoint );
  this->m_DisplacementImage->SetPixel( index, centerPoint - point );
}

template < class TFixedImage, class TMovingImage, class TMetricImage,
         class TDisplacementImage >
void
MultiResolutionThresholdBoundingBoxSearchRegionImageSource<
TFixedImage, TMovingImage, TMetricImage, TDisplacementImage >
::BeforeThreadedGenerateData()
{
  if( this->m_CurrentLevel != 0 )
    {
    //! @todo these resampler should be replaced by resamplers for each
    //component that can specify a neumann boundary condition
    // ditto with FixedSearchRegionImageSource
    m_DisplacementResampler->SetInput( this->m_PreviousDisplacements );
    typename OutputImageType::Pointer outputPtr = this->GetOutput();
    if( !outputPtr )
      return;
    m_DisplacementResampler->SetSize( outputPtr->GetRequestedRegion().GetSize() );
    m_DisplacementResampler->SetOutputStartIndex( outputPtr->GetRequestedRegion().GetIndex() );
    m_DisplacementResampler->SetOutputSpacing( outputPtr->GetSpacing() );
    m_DisplacementResampler->SetOutputOrigin( outputPtr->GetOrigin() );
    m_DisplacementResampler->SetOutputDirection( outputPtr->GetDirection() );
    m_DisplacementResampler->UpdateLargestPossibleRegion();

    m_SearchRegionRadiusResampler->SetInput( this->m_SearchRegionRadiusImage );
    m_SearchRegionRadiusResampler->SetSize( outputPtr->GetRequestedRegion().GetSize() );
    m_SearchRegionRadiusResampler->SetOutputStartIndex( outputPtr->GetRequestedRegion().GetIndex() );
    m_SearchRegionRadiusResampler->SetOutputSpacing( outputPtr->GetSpacing() );
    m_SearchRegionRadiusResampler->SetOutputOrigin( outputPtr->GetOrigin() );
    m_SearchRegionRadiusResampler->SetOutputDirection( outputPtr->GetDirection() );
    m_SearchRegionRadiusResampler->UpdateLargestPossibleRegion();
    }

  if( ! m_TopLevelRadiusSpecified )
    {
    itkExceptionMacro( << "The top level radius has not been specified." );
    }
}

template < class TFixedImage, class TMovingImage, class TMetricImage,
         class TDisplacementImage >
void
MultiResolutionThresholdBoundingBoxSearchRegionImageSource<
TFixedImage, TMovingImage, TMetricImage, TDisplacementImage >
::ThreadedGenerateData( const OutputRegionType& outputRegion, ThreadIdType threadID )
{
  OutputImageType * outputPtr = this->GetOutput();
  if( !outputPtr )
    {
    return;
    }

  OutputRegionType region;
  OutputRegionType movingLargestRegion = this->m_MovingImage->GetLargestPossibleRegion();
  typename MovingImageType::PointType point;
  typename MovingImageType::IndexType index;
  typename MovingImageType::SizeType unitySize;
  unitySize.Fill( 1 );
  RadiusType radius;

  if( this->m_CurrentLevel == 0 )
    {
    ImageRegionIteratorWithIndex< OutputImageType > it( outputPtr, outputRegion );
    for( it.GoToBegin(); !it.IsAtEnd(); ++it )
      {
      index = it.GetIndex();
      outputPtr->TransformIndexToPhysicalPoint( index, point );
      this->m_MovingImage->TransformPhysicalPointToIndex( point, index );
      region.SetIndex( index );
      region.SetSize( unitySize );
      region.PadByRadius( m_TopLevelRadius );
      if( !region.Crop( movingLargestRegion ) )
        {
        itkExceptionMacro( << "Attempted to create a search region entirely outside the moving image." );
        }
      it.Set( region );
      }
    }
  else
    {
    // Cached for determining the closest index when outside the input.
    IndexType startIndex = this->m_MovingImage->GetBufferedRegion().GetIndex();
    IndexType endIndex;
    IndexType closestIndex;
    for( unsigned int i = 0; i < ImageDimension; ++i )
      {
      endIndex[i] = startIndex[i] + this->m_MovingImage->GetBufferedRegion().GetSize()[i] - 1;
      }

    ImageRegionIteratorWithIndex< OutputImageType > it( outputPtr, outputRegion );
    ImageRegionConstIterator< DisplacementImageType > dispIt( m_DisplacementResampler->GetOutput(),
      outputRegion );
    ImageRegionConstIterator< SearchRegionRadiusImageType > radiusIt( m_SearchRegionRadiusResampler->GetOutput(), outputRegion );

    for( it.GoToBegin(), dispIt.GoToBegin(), radiusIt.GoToBegin();
      !it.IsAtEnd();
      ++it, ++dispIt, ++radiusIt )
      {
      index = it.GetIndex();
      outputPtr->TransformIndexToPhysicalPoint( index, point );
      // resample displacement image
      this->m_MovingImage->TransformPhysicalPointToIndex( point + dispIt.Get(),
        index );
      region.SetIndex( index );
      region.SetSize( unitySize );
      // Expand the search region radius by the scaling that occured between
      // levels.
      for( unsigned int i = 0; i < ImageDimension; ++i )
        {
        radius[i] = std::ceil( static_cast< float >( radiusIt.Get()[i] ) *
          static_cast< float >( this->m_PyramidSchedule( this->m_CurrentLevel - 1, i )) /
          static_cast< float >( this->m_PyramidSchedule( this->m_CurrentLevel, i )));
        if( radius[i] < m_MinimumSearchRegionRadius[i] )
          radius[i] = m_MinimumSearchRegionRadius[i];
        }
      region.PadByRadius( radius );
      if( !region.Crop( movingLargestRegion ) )
        {
        // Set to the closest index and size one.
        for( unsigned int i = 0; i < ImageDimension; ++i )
          {
          closestIndex[i] = index[i];
          }
        for( unsigned int i = 0; i < ImageDimension; ++i )
          {
          if( index[i] < startIndex[i] )
            {
            closestIndex[i] = startIndex[i];
            }
          if( index[i] > endIndex[i] )
            {
            closestIndex[i] = endIndex[i];
            }
          }
        region.SetIndex( closestIndex );
        region.SetSize( unitySize );
        }
      it.Set( region );
      }
    }
}

} // end namespace itk
} // end namespace BlockMatching

#endif
