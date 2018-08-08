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
#ifndef itkBlockMatchingMultiResolutionSimilarityFunctionSearchRegionImageSource_hxx
#define itkBlockMatchingMultiResolutionSimilarityFunctionSearchRegionImageSource_hxx

#include "itkBlockMatchingMultiResolutionSimilarityFunctionSearchRegionImageSource.h"
#include "itkBlockMatchingMaximumPixelDisplacementCalculator.h"

#include "itkImageConstIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkLinearInterpolateImageFunction.h"

namespace itk
{
namespace BlockMatching
{

template <class TFixedImage, class TMovingImage, class TMetricImage,
          class TDisplacementImage, class TFunctor, class TInterpolatorPrecisionType>
MultiResolutionSimilarityFunctionSearchRegionImageSource<
  TFixedImage, TMovingImage, TMetricImage, TDisplacementImage, TFunctor,
  TInterpolatorPrecisionType>
::MultiResolutionSimilarityFunctionSearchRegionImageSource() :
  m_TopLevelRadiusSpecified( false )
{
  this->m_CacheMetricImage = true;
  this->m_SearchRegionRadiusImage = nullptr;

  m_SearchRegionRadiusResampler = SearchRegionRadiusResamplerType::New();

  // Sane default.
  m_MinimumSearchRegionRadiusFactor.Fill( 1.0 );

  m_DisplacementCalculator =
    MaximumPixelDisplacementCalculator<TMetricImage,
                                       TDisplacementImage>::New();

  m_Interpolator = LinearInterpolateImageFunction<MetricImageType,
                                                  TInterpolatorPrecisionType>::New();
}


template <class TFixedImage, class TMovingImage, class TMetricImage,
          class TDisplacementImage, class TFunctor, class TInterpolatorPrecisionType>
void
MultiResolutionSimilarityFunctionSearchRegionImageSource<
  TFixedImage, TMovingImage, TMetricImage, TDisplacementImage, TFunctor,
  TInterpolatorPrecisionType>
::SetDisplacementImage( DisplacementImageType *image )
{
  DisplacementCalculatorSuperclass::SetDisplacementImage( image );

  if( this->m_SearchRegionRadiusImage.GetPointer() == nullptr )
    {
    this->m_SearchRegionRadiusImage = SearchRegionRadiusImageType::New();
    }
  m_SearchRegionRadiusImage->CopyInformation( image );
  m_SearchRegionRadiusImage->SetRegions( image->GetLargestPossibleRegion() );
  m_SearchRegionRadiusImage->Allocate();

  m_DisplacementCalculator->SetDisplacementImage( image );
}


template <class TFixedImage, class TMovingImage, class TMetricImage,
          class TDisplacementImage, class TFunctor, class TInterpolatorPrecisionType>
void
MultiResolutionSimilarityFunctionSearchRegionImageSource<
  TFixedImage, TMovingImage, TMetricImage, TDisplacementImage, TFunctor,
  TInterpolatorPrecisionType>
::SetMetricImagePixel( const PointType & point, const IndexType & index, MetricImageType *metricImage )
{
  DisplacementCalculatorSuperclass::SetMetricImagePixel( point, index, metricImage );

  m_DisplacementCalculator->SetMetricImagePixel( point, index, metricImage );
}


template <class TFixedImage, class TMovingImage, class TMetricImage,
          class TDisplacementImage, class TFunctor, class TInterpolatorPrecisionType>
void
MultiResolutionSimilarityFunctionSearchRegionImageSource<
  TFixedImage, TMovingImage, TMetricImage, TDisplacementImage, TFunctor,
  TInterpolatorPrecisionType>
::Compute()
{
  m_DisplacementCalculator->Compute();

  OutputImageType * outputPtr = this->GetOutput();
  if( !outputPtr )
    {
    return;
    }

  PointType                                       point;
  typename MetricImageType::PixelType             metric;
  typename SearchRegionRadiusImageType::PixelType radius;
  typename OutputRegionType::SizeType             size;

  // Iterate over the displacements and get the corresponding metric image value
  // at the given point.
  ImageRegionConstIteratorWithIndex<MetricImageImageType>
  metricImageImageConstIt( this->m_MetricImageImage,
                           this->m_MetricImageImage->GetBufferedRegion() );
  ImageRegionConstIterator<CenterPointsImageType>
  centerPointsConstIt( this->m_CenterPointsImage,
                       this->m_CenterPointsImage->GetBufferedRegion() );
  ImageRegionConstIterator<DisplacementImageType>
  displacementIt( this->m_DisplacementImage,
                  this->m_DisplacementImage->GetBufferedRegion() );
  ImageRegionIterator<OutputImageType>
  previousSearchRegionIt( outputPtr,
                          outputPtr->GetBufferedRegion() );
  ImageRegionIterator<SearchRegionRadiusImageType>
  searchRegionRadiusIt( this->m_SearchRegionRadiusImage,
                        this->m_SearchRegionRadiusImage->GetBufferedRegion() );
  for( metricImageImageConstIt.GoToBegin(),
       centerPointsConstIt.GoToBegin(),
       displacementIt.GoToBegin(),
       previousSearchRegionIt.GoToBegin(),
       searchRegionRadiusIt.GoToBegin();
       !metricImageImageConstIt.IsAtEnd();
       ++metricImageImageConstIt,
       ++centerPointsConstIt,
       ++displacementIt,
       ++previousSearchRegionIt,
       ++searchRegionRadiusIt )
    {
    point = centerPointsConstIt.Get() + displacementIt.Get();
    m_Interpolator->SetInputImage( metricImageImageConstIt.Get() );
    metric = m_Interpolator->Evaluate( point );
    size = previousSearchRegionIt.Get().GetSize();
    for( unsigned int i = 0; i < ImageDimension; ++i )
      {
      radius[i] = ( size[i] - 1 ) / 2 * m_Functor( metric );
      }
    searchRegionRadiusIt.Set( radius );
    }
}


template <class TFixedImage, class TMovingImage, class TMetricImage,
          class TDisplacementImage, class TFunctor, class TInterpolatorPrecisionType>
void
MultiResolutionSimilarityFunctionSearchRegionImageSource<
  TFixedImage, TMovingImage, TMetricImage, TDisplacementImage, TFunctor,
  TInterpolatorPrecisionType>
::BeforeThreadedGenerateData()
{
  Superclass::BeforeThreadedGenerateData();

  if( this->m_CurrentLevel != 0 )
    {
    typename OutputImageType::Pointer outputPtr = this->GetOutput();
    if( !outputPtr )
      {
      return;
      }
    m_SearchRegionRadiusResampler->SetInput( this->m_SearchRegionRadiusImage );
    m_SearchRegionRadiusResampler->SetSize( outputPtr->GetRequestedRegion().GetSize() );
    m_SearchRegionRadiusResampler->SetOutputStartIndex( outputPtr->GetRequestedRegion().GetIndex() );
    m_SearchRegionRadiusResampler->SetOutputSpacing( outputPtr->GetSpacing() );
    m_SearchRegionRadiusResampler->SetOutputOrigin( outputPtr->GetOrigin() );
    m_SearchRegionRadiusResampler->SetOutputDirection( outputPtr->GetDirection() );
    m_SearchRegionRadiusResampler->UpdateLargestPossibleRegion();
    }

  if( !m_TopLevelRadiusSpecified )
    {
    itkExceptionMacro( << "The top level radius has not been specified." );
    }
}


template <class TFixedImage, class TMovingImage, class TMetricImage,
          class TDisplacementImage, class TFunctor, class TInterpolatorPrecisionType>
void
MultiResolutionSimilarityFunctionSearchRegionImageSource<
  TFixedImage, TMovingImage, TMetricImage, TDisplacementImage, TFunctor,
  TInterpolatorPrecisionType>
::ThreadedGenerateData( const OutputRegionType & outputRegion, ThreadIdType threadID )
{
  OutputImageType * outputPtr = this->GetOutput();
  if( !outputPtr )
    {
    return;
    }

  OutputRegionType                    region;
  OutputRegionType                    movingLargestRegion = this->m_MovingImage->GetLargestPossibleRegion();
  typename MovingImageType::PointType point;
  typename MovingImageType::IndexType index;
  typename MovingImageType::SizeType  unitySize;
  unitySize.Fill( 1 );
  RadiusType   radius;

  if( this->m_CurrentLevel == 0 )
    {
    ImageRegionIteratorWithIndex<OutputImageType> it( outputPtr, outputRegion );
    for( it.GoToBegin();
         !it.IsAtEnd();
         ++it )
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
    IndexType                          startIndex = this->m_MovingImage->GetLargestPossibleRegion().GetIndex();
    IndexType                          endIndex;
    IndexType                          closestIndex;
    RadiusType                         minimumSearchRegionRadius;
    typename MovingImageType::SizeType minimumRegionSize;
    for( unsigned int i = 0; i < ImageDimension; ++i )
      {
      minimumSearchRegionRadius[i] =
        static_cast<unsigned int>( Math::Ceil( static_cast< double >( this->m_FixedBlockRadius[i] ) *
                                             this->m_MinimumSearchRegionRadiusFactor[i] ) );
      // @todo figure out what the theoretical minimum required here without
      // getting weird out-of-region bounds issues.
      minimumRegionSize[i] = 2 * minimumSearchRegionRadius[i] + 1;
      endIndex[i] = startIndex[i] + this->m_MovingImage->GetLargestPossibleRegion().GetSize()[i] - 1 - minimumRegionSize[i];
      }

    ImageRegionIteratorWithIndex<OutputImageType>   it( outputPtr, outputRegion );
    ImageRegionConstIterator<DisplacementImageType> dispIt( this->m_DisplacementResampler->GetOutput(),
                                                            outputRegion );
    ImageRegionConstIterator<SearchRegionRadiusImageType> radiusIt(
      m_SearchRegionRadiusResampler->GetOutput(), outputRegion );
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
        radius[i] = vcl_ceil( static_cast<float>( radiusIt.Get()[i] )
                              * static_cast<float>( this->m_PyramidSchedule( this->m_CurrentLevel - 1, i ) )
                              / static_cast<float>( this->m_PyramidSchedule( this->m_CurrentLevel, i ) ) );
        if( radius[i] < minimumSearchRegionRadius[i] )
          {
          radius[i] = minimumSearchRegionRadius[i];
          }
        else if( radius[i] > m_TopLevelRadius[i] )
          {
          radius[i] = m_TopLevelRadius[i] / 2;
          }
        }
      region.PadByRadius( radius );
      if( !region.Crop( movingLargestRegion ) )
        {
        // Set to the closest index and with a valid minimumsearch region
        // radius. .
        for( unsigned int i = 0; i < ImageDimension; ++i )
          {
          closestIndex[i] = index[i];
          if( index[i] < startIndex[i] )
            {
            closestIndex[i] = startIndex[i];
            }
          else if( index[i] > endIndex[i] )
            {
            closestIndex[i] = endIndex[i] - minimumRegionSize[i];
            }
          }
        region.SetIndex( closestIndex );
        region.SetSize( minimumRegionSize );
        }
      it.Set( region );
      }
    }
}

} // end namespace itk
} // end namespace BlockMatching

#endif
