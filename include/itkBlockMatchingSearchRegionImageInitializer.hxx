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
#ifndef itkBlockMatchingSearchRegionImageInitializer_hxx
#define itkBlockMatchingSearchRegionImageInitializer_hxx

#include "itkBlockMatchingSearchRegionImageInitializer.h"

#include "itkImageRegionIteratorWithIndex.h"

namespace itk
{
namespace BlockMatching
{

template < typename TFixedImage, typename TMovingImage >
SearchRegionImageInitializer< TFixedImage, TMovingImage >
::SearchRegionImageInitializer():
  m_Overlap( 1.0 )
{
  m_FixedImage  = nullptr;
  m_MovingImage = nullptr;

  m_FixedBlockRadius.Fill( 0 );
  m_SearchRegionRadius.Fill( 0 );
}


template < typename TFixedImage, typename TMovingImage >
void
SearchRegionImageInitializer< TFixedImage, TMovingImage >
::SetFixedImage( FixedImageType * fixedImage )
{
  if (this->m_FixedImage.GetPointer() != fixedImage )
    {
    this->m_FixedImage = fixedImage;
    this->Modified();
    }
}


template < typename TFixedImage, typename TMovingImage >
void
SearchRegionImageInitializer< TFixedImage, TMovingImage >
::SetMovingImage( MovingImageType * movingImage )
{
  if (this->m_MovingImage.GetPointer() != movingImage )
    {
    this->m_MovingImage = movingImage;
    this->Modified();
    }
}


template < typename TFixedImage, typename TMovingImage >
void
SearchRegionImageInitializer< TFixedImage, TMovingImage >
::BeforeThreadedGenerateData()
{
  if( m_MovingImage.GetPointer() == nullptr )
    {
    itkExceptionMacro( << "Moving Image is not present." );
    }
  m_MovingImage->UpdateOutputInformation();

  RadiusType nullRadius;
  nullRadius.Fill( 0 );
  if( m_SearchRegionRadius == nullRadius )
    {
    itkExceptionMacro( << "The SearchRegionRadius has not been set." );
    }
}


template < typename TFixedImage, typename TMovingImage >
void
SearchRegionImageInitializer< TFixedImage, TMovingImage >
::GenerateOutputInformation()
{
  OutputImageType * outputPtr = this->GetOutput();
  if( !outputPtr )
    {
    return;
    }

  if( m_FixedImage.GetPointer() == nullptr )
    {
    itkExceptionMacro( << "Fixed Image is not present." );
    }
  m_FixedImage->UpdateOutputInformation();
  outputPtr->SetDirection( m_FixedImage->GetDirection() );

  // Set origin.  The first block is in the corner of the fixed image.
  typename FixedImageType::PointType  fixedOrigin = m_FixedImage->GetOrigin();
  typename OutputImageType::PointType origin;

  typename FixedImageType::SpacingType fixedSpacing = m_FixedImage->GetSpacing();

  RadiusType nullRadius;
  nullRadius.Fill( 0 );
  if( m_FixedBlockRadius == nullRadius )
    {
    itkExceptionMacro( << "The FixedBlockRadius has not been set." );
    }

  typename FixedImageType::IndexType fixedIndex = m_FixedImage->GetLargestPossibleRegion().GetIndex();
  for( unsigned int i = 0; i < ImageDimension; ++i )
    {
    origin[i] = fixedOrigin[i] +
      + ( fixedIndex[i] + m_FixedBlockRadius[i] ) * fixedSpacing[i];
    }
  outputPtr->SetOrigin( origin );

  typename OutputImageType::SpacingType spacing;
  for( unsigned int i = 0; i < ImageDimension; ++i )
    {
    spacing[i] = 2 * m_FixedBlockRadius[i] * fixedSpacing[i] * m_Overlap;
    }
  outputPtr->SetSpacing( spacing );

  OutputRegionType region;
  typename OutputRegionType::IndexType index;
  index.Fill( 0 );
  region.SetIndex( index );
  typename OutputRegionType::SizeType size;
  typename FixedImageType::SizeType fixedSize = m_FixedImage->GetLargestPossibleRegion().GetSize();
  for( unsigned int i = 0; i < ImageDimension; ++i )
    {
    size[i] = static_cast< typename OutputRegionType::SizeType::SizeValueType >( vcl_floor(
    ( fixedSize[i] - m_FixedBlockRadius[i] - 1 ) * fixedSpacing[i] / spacing[i] ));
    }
  region.SetSize( size );
  outputPtr->SetLargestPossibleRegion( region );
}


template < typename TFixedImage, typename TMovingImage >
void
SearchRegionImageInitializer< TFixedImage, TMovingImage >
::DynamicThreadedGenerateData( const OutputRegionType& outputRegion )
{
  OutputImageType * outputPtr = this->GetOutput();
  if( !outputPtr )
    {
    return;
    }

  OutputRegionType region;
  OutputRegionType movingLargestRegion = m_MovingImage->GetLargestPossibleRegion();
  typename MovingImageType::PointType point;
  typename MovingImageType::IndexType index;
  typename MovingImageType::SizeType unitySize;
  unitySize.Fill( 1 );

  ImageRegionIteratorWithIndex< OutputImageType > it( outputPtr, outputRegion );
  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    index = it.GetIndex();
    outputPtr->TransformIndexToPhysicalPoint( index, point );
    m_MovingImage->TransformPhysicalPointToIndex( point, index );
    region.SetIndex( index );
    region.SetSize( unitySize );
    region.PadByRadius( m_SearchRegionRadius );
    if( !region.Crop( movingLargestRegion ) )
      {
      itkExceptionMacro( << "The all of the fixed image must overlap with the moving image for this initializer." );
      }
    it.Set( region );
    }
}

} // end namespace BlockMatching
} // end namespace itk

#endif
