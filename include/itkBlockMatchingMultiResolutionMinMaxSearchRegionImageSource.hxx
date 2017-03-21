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
#ifndef itkBlockMatchingMultiResolutionMinMaxSearchRegionImageSource_hxx
#define itkBlockMatchingMultiResolutionMinMaxSearchRegionImageSource_hxx

#include "itkBlockMatchingMultiResolutionMinMaxSearchRegionImageSource.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"

namespace itk
{
namespace BlockMatching
{

template < class TFixedImage, class TMovingImage, class TDisplacementImage >
void
MultiResolutionMinMaxSearchRegionImageSource< TFixedImage, TMovingImage, TDisplacementImage >
::ThreadedGenerateData( const OutputRegionType& outputRegion, ThreadIdType threadID )
{
  OutputImageType * outputPtr = this->GetOutput();
  if( !outputPtr )
    {
    return;
    }

  typedef typename MovingImageType::IndexType IndexType;

  OutputRegionType                    region;
  OutputRegionType                    movingLargestRegion = this->m_MovingImage->GetLargestPossibleRegion();
  typename MovingImageType::PointType point;
  IndexType                           index;
  typename MovingImageType::SizeType  unitySize;
  typename MovingImageType::SizeType  minimumRegionSize;
  unitySize.Fill( 1 );
  IndexType startIndex = this->m_MovingImage->GetLargestPossibleRegion().GetIndex();
  IndexType endIndex;
  IndexType closestIndex;

  double slope;
  RadiusType radius;
  for( unsigned int i = 0; i < ImageDimension; ++i )
    {
    slope = (this->m_MinFactor[i] - this->m_MaxFactor[i]) / (this->m_PyramidSchedule.rows() - 1.0);
    radius[i] = Math::Ceil< typename RadiusType::SizeValueType >( this->m_FixedBlockRadius[i] *
                                                        ( slope * this->m_CurrentLevel + this->m_MaxFactor[i] ) );
    minimumRegionSize[i] = 2 * radius[i] + 1;
    endIndex[i] = startIndex[i] + this->m_MovingImage->GetLargestPossibleRegion().GetSize()[i] - 1 -
      minimumRegionSize[i];
    }

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
      region.PadByRadius( radius );
      if( !region.Crop( movingLargestRegion ) )
        {
        itkExceptionMacro( << "Attempted to create a search region entirely outside the moving image." );
        }
      it.Set( region );
      }
    }
  else
    {
    ImageRegionIteratorWithIndex< OutputImageType > it( outputPtr, outputRegion );
    ImageRegionConstIterator< DisplacementImageType > dispIt( this->m_DisplacementResampler->GetOutput(),
                                                              outputRegion );

    for( it.GoToBegin(), dispIt.GoToBegin();
         !it.IsAtEnd();
         ++it,
         ++dispIt )
      {
      index = it.GetIndex();
      outputPtr->TransformIndexToPhysicalPoint( index, point );
      // resample displacement image
      this->m_MovingImage->TransformPhysicalPointToIndex( point + dispIt.Get(),
                                                          index );
      region.SetIndex( index );
      region.SetSize( unitySize );
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
