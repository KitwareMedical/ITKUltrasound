#ifndef __itkBlockMatchingMultiResolutionMinMaxSearchRegionImageSource_hxx
#define __itkBlockMatchingMultiResolutionMinMaxSearchRegionImageSource_hxx

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
::ThreadedGenerateData( const OutputRegionType& outputRegion, int threadID )
{
  typename OutputImageType::Pointer outputPtr = this->GetOutput();
  if( !outputPtr )
    return;

  typedef typename MovingImageType::IndexType IndexType;

  OutputRegionType                    region;
  OutputRegionType                    movingLargestRegion = this->m_MovingImage->GetLargestPossibleRegion();
  typename MovingImageType::PointType point;
  IndexType                           index;
  typename MovingImageType::SizeType  unitySize;
  typename MovingImageType::SizeType  minimumRegionSize;
  unitySize.Fill( 1 );
  unsigned int i;
  IndexType startIndex = this->m_MovingImage->GetLargestPossibleRegion().GetIndex();
  IndexType endIndex;
  IndexType closestIndex;

  double slope;
  RadiusType radius;
  for( i = 0; i < ImageDimension; ++i )
    {
    slope = (this->m_MinFactor[i] - this->m_MaxFactor[i]) / (this->m_PyramidSchedule.rows() - 1.0);
    radius[i] =
      Math::Ceil< typename RadiusType::SizeValueType >( this->m_FixedBlockRadius[i] *
                                                        ( slope * this->m_CurrentLevel + this->m_MaxFactor[i] ) );
    minimumRegionSize[i] = 2 * radius[i] + 1;
    endIndex[i] = startIndex[i] + this->m_MovingImage->GetLargestPossibleRegion().GetSize()[i] - 1 -
      minimumRegionSize[i];
    }

  if( this->m_CurrentLevel == 0 )
    {
    ::itk::ImageRegionIteratorWithIndex< OutputImageType > it( outputPtr, outputRegion );
    for( it.GoToBegin();
         !it.IsAtEnd();
         ++it )
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
    ::itk::ImageRegionIteratorWithIndex< OutputImageType > it( outputPtr, outputRegion );
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
        for( i = 0; i < ImageDimension; ++i )
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
