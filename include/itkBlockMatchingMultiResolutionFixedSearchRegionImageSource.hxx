#ifndef itkBlockMatchingMultiResolutionFixedSearchRegionImageSource_hxx
#define itkBlockMatchingMultiResolutionFixedSearchRegionImageSource_hxx

#include "itkBlockMatchingMultiResolutionFixedSearchRegionImageSource.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkResampleImageFilter.h"

namespace itk
{
namespace BlockMatching
{

template < class TFixedImage, class TMovingImage, class TDisplacementImage >
MultiResolutionFixedSearchRegionImageSource< TFixedImage, TMovingImage,
                                             TDisplacementImage >
::MultiResolutionFixedSearchRegionImageSource():
  m_SearchRegionRadiusSet( false )
{
}

template < class TFixedImage, class TMovingImage, class TDisplacementImage >
void
MultiResolutionFixedSearchRegionImageSource< TFixedImage, TMovingImage,
                                             TDisplacementImage >
::SetSearchRegionRadiusSchedule( const unsigned int rad )
{
  RadiusType radius;
  radius.Fill( rad );
  this->SetSearchRegionRadiusSchedule( radius );
}

template < class TFixedImage, class TMovingImage, class TDisplacementImage >
void
MultiResolutionFixedSearchRegionImageSource< TFixedImage, TMovingImage,
                                             TDisplacementImage >
::SetSearchRegionRadiusSchedule( const RadiusType& radius )
{
  // Check to make sure the PyramidSchedule has been set.
  if( this->m_PyramidSchedule.size() == 0 )
    {
    itkExceptionMacro(<<"The PyramidSchedule must be set before calling this method.");
    }

  m_SearchRegionRadiusSchedule.SetSize( this->m_PyramidSchedule.rows(), this->m_PyramidSchedule.cols() );
  for( unsigned int i = 0; i < this->m_PyramidSchedule.rows(); ++i )
    {
    for( unsigned int j = 0; j < this->m_PyramidSchedule.cols(); ++j )
      {
      m_SearchRegionRadiusSchedule( i, j ) = radius[j];
      }
    }
  this->Modified();
  m_SearchRegionRadiusSet = true;
}

template < class TFixedImage, class TMovingImage, class TDisplacementImage >
void
MultiResolutionFixedSearchRegionImageSource< TFixedImage, TMovingImage,
                                             TDisplacementImage >
::BeforeThreadedGenerateData()
{
  Superclass::BeforeThreadedGenerateData();

  if( !m_SearchRegionRadiusSet )
    {
    itkExceptionMacro( << "The SearchRegionRadius has not been set." );
    }
}

template < class TFixedImage, class TMovingImage, class TDisplacementImage >
void
MultiResolutionFixedSearchRegionImageSource< TFixedImage, TMovingImage,
                                             TDisplacementImage >
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
  IndexType    startIndex = this->m_MovingImage->GetLargestPossibleRegion().GetIndex();
  IndexType    endIndex;
  IndexType    closestIndex;

  RadiusType radius;
  for( i = 0; i < this->m_SearchRegionRadiusSchedule.cols(); ++i )
    {
    radius[i] = m_SearchRegionRadiusSchedule( this->m_CurrentLevel, i );
    minimumRegionSize[i] = 2 * radius[i] + 1;
    endIndex[i] = startIndex[i] + this->m_MovingImage->GetLargestPossibleRegion().GetSize()[i] - 1 - minimumRegionSize[i];
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
