#ifndef itkBlockMatchingCosineInterpolationDisplacementCalculator_hxx
#define itkBlockMatchingCosineInterpolationDisplacementCalculator_hxx

#include "itkBlockMatchingCosineInterpolationDisplacementCalculator.h"

#include "itkImageRegionConstIteratorWithIndex.h"

#include "vcl_cmath.h"

namespace itk
{
namespace BlockMatching
{

template < class TMetricImage, class TDisplacementImage, class TCoordRep >
CosineInterpolationDisplacementCalculator< TMetricImage, TDisplacementImage,
  TCoordRep >
::CosineInterpolationDisplacementCalculator()
{
}

template < class TMetricImage, class TDisplacementImage, class TCoordRep >
void
CosineInterpolationDisplacementCalculator< TMetricImage, TDisplacementImage,
  TCoordRep >
::SetMetricImagePixel( const PointType& centerPoint,
  const IndexType& displacementIndex,
  MetricImageType* metricImage )
{
  Superclass::SetMetricImagePixel( centerPoint, displacementIndex, metricImage );

  // Find index of the maximum value.
  PixelType max = NumericTraits< PixelType >::min();
  IndexType maxIndex;
  maxIndex.Fill( 0 );

  const typename MetricImageType::RegionType region = metricImage->GetBufferedRegion();
  itk::ImageRegionConstIteratorWithIndex< MetricImageType >
    it( metricImage, region );
  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    if( it.Get() > max )
      {
      max = it.Get();
      maxIndex = it.GetIndex();
      }
    }

  TCoordRep y1 = max;
  TCoordRep y0;
  TCoordRep y2;
  TCoordRep omega;
  TCoordRep theta;

  PointType maxPoint;
  metricImage->TransformIndexToPhysicalPoint( maxIndex, maxPoint );

  SpacingType spacing = metricImage->GetSpacing();
  IndexType tempIndex = maxIndex;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    tempIndex[i] = maxIndex[i] - 1;
    if( ! region.IsInside( tempIndex ) )
      {
      tempIndex[i] = maxIndex[i];
      continue;
      }
    y0 = metricImage->GetPixel( tempIndex );
    tempIndex[i] = maxIndex[i] + 1;
    if( ! region.IsInside( tempIndex ) )
      {
      tempIndex[i] = maxIndex[i];
      continue;
      }
    y2 = metricImage->GetPixel( tempIndex );
    tempIndex[i] = maxIndex[i];

    omega = vcl_acos( ( y0 + y2 ) / ( 2 * y1 ) );
    theta = vcl_atan( ( y0 - y2 ) / ( 2 * y1 * vcl_sin( omega ) ) ); 
    // @todo is this right ?
    maxPoint[i] += spacing[i] / itk::Math::pi * -1 * theta / omega;
    }

  this->m_DisplacementImage->SetPixel( displacementIndex, maxPoint - centerPoint );
}

} // end namespace BlockMatching
} // end namespace itk

#endif
