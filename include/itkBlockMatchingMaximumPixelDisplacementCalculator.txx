#ifndef __itkBlockMatchingMaximumPixelDisplacementCalculator_txx
#define __itkBlockMatchingMaximumPixelDisplacementCalculator_txx

#include "itkBlockMatchingMaximumPixelDisplacementCalculator.h"

#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk
{
namespace BlockMatching
{

template < class TMetricImage, class TDisplacementImage >
void
MaximumPixelDisplacementCalculator< TMetricImage, TDisplacementImage >
::SetMetricImagePixel( const PointType& point, const IndexType& index,
  MetricImageType* metricImage )
{
  Superclass::SetMetricImagePixel( point, index, metricImage );

  PixelType max = NumericTraits< PixelType >::min();
  IndexType maxIndex;
  maxIndex.Fill( 0 );

  itk::ImageRegionConstIteratorWithIndex< MetricImageType >
    it( metricImage, metricImage->GetBufferedRegion() );
  for( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    if( it.Get() > max )
      {
      max = it.Get();
      maxIndex = it.GetIndex();
      }
    }

  PointType maxPoint;
  metricImage->TransformIndexToPhysicalPoint( maxIndex, maxPoint );
  this->m_DisplacementImage->SetPixel( index, maxPoint - point );
}

} // end namespace BlockMatching
} // end namespace itk

#endif
