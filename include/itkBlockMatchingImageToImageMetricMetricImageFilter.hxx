#ifndef itkBlockMatchingImageToImageMetricMetricImageFilter_hxx
#define itkBlockMatchingImageToImageMetricMetricImageFilter_hxx

#include "itkBlockMatchingImageToImageMetricMetricImageFilter.h"

namespace itk
{
namespace BlockMatching
{

template <class TFixedImage, class TMovingImage,
	 class TMetricImage >
ImageToImageMetricMetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
::ImageToImageMetricMetricImageFilter():
  m_MetricImageSpacingDefined( false )
{
}


template <class TFixedImage, class TMovingImage,
	 class TMetricImage >
void
ImageToImageMetricMetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
::SetMetricImageSpacing( const MetricImageSpacingType & spacing )
{
  m_MetricImageSpacing = spacing;
  m_MetricImageSpacingDefined = true;
  this->Modified();
}


template <class TFixedImage, class TMovingImage,
	 class TMetricImage >
void
ImageToImageMetricMetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
::GenerateOutputInformation()
{
  // get origin and direction from fixed image.
  Superclass::Superclass::GenerateOutputInformation();

  typename MovingImageType::ConstPointer movingPtr = this->GetInput(1);
  if( !movingPtr )
    {
    itkExceptionMacro( << "MovingImage input has not been set" );
    }

  typename MetricImageType::Pointer outputPtr = this->GetOutput(0);
  if( !outputPtr )
    {
    return;
    }

  if( !this->m_MovingImageRegionDefined )
    {
    itkExceptionMacro( << "MovingImageRegion has not been set" );
    }

  typename MovingImageType::SpacingType movingSpacing = movingPtr->GetSpacing();

  MetricImageRegionType metricRegion;
  typename MetricImageRegionType::IndexType metricIndex;
  metricIndex.Fill( 0 );
  metricRegion.SetIndex( metricIndex );
  typename MetricImageRegionType::SizeType  metricSize;

  typename MovingImageRegionType::SizeType movingSize = this->m_MovingImageRegion.GetSize();

  if( m_MetricImageSpacingDefined )
    {
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      metricSize[i] = vcl_ceil( movingSize[i] * movingSpacing[i] / m_MetricImageSpacing[i] );
      }
    outputPtr->SetSpacing( m_MetricImageSpacing );
    }
  else
    {
    metricSize = movingSize;
    outputPtr->SetSpacing( movingSpacing );
    }
  metricRegion.SetSize( metricSize );
  outputPtr->SetLargestPossibleRegion( metricRegion );
}

} // end namespace BlockMatching
} // end namespace itk

#endif

