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
#ifndef itkBlockMatchingImageToImageMetricMetricImageFilter_hxx
#define itkBlockMatchingImageToImageMetricMetricImageFilter_hxx

#include "itkBlockMatchingImageToImageMetricMetricImageFilter.h"

namespace itk
{
namespace BlockMatching
{

template< typename TFixedImage, typename TMovingImage, typename TMetricImage >
ImageToImageMetricMetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
::ImageToImageMetricMetricImageFilter():
  m_MetricImageSpacingDefined( false )
{
}


template< typename TFixedImage, typename TMovingImage, typename TMetricImage >
void
ImageToImageMetricMetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
::SetMetricImageSpacing( const MetricImageSpacingType & spacing )
{
  m_MetricImageSpacing = spacing;
  m_MetricImageSpacingDefined = true;
  this->Modified();
}


template< typename TFixedImage, typename TMovingImage, typename TMetricImage >
void
ImageToImageMetricMetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
::GenerateOutputInformation()
{
  // get origin and direction from fixed image.
  Superclass::Superclass::GenerateOutputInformation();

  const MovingImageType * moving = this->GetInput(1);
  if( !moving )
    {
    itkExceptionMacro( << "MovingImage input has not been set" );
    }

  MetricImageType * output = this->GetOutput(0);
  if( !output )
    {
    return;
    }

  if( !this->m_MovingImageRegionDefined )
    {
    itkExceptionMacro( << "MovingImageRegion has not been set" );
    }

  typename MovingImageType::SpacingType movingSpacing = moving->GetSpacing();

  MetricImageRegionType metricRegion;
  typename MetricImageRegionType::IndexType metricIndex;
  metricIndex.Fill( 0 );
  metricRegion.SetIndex( metricIndex );
  typename MetricImageRegionType::SizeType  metricSize;

  typename MovingImageRegionType::SizeType movingSize = this->m_MovingImageRegion.GetSize();

  if( m_MetricImageSpacingDefined )
    {
    for( unsigned int i = 0; i < ImageDimension; ++i )
      {
      metricSize[i] = std::ceil( movingSize[i] * movingSpacing[i] / m_MetricImageSpacing[i] );
      }
    output->SetSpacing( m_MetricImageSpacing );
    }
  else
    {
    metricSize = movingSize;
    output->SetSpacing( movingSpacing );
    }
  metricRegion.SetSize( metricSize );
  output->SetLargestPossibleRegion( metricRegion );
}

} // end namespace BlockMatching
} // end namespace itk

#endif
