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
#ifndef itkBlockMatchingMetricImageToDisplacementCalculator_hxx
#define itkBlockMatchingMetricImageToDisplacementCalculator_hxx

#include "itkBlockMatchingMetricImageToDisplacementCalculator.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

namespace itk
{
namespace BlockMatching
{

template < typename TMetricImage, typename TDisplacementImage >
MetricImageToDisplacementCalculator< TMetricImage, TDisplacementImage >
::MetricImageToDisplacementCalculator():
  m_CacheMetricImage( false ),
  m_RegionsDefined( false )
{
  m_MetricImageImage  = nullptr;
  m_DisplacementImage = nullptr;
  m_MetricImageDuplicator = MetricImageDuplicatorType::New();
  m_MultiThreader = MultiThreaderBase::New();
}


template < typename TMetricImage, typename TDisplacementImage >
void
MetricImageToDisplacementCalculator< TMetricImage, TDisplacementImage >
::SetMetricImagePixel( const PointType& point, const IndexType& index,
                       MetricImageType * metricImage )
{
  if( m_CacheMetricImage )
    {
    // Copy the image contents to a new image.
    m_MetricImageDuplicator->SetInputImage( metricImage );
    m_MetricImageDuplicator->Update();
    m_MetricImageImage->SetPixel ( index, m_MetricImageDuplicator->GetOutput() );
    m_CenterPointsImage->SetPixel( index, point );
    }
}


template < typename TMetricImage, typename TDisplacementImage >
void
MetricImageToDisplacementCalculator< TMetricImage, TDisplacementImage >
::SetDisplacementImage( DisplacementImageType * image )
{
  if( this->m_DisplacementImage.GetPointer() != image )
    {
    this->m_DisplacementImage = image;
    this->Modified();
    }

  if( m_CacheMetricImage && ( !m_MetricImageImage.GetPointer() ||
      m_MetricImageImage->GetLargestPossibleRegion() != image->GetLargestPossibleRegion() ))
    {
    // We create a new image and reallocate here every time because if we try to
    // keep the same image reallocation appears to do something funny with the
    // SmartPointer pixels.
    m_MetricImageImage = MetricImageImageType::New();
    m_CenterPointsImage = CenterPointsImageType::New();
    m_MetricImageImage->CopyInformation( image );
    m_MetricImageImage->SetRegions( image->GetLargestPossibleRegion() );
    m_MetricImageImage->Allocate();
    m_CenterPointsImage->CopyInformation( image );
    m_CenterPointsImage->SetRegions( image->GetLargestPossibleRegion() );
    m_CenterPointsImage->Allocate();
    this->Modified();
    }
}

} // end namespace BlockMatching
} // end namespace itk

#endif
