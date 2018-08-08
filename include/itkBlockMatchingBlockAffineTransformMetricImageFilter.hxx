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
#ifndef itkBlockMatchingBlockAffineTransformMetricImageFilter_hxx
#define itkBlockMatchingBlockAffineTransformMetricImageFilter_hxx

#include "itkBlockMatchingBlockAffineTransformMetricImageFilter.h"

#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionIterator.h"

namespace itk
{
namespace BlockMatching
{

template< typename TFixedImage,  typename TMovingImage,
          typename TMetricImage, typename TStrainValueType >
BlockAffineTransformMetricImageFilter< TFixedImage, TMovingImage,
                                       TMetricImage, TStrainValueType >
::BlockAffineTransformMetricImageFilter()
{
  m_MetricImageFilter = nullptr;

  m_StrainImage = nullptr;

  m_Transform    = TransformType::New();
  m_Interpolator = InterpolatorType::New();
  m_StrainInterpolator = StrainInterpolatorType::New();

  m_TransformedFixedImage = FixedImageType::New();
}

template< typename TFixedImage,  typename TMovingImage,
          typename TMetricImage, typename TStrainValueType >
void
BlockAffineTransformMetricImageFilter< TFixedImage, TMovingImage,
                                       TMetricImage, TStrainValueType >
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();

  // const cast so we can set the requested region.
  FixedImageType * fixedPtr = const_cast< TFixedImage* >( this->GetInput(0) );
  if( !fixedPtr )
    {
    return;
    }

  fixedPtr->SetRequestedRegionToLargestPossibleRegion();
}

template< typename TFixedImage,  typename TMovingImage,
          typename TMetricImage, typename TStrainValueType >
void
BlockAffineTransformMetricImageFilter< TFixedImage, TMovingImage,
                                       TMetricImage, TStrainValueType >
::GenerateData()
{
  if( m_MetricImageFilter.GetPointer() == nullptr )
    {
    itkExceptionMacro(<< "The internal MetricImageFilter has not been set.");
    }

  FixedImageType * fixedPtr = const_cast< TFixedImage* >( this->GetInput(0) );
  if( !fixedPtr )
    {
    return;
    }
  MovingImageType * movingPtr = const_cast< TMovingImage* >( this->GetInput(1) );
  if( !movingPtr )
    {
    return;
    }
  if( m_StrainImage.GetPointer() != nullptr )
    {
    m_TransformedFixedImage->CopyInformation( fixedPtr );
    m_TransformedFixedImage->SetRequestedRegion( this->m_FixedImageRegion );
    m_TransformedFixedImage->SetBufferedRegion(  this->m_FixedImageRegion );
    m_TransformedFixedImage->SetLargestPossibleRegion( fixedPtr->GetLargestPossibleRegion() );
    m_TransformedFixedImage->Allocate();

    m_Transform->SetIdentity();
    typename TransformType::InputPointType center;
    typename FixedImageType::PointType     point;
    fixedPtr->TransformIndexToPhysicalPoint( this->m_FixedImageRegion.GetIndex(), point );
    const typename FixedImageType::SpacingType spacing = fixedPtr->GetSpacing();

    for(unsigned int i = 0; i < ImageDimension; ++i )
      {
      center[i] = point[i] + this->m_FixedRadius[i] * spacing[i];
      }
    m_Transform->SetCenter( center );
    m_StrainInterpolator->SetInputImage( m_StrainImage );
    ContinuousIndex< double, ImageDimension > strainIndex;
    m_StrainImage->TransformPhysicalPointToContinuousIndex( center, strainIndex );
    const typename StrainImageType::PixelType strainPixel = m_StrainInterpolator->EvaluateAtContinuousIndex( strainIndex );
    typename TransformType::OutputVectorType scaling;
    for(unsigned int i = 0; i < ImageDimension; ++i )
      {
      // Note that we use a minus here because of the direction of the transform
      // -- transformed fixed image to original fixed image.
      scaling[i] = 1.0 - strainPixel( i, i );
      }
    m_Transform->Scale( scaling );

    m_Interpolator->SetInputImage( fixedPtr );
    ImageRegionIterator< FixedImageType > it( m_TransformedFixedImage, this->m_FixedImageRegion );
    typename FixedImageType::PointType originalPoint;
    for( it.GoToBegin(); !it.IsAtEnd(); ++it )
      {
      m_TransformedFixedImage->TransformIndexToPhysicalPoint( it.GetIndex(), point );
      originalPoint = m_Transform->TransformPoint( point );
      if ( m_Interpolator->IsInsideBuffer( originalPoint ) )
        {
        it.Set( static_cast< typename FixedImageType::PixelType >(
        m_Interpolator->Evaluate( originalPoint ) ));
        }
      else
        {
        // \todo use extrapolator here
        it.Set( NumericTraits< typename FixedImageType::PixelType >::Zero );
        }
      }
    m_MetricImageFilter->SetFixedImage( m_TransformedFixedImage );
    }
  else // when no strain image is available...
    {
    m_MetricImageFilter->SetFixedImage( fixedPtr );
    }

  m_MetricImageFilter->SetMovingImage( movingPtr );
  m_MetricImageFilter->SetFixedImageRegion(  this->m_FixedImageRegion );
  m_MetricImageFilter->SetMovingImageRegion( this->m_MovingImageRegion );
  m_MetricImageFilter->GraftOutput( this->GetOutput() );
  m_MetricImageFilter->Update();
  this->GraftOutput( m_MetricImageFilter->GetOutput() );
}

} // end namespace BlockMatching
} // end namespace itk

#endif
