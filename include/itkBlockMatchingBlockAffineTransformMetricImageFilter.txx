#ifndef __itkBlockMatchingBlockAffineTransformMetricImageFilter_txx
#define __itkBlockMatchingBlockAffineTransformMetricImageFilter_txx

#include "itkBlockMatchingBlockAffineTransformMetricImageFilter.h"

#include "itkImageRegionIteratorWithIndex.h"

namespace itk
{
namespace BlockMatching
{

template <class TFixedImage,  class TMovingImage,
          class TMetricImage, class TStrainValueType >
BlockAffineTransformMetricImageFilter< TFixedImage, TMovingImage,
                                       TMetricImage, TStrainValueType >
::BlockAffineTransformMetricImageFilter()
{
  m_MetricImageFilter = NULL;

  m_StrainImage = NULL;

  m_Transform    = TransformType::New();
  m_Interpolator = InterpolatorType::New();
  m_StrainInterpolator = StrainInterpolatorType::New();

  m_TransformedFixedImage = FixedImageType::New();
}

template <class TFixedImage,  class TMovingImage,
          class TMetricImage, class TStrainValueType >
void
BlockAffineTransformMetricImageFilter< TFixedImage, TMovingImage,
                                       TMetricImage, TStrainValueType >
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();

  // const cast so we can set the requested region.
  typename FixedImageType::Pointer fixedPtr = const_cast< TFixedImage* >( this->GetInput(0) );
  if( !fixedPtr )
    {
    return;
    }

  fixedPtr->SetRequestedRegionToLargestPossibleRegion();
}

template <class TFixedImage,  class TMovingImage,
          class TMetricImage, class TStrainValueType >
void
BlockAffineTransformMetricImageFilter< TFixedImage, TMovingImage,
                                       TMetricImage, TStrainValueType >
::GenerateData()
{
  if( m_MetricImageFilter.GetPointer() == NULL )
    {
    itkExceptionMacro(<< "The internal MetricImageFilter has not been set.");
    }

  typename FixedImageType::Pointer fixedPtr = const_cast< TFixedImage* >( this->GetInput(0) );
  if( !fixedPtr.GetPointer() )
    {
    return;
    }
  typename MovingImageType::Pointer movingPtr = const_cast< TMovingImage* >( this->GetInput(1) );
  if( !movingPtr.GetPointer() )
    {
    return;
    }
  if( m_StrainImage.GetPointer() != NULL )
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

    unsigned int i;
    for( i = 0; i < ImageDimension; ++i )
      {
      center[i] = point[i] + this->m_FixedRadius[i] * spacing[i];
      }
    m_Transform->SetCenter( center );
    m_StrainInterpolator->SetInputImage( m_StrainImage );
    ContinuousIndex< double, ImageDimension > strainIndex;
    m_StrainImage->TransformPhysicalPointToContinuousIndex( center, strainIndex );
    const typename StrainImageType::PixelType strainPixel = m_StrainInterpolator->EvaluateAtContinuousIndex( strainIndex );
    typename TransformType::OutputVectorType scaling;
    for( i = 0; i < ImageDimension; ++i )
      {
      // Note that we use a minus here because of the direction of the transform
      // -- transformed fixed image to original fixed image.
      scaling[i] = 1.0 - strainPixel( i, i );
      }
    m_Transform->Scale( scaling );

    m_Interpolator->SetInputImage( fixedPtr.GetPointer() );
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
