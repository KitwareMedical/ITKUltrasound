#ifndef itkBlockMatchingMetricImageFilter_hxx
#define itkBlockMatchingMetricImageFilter_hxx

#include "itkBlockMatchingMetricImageFilter.h"

namespace itk
{
namespace BlockMatching
{

template <class TFixedImage, class TMovingImage,
          class TMetricImage >
MetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
::MetricImageFilter() :
  m_FixedImageRegionDefined( false ),
  m_MovingImageRegionDefined( false ),
  m_MinimumSplitSize( 6 )
{
}

template <class TFixedImage, class TMovingImage,
          class TMetricImage >
void
MetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
::SetFixedImage( FixedImageType * fixedImage )
{
  this->SetInput( 0, fixedImage );
}

template <class TFixedImage, class TMovingImage,
          class TMetricImage >
void
MetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
::SetMovingImage( MovingImageType * movingImage )
{
  this->SetInput( 1, movingImage );
}

template <class TFixedImage, class TMovingImage,
          class TMetricImage >
void
MetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
::SetFixedImageRegion( const FixedImageRegionType & region )
{
  typename FixedImageType::Pointer fixedPtr = const_cast< TFixedImage* >( this->GetInput(0) );
  if( !fixedPtr )
    {
    itkExceptionMacro( << "The FixedImage must be set before specifying the fixed image region." );
    }
  fixedPtr->UpdateOutputInformation();
  m_FixedImageRegion = region;
  if( !m_FixedImageRegion.Crop( fixedPtr->GetLargestPossibleRegion() ) )
    {
    itkExceptionMacro( << "Requested block is outside of the fixed image." );
    }
  typename FixedImageRegionType::SizeType fixedSize     = m_FixedImageRegion.GetSize();
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    // The radius may have been truncated if the fixed region was outside the
    // fixed image's LargestPossibleRegion.
    if( fixedSize[i] % 2 == 0 )
      fixedSize[i]--;
    m_FixedRadius[i] = (fixedSize[i] - 1) / 2;
    }
  m_FixedImageRegion.SetSize( fixedSize );
  m_FixedImageRegionDefined = true;

  typename MovingImageType::Pointer movingPtr = const_cast< TMovingImage* >( this->GetInput(1) );
  if( !movingPtr )
    {
    itkExceptionMacro( << "The MovingImage must be set before specifying the fixed image region." );
    }
  movingPtr->UpdateOutputInformation();
  m_MovingRadius = m_FixedRadius;
  typename FixedImageType::SpacingType  fixedSpacing  = fixedPtr->GetSpacing();
  typename MovingImageType::SpacingType movingSpacing = movingPtr->GetSpacing();
  if( !( fixedSpacing == movingSpacing ) )
    {
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      m_MovingRadius[i] = static_cast< typename RadiusType::SizeValueType >( vcl_ceil((
      fixedSpacing[i] * m_FixedRadius[i] / movingSpacing[i] )));
      }
    }
  this->Modified();
}

template <class TFixedImage, class TMovingImage,
          class TMetricImage >
void
MetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
::SetMovingImageRegion( const MovingImageRegionType & region )
{
  m_MovingImageRegion = region;
  m_MovingImageRegionDefined = true;
  this->Modified();
}

template <class TFixedImage, class TMovingImage,
          class TMetricImage >
void
MetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
::GenerateOutputInformation()
{
  typename MovingImageType::ConstPointer movingPtr = this->GetInput( 1 );
  if( !movingPtr )
    return;

  typename MetricImageType::Pointer outputPtr = this->GetOutput();
  if( !outputPtr )
    return;

  if( !m_MovingImageRegionDefined )
    {
    itkExceptionMacro( << "MovingImageRegion has not been set" );
    }

  MetricImageRegionType                     metricRegion;
  typename MetricImageRegionType::IndexType metricIndex;
  metricIndex.Fill( 0 );
  metricRegion.SetIndex( metricIndex );
  // the Default is to to use the moving image size and spacing.
  metricRegion.SetSize( m_MovingImageRegion.GetSize() );
  outputPtr->SetLargestPossibleRegion( metricRegion );
  outputPtr->SetSpacing( movingPtr->GetSpacing() );

  typename MetricImageType::IndexType metricStart( m_MovingImageRegion.GetIndex() );

  typename MetricImageType::PointType origin;
  movingPtr->TransformIndexToPhysicalPoint( metricStart, origin );
  outputPtr->SetOrigin( origin );

  // The metric image direction is the same as the fixed image direction.
  outputPtr->SetDirection( movingPtr->GetDirection() );
}

template <class TFixedImage, class TMovingImage,
          class TMetricImage >
const int &
MetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
::GetNumberOfThreads()
{
  typename MetricImageType::Pointer outputPtr = this->GetOutput();

  if( !outputPtr )
    return Superclass::GetNumberOfThreads();

  typename MetricImageRegionType::SizeType requestedRegionSize = outputPtr->GetRequestedRegion().GetSize();
  // split on the outermost dimension available
  int splitAxis = outputPtr->GetImageDimension() - 1;
  while ( requestedRegionSize[splitAxis] == 1 )
    {
    --splitAxis;
    if ( splitAxis < 0 )
      { // cannot split
      itkDebugMacro("  Cannot Split");
      return Superclass::GetNumberOfThreads();
      }
    }

  // determine the actual number of pieces that will be generated
  typename MetricImageRegionType::SizeType::SizeValueType range = requestedRegionSize[splitAxis];
  int valuesPerThread = Math::Ceil< int >(range / (double)Superclass::GetNumberOfThreads());
  if( valuesPerThread < m_MinimumSplitSize )
    {
    m_SpecialThreadCount = Math::Floor< int >( range / (double)m_MinimumSplitSize );
    return m_SpecialThreadCount;
    }

  return Superclass::GetNumberOfThreads();
}

template <class TFixedImage, class TMovingImage,
          class TMetricImage >
void
MetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();

  // const cast so we can set the requested region.
  typename FixedImageType::Pointer fixedPtr = const_cast< TFixedImage* >( this->GetInput(0) );
  if( !fixedPtr )
    {
    return;
    }

  typename MovingImageType::Pointer movingPtr = const_cast< TMovingImage* >( this->GetInput(1) );
  if( !movingPtr )
    {
    return;
    }

  if( !m_FixedImageRegionDefined )
    {
    itkExceptionMacro( << "FixedImageRegion has not been set" );
    }

  if( !m_MovingImageRegionDefined )
    {
    itkExceptionMacro( << "MovingImageRegion has not been set" );
    }

  fixedPtr->SetRequestedRegion( m_FixedImageRegion );
  MovingImageRegionType movingImageRequestedRegion = m_MovingImageRegion;
  movingImageRequestedRegion.PadByRadius( m_MovingRadius );
  // make sure the requested region is within the largest possible.
  if ( movingImageRequestedRegion.Crop( movingPtr->GetLargestPossibleRegion() ) )
    {
    movingPtr->SetRequestedRegion( movingImageRequestedRegion );
    return;
    }
  else
    {
    // store what we tried( prior to try to crop )
    movingPtr->SetRequestedRegion( movingImageRequestedRegion );

    itkExceptionMacro(<< "Moving image requested region is at least partially outside the LargestPossibleRegion.");
    }
}

template <class TFixedImage, class TMovingImage,
          class TMetricImage >
void
MetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
::EnlargeOutputRequestedRegion(DataObject *data)
{
  Superclass::EnlargeOutputRequestedRegion( data );
  data->SetRequestedRegionToLargestPossibleRegion();
}

} // end namespace BlockMatching
} // end namespace itk

#endif
