#ifndef __itkBlockMatchingOptimizingInterpolationDisplacementCalculator_hxx
#define __itkBlockMatchingOptimizingInterpolationDisplacementCalculator_hxx

#include "itkBlockMatchingOptimizingInterpolationDisplacementCalculator.h"

#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk
{
namespace BlockMatching
{

template < class TMetricImage, class TDisplacementImage, class TCoordRep >
OptimizingInterpolationDisplacementCalculator< TMetricImage,
  TDisplacementImage, TCoordRep >
::OptimizingInterpolationDisplacementCalculator()
{
  m_CostFunction = OptimizingInterpolationCostFunction::New();

  // @todo sensible default optimizer and interpolator
}

template < class TMetricImage, class TDisplacementImage, class TCoordRep >
void
OptimizingInterpolationDisplacementCalculator< TMetricImage,
  TDisplacementImage, TCoordRep >
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

  unsigned int i;

  // If the maxIndex is on the edge of the image we don't try interpolation.
  bool onEdge = false;
  IndexType metricIndex = region.GetIndex();
  const typename MetricImageType::SizeType metricSize   = region.GetSize();
  for( i = 0; i < ImageDimension; ++i )
    {
    if( maxIndex[i] == metricIndex[i] || maxIndex[i] == metricIndex[i] + static_cast< typename IndexType::IndexValueType >( metricSize[i] ) )
      {
      onEdge = true;
      break;
      }
    }

  PointType maxPoint;
  if( onEdge )
    {
    metricImage->TransformIndexToPhysicalPoint( maxIndex, maxPoint );
    }
  else
    {
    typename OptimizingInterpolationCostFunction::ParametersType parameters( ImageDimension );
    for( i = 0; i < ImageDimension; i++ )
      {
      parameters[i] = static_cast< typename OptimizingInterpolationCostFunction::ParametersValueType >( maxIndex[i] );
      }
    m_Optimizer->SetInitialPosition( parameters );
    // Is this the right offset?
    for( i = 0; i < ImageDimension; i++ )
      parameters[i] += 0.1;
    m_CostFunction->Initialize( parameters );
    m_CostFunction->GetInterpolator()->SetInputImage( metricImage );
    m_Optimizer->StartOptimization();

    parameters = m_Optimizer->GetCurrentPosition();

    ContinuousIndexType continuousIndex;
    for( i = 0; i < ImageDimension; i++ )
      {
      continuousIndex[i] = parameters[i];
      }
    metricImage->TransformContinuousIndexToPhysicalPoint( continuousIndex, maxPoint );
    }

  this->m_DisplacementImage->SetPixel( displacementIndex, maxPoint - centerPoint );
}


template < class TMetricImage, class TDisplacementImage, class TCoordRep >
OptimizingInterpolationDisplacementCalculator< TMetricImage,
  TDisplacementImage, TCoordRep >
::OptimizingInterpolationCostFunction
::OptimizingInterpolationCostFunction()
{
  // This is needed to compile :S
  this->m_ContinuousIndexPtr = this->m_ContinuousIndex.GetDataPointer();
}

template < class TMetricImage, class TDisplacementImage, class TCoordRep >
typename OptimizingInterpolationDisplacementCalculator< TMetricImage, TDisplacementImage,
  TCoordRep >::OptimizingInterpolationCostFunction::MeasureType
OptimizingInterpolationDisplacementCalculator< TMetricImage,
  TDisplacementImage, TCoordRep >
::OptimizingInterpolationCostFunction::GetValue( const ParametersType& parameters ) const
{
  for( unsigned int i=0; i < ImageDimension; i++ )
    {
    m_ContinuousIndexPtr[i] = parameters[i];
    }

  MeasureType measure = -1 * m_Interpolator->EvaluateAtContinuousIndex( m_ContinuousIndex );
  return measure;
}


template < class TMetricImage, class TDisplacementImage, class TCoordRep >
void
OptimizingInterpolationDisplacementCalculator< TMetricImage,
  TDisplacementImage, TCoordRep >
::OptimizingInterpolationCostFunction::GetDerivative(
  const ParametersType& parameters,
  DerivativeType& derivative ) const
{
  MeasureType value;
  MeasureType centerValue;
  TCoordRep delta;
  unsigned int i;
  for( i = 0; i < ImageDimension; i++ )
    {
    m_ContinuousIndexPtr[i] = parameters[i];
    }
  centerValue = -1 * m_Interpolator->EvaluateAtContinuousIndex( m_ContinuousIndex );
  for( i = 0; i < ImageDimension; i++ )
    {
    const static double acceptibleValue = 100 * NumericTraits< double >::epsilon();
    delta = parameters[i] - m_PreviousParameters[i];
    if( vcl_abs( delta ) > acceptibleValue )
      {
      m_ContinuousIndexPtr[i] += delta;
      value = -1 * m_Interpolator->EvaluateAtContinuousIndex( m_ContinuousIndex );
      derivative[i] = ( value - centerValue ) / delta;
      }
    else
      {
      derivative[i] = 0.0;
      }
    m_ContinuousIndexPtr[i] = parameters[i];
    }

  for( i = 0; i < ImageDimension; ++i )
    {
    m_PreviousParametersPtr[i] = parameters[i];
    }
}

} // end namespace BlockMatching
} // end namespace itk

#endif
