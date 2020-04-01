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
#ifndef itkBlockMatchingOptimizingInterpolationDisplacementCalculator_hxx
#define itkBlockMatchingOptimizingInterpolationDisplacementCalculator_hxx

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

  // If the maxIndex is on the edge of the image we don't try interpolation.
  bool onEdge = false;
  IndexType metricIndex = region.GetIndex();
  const typename MetricImageType::SizeType metricSize   = region.GetSize();
  for( unsigned int i = 0; i < ImageDimension; ++i )
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
    for( unsigned int i = 0; i < ImageDimension; ++i )
      {
      parameters[i] = static_cast< typename OptimizingInterpolationCostFunction::ParametersValueType >( maxIndex[i] );
      }
    m_Optimizer->SetInitialPosition( parameters );
    // Is this the right offset?
    for( unsigned int i = 0; i < ImageDimension; i++ )
      parameters[i] += 0.1;
    m_CostFunction->Initialize( parameters );
    m_CostFunction->GetInterpolator()->SetInputImage( metricImage );
    m_Optimizer->StartOptimization();

    parameters = m_Optimizer->GetCurrentPosition();

    ContinuousIndexType continuousIndex;
    for( unsigned int i = 0; i < ImageDimension; i++ )
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
  for( unsigned int i = 0; i < ImageDimension; ++i )
    {
    m_ContinuousIndexPtr[i] = parameters[i];
    }
  centerValue = -1 * m_Interpolator->EvaluateAtContinuousIndex( m_ContinuousIndex );
  for( unsigned int i = 0; i < ImageDimension; ++i )
    {
    const static double acceptibleValue = 100 * NumericTraits< double >::epsilon();
    delta = parameters[i] - m_PreviousParameters[i];
    if( std::abs( delta ) > acceptibleValue )
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

  for( unsigned int i = 0; i < ImageDimension; ++i )
    {
    m_PreviousParametersPtr[i] = parameters[i];
    }
}

} // end namespace BlockMatching
} // end namespace itk

#endif
