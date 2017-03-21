#ifndef itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter_hxx
#define itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter_hxx

#include "itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter.h"

#include "itkConstantBoundaryCondition.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImageKernelOperator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkNeighborhoodInnerProduct.h"

namespace itk
{
namespace BlockMatching
{

template <class TFixedImage, class TMovingImage,
          class TMetricImage >
void
NormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
::BeforeThreadedGenerateData()
{
  this->GenerateHelperImages();
}

template <class TFixedImage, class TMovingImage,
          class TMetricImage >
void
NormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
::ThreadedGenerateData( const MetricImageRegionType& outputRegion, int threadId )
{
  FixedImageConstPointerType  fixedPtr   = this->GetInput( 0 );
  MovingImageConstPointerType movingPtr = this->GetInput( 1 );

  if( !fixedPtr || !movingPtr )
    return;

  // Our output and helper images that have been pre-computed.
  MetricImagePointerType      metricPtr       = this->GetOutput( 0 );
  MetricImageConstPointerType denom           = this->GetOutput( 1 );
  MetricImagePointerType      fixedMinusMean  = this->GetOutput( 2 );
  MetricImageConstPointerType movingMinusMean = this->GetOutput( 3 );

  // This is in case a truncated radius does not cover everything or the denom
  // is zero.
  metricPtr->FillBuffer( NumericTraits< MetricImagePixelType >::Zero );
  typename MovingImageType::SizeType radius = this->m_MovingRadius;
  typename MovingImageType::SizeType bufferedMeanSize = movingMinusMean->GetBufferedRegion().GetSize();
  typename MovingImageType::SizeType movingImageRegionSize = this->m_MovingImageRegion.GetSize();
  unsigned int i;
  for( i = 0; i < ImageDimension; ++i )
    {
    if( 2 * radius[i] + 1 > bufferedMeanSize[i] )
      {
      radius[i] = Math::Floor< typename MovingImageType::SizeType::SizeValueType >( ( bufferedMeanSize[i] - 1.0 ) / 2.0 );
      }
    if( 2 * radius[i] + 1 > movingImageRegionSize[i] )
      {
      radius[i] = Math::Floor< typename MovingImageType::SizeType::SizeValueType >( ( movingImageRegionSize[i] - 1.0 ) / 2.0 );
      }
    }
  //// The outputRegion's index corresponds to the metric image index -- which
  //// starts from 0 at the metric and (usually) somewhere else in the moving
  //// image.
  typename MetricImageRegionType::IndexType metricIndex = outputRegion.GetIndex();
  typename MovingImageType::IndexType movingRegionIndex = this->m_MovingImageRegion.GetIndex();

  // The types of our iterators.
  typedef ImageRegionIterator< MetricImageType > MetricIteratorType;
  MetricImageRegionType metricRegion;
  typename MovingImageType::IndexType fitIndex;
  typedef ImageRegionConstIterator< MetricImageType > MetricConstIteratorType;
  typedef ConstantBoundaryCondition< MetricImageType > BoundaryConditionType;
  BoundaryConditionType boundaryCondition;
  boundaryCondition.SetConstant( NumericTraits< MetricImagePixelType >::Zero );
  typedef ConstNeighborhoodIterator< MetricImageType > NeighborhoodIteratorType;

  ImageKernelOperator< typename MetricImageType::PixelType, ImageDimension > fixedKernelOperator;
  // in case the radius was truncated we need to truncate the fixedMinusMean
  // size because the ImageKernelOperator is currently broken if CreateToRadius
  // does not have a radius corresponding to the image size.
  typename MovingImageType::SizeType fixedMinusMeanSize;
  for( i = 0; i < ImageDimension; ++i )
    {
    fixedMinusMeanSize[i] = 2 * radius[i] + 1;
    }
  typename MovingImageType::RegionType fixedMinusMeanRegion;
  fixedMinusMeanRegion.SetSize( fixedMinusMeanSize );
  fixedMinusMean->SetRegions( fixedMinusMeanRegion );
  fixedKernelOperator.SetImageKernel( fixedMinusMean );
  fixedKernelOperator.CreateToRadius( radius );
  NeighborhoodInnerProduct< MetricImageType > innerProduct;

  typedef typename NeighborhoodAlgorithm::
  ImageBoundaryFacesCalculator< MovingImageType > FaceCalculatorType;
  FaceCalculatorType                        faceCalculator;
  typename FaceCalculatorType::FaceListType faceList = faceCalculator(
      movingPtr, this->m_MovingImageRegion, radius );
  typename FaceCalculatorType::FaceListType::iterator fit;
  const MetricImagePixelType                          negativeOne = -1 * NumericTraits< MetricImagePixelType >::One;
  const MetricImagePixelType                          positiveOne = NumericTraits< MetricImagePixelType >::One;
  MetricImagePixelType                                normXcorr;
  for( fit = faceList.begin(); fit != faceList.end(); ++fit )
    {
    (*fit).Crop( outputRegion );

    // Unlike the other regions, the Index on the metric image was set to 0 at
    // the beginning of the metric image.
    metricRegion = *fit;
    fitIndex = metricRegion.GetIndex();
    for( i = 0; i < ImageDimension; i++ )
      {
      metricIndex[i] = fitIndex[i] - movingRegionIndex[i];
      }
    metricRegion.SetIndex( metricIndex );

    MetricIteratorType metricIt( metricPtr, metricRegion );

    MetricConstIteratorType denomIt( denom, *fit );
    NeighborhoodIteratorType movingNeighborIt( radius, movingMinusMean, *fit );
    movingNeighborIt.OverrideBoundaryCondition( &boundaryCondition );

    // for every pixel in the output
    for( metricIt.GoToBegin(), denomIt.GoToBegin(), movingNeighborIt.GoToBegin();
         !metricIt.IsAtEnd();
         ++metricIt, ++denomIt, ++movingNeighborIt )
      {
      if( !(denomIt.Get() == NumericTraits< MetricImagePixelType >::Zero) )
        {
        normXcorr = innerProduct( movingNeighborIt, fixedKernelOperator ) / denomIt.Get();
        // Why does this happen?  Bug?  Funky floating point behavior?
        if( normXcorr < negativeOne )
          metricIt.Set( negativeOne );
        else if ( normXcorr > positiveOne )
          metricIt.Set( positiveOne );
        else
          metricIt.Set( normXcorr );
        }
      }
    }
}


} // end namespace BlockMatching
} //end namespace itk

#endif
