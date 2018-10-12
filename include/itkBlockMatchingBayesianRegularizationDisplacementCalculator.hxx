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
#ifndef itkBlockMatchingBayesianRegularizationDisplacementCalculator_hxx
#define itkBlockMatchingBayesianRegularizationDisplacementCalculator_hxx

#include "itkBlockMatchingBayesianRegularizationDisplacementCalculator.h"
#include "itkBlockMatchingMaximumPixelDisplacementCalculator.h"

#include "itkNeighborhoodAlgorithm.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkCompensatedSummation.h"

namespace
{

template <typename T>
struct identity {
  using type = T;
};

template <typename T>
using IdentityType = typename identity<T>::type;

template< typename TFloat >
TFloat
AtomicCompensatedSummation( std::atomic< TFloat > & sum, IdentityType< TFloat > addition )
{
  if( addition == static_cast< TFloat >( 0.0 ) )
    {
    return sum;
    }
  TFloat old = sum;
  using CompensatedSummationType = itk::CompensatedSummation< TFloat >;
  CompensatedSummationType desired;
  do {
    desired.ResetToZero();
    desired += old;
    desired += addition;
  } while ( !sum.compare_exchange_weak( old, desired.GetSum() ) );
  return desired.GetSum();
}

}

namespace itk
{
namespace BlockMatching
{

template< typename TMetricImage, typename TDisplacementImage >
BayesianRegularizationDisplacementCalculator< TMetricImage, TDisplacementImage >
::BayesianRegularizationDisplacementCalculator():
  m_MaximumIterations( 3 ),
  m_MeanChangeThreshold( 2.0 ),
  m_MeanChangeThresholdDefined( false ),
  m_MeanChange( 3.0 ),
  m_MetricLowerBoundDefined( false ),
  m_StrainSigma( 0.1 ),
  m_MaximumStrain( 0.0 ),
  m_CurrentIteration( 0 )
{
  this->m_CacheMetricImage = true;

  m_DisplacementCalculator =
    MaximumPixelDisplacementCalculator< TMetricImage,
                                        TDisplacementImage >::New();

  m_PriorPr = nullptr;

  for( unsigned int dim = 0; dim < ImageDimension; dim++ )
    m_GaussianKernels.push_back( MetricImageType::New() );

  m_GaussianKernelRadii.resize( ImageDimension );
}


template< typename TMetricImage, typename TDisplacementImage >
void
BayesianRegularizationDisplacementCalculator< TMetricImage, TDisplacementImage >
::AllocatePriorPrImage()
{
  typedef typename itk::ImageRegionIterator< MetricImageImageType >
    PrIteratorType;
  typedef typename itk::ImageRegionConstIterator< MetricImageImageType >
    PrConstIteratorType;
  MetricImagePointerType imagePr;
  m_PriorPr = MetricImageImageType::New();
  m_PriorPr->CopyInformation( this->m_MetricImageImage );
  m_PriorPr->SetRegions(
    this->m_MetricImageImage->GetLargestPossibleRegion() );
  m_PriorPr->Allocate();
  PrIteratorType priorIt( m_PriorPr, m_PriorPr->GetLargestPossibleRegion() );
  PrConstIteratorType postIt( this->m_MetricImageImage,
                              this->m_MetricImageImage->GetLargestPossibleRegion() );

  for( priorIt.GoToBegin(), postIt.GoToBegin();
       !priorIt.IsAtEnd();
       ++priorIt, ++postIt )
    {
    imagePr = MetricImageType::New();
    imagePr->CopyInformation( postIt.Get() );
    imagePr->SetRegions( postIt.Get()->GetLargestPossibleRegion() );
    imagePr->Allocate();
    priorIt.Set( imagePr );
    }
}


template< typename TMetricImage, typename TDisplacementImage >
void
BayesianRegularizationDisplacementCalculator< TMetricImage, TDisplacementImage >
::ScaleToUnity()
{
  this->m_MultiThreader->template ParallelizeImageRegion<ImageDimension>(
    this->m_DisplacementImage->GetRequestedRegion(),
    [this]( const RegionType & outputRegion )
        {
        this->ThreadedScaleToUnity( outputRegion );
        },
    nullptr);

  if( m_MeanChangeThresholdDefined && m_CurrentIteration > 0 )
    {
    m_ChangeSum.store( 0.0 );
    m_ChangeCount.store( 0 );
    // This adds up m_ChangeSum and m_ChangeCount.
    this->m_MultiThreader->template ParallelizeImageRegion<ImageDimension>(
      this->m_DisplacementImage->GetRequestedRegion(),
      [this]( const RegionType & outputRegion )
          {
          this->ThreadedMeanChange( outputRegion );
          },
      nullptr);
    m_MeanChange = m_ChangeSum.load() / static_cast< double >( m_ChangeCount.load() );
    }
}


template< typename TMetricImage, typename TDisplacementImage >
void
BayesianRegularizationDisplacementCalculator< TMetricImage, TDisplacementImage >
::GenerateGaussianLikeKernels()
{
  typename GaussianKernelArrayType::iterator it;
  typename GaussianKernelType::RegionType    region;
  typename GaussianKernelType::IndexType     index;
  typename GaussianKernelType::SizeType      size;
  typename GaussianKernelType::Pointer       image;
  GaussianKernelRadiusType                   radius;
  // This is where the equal metric image spacing assumption comes in.  We could
  // avoid the assumption, but then we may have to regenerate the gaussian
  // kernels every time.
  typename MetricImageType::IndexType dummyIndex;
  dummyIndex.Fill( 0 );
  SpacingType maxStrain;
  SpacingType metricSpacing =
    this->m_MetricImageImage->GetPixel( dummyIndex )->GetSpacing();
  SpacingType displacementSpacing =
    this->m_DisplacementImage->GetSpacing();
  unsigned int dimIt;
  unsigned int dimIt2;
  PixelType    temp;
  PixelType    temp2;
  PixelType    neighborValue;

  // Set the maximum strain if it has not been specified.
  for( dimIt = 0; dimIt < ImageDimension; dimIt++ )
    {
    if( m_MaximumStrain[dimIt] == 0.0 )
      maxStrain[dimIt] = 3.0 * m_StrainSigma[dimIt];
    else
      maxStrain[dimIt] = m_MaximumStrain[dimIt];
    }

  for( it = m_GaussianKernels.begin(), dimIt = 0;
       it != m_GaussianKernels.end();
       ++it, ++dimIt )
    {
    image = *it;
    for( dimIt2 = 0; dimIt2 < ImageDimension; dimIt2++ )
      {
      radius[dimIt2] = static_cast< typename GaussianKernelRadiusType::SizeValueType >( std::ceil(
      displacementSpacing[dimIt] * maxStrain[dimIt2] / metricSpacing[dimIt2] ));
      size[dimIt2] = radius[dimIt2] * 2 + 1;
      index[dimIt2] = -1 * radius[dimIt2];
      }
    m_GaussianKernelRadii[dimIt] = radius;
    region.SetSize( size );
    region.SetIndex( index );
    image->SetRegions( region );
    image->Allocate();

    itk::ImageRegionIteratorWithIndex< GaussianKernelType >
      imageIt( image, region );
    for( imageIt.GoToBegin(); !imageIt.IsAtEnd(); ++imageIt )
      {
      neighborValue = 0.0;
      index = imageIt.GetIndex();
      for( dimIt2 = 0; dimIt2 < ImageDimension; dimIt2++ )
        {
        temp = std::abs( index[dimIt2] ) * metricSpacing[dimIt2];
        temp2 = displacementSpacing[dimIt] * m_StrainSigma[dimIt2];
        neighborValue += temp * temp / ( temp2 * temp2 );
        }
      neighborValue *= -0.5;
      imageIt.Set( std::exp( neighborValue ) );
      }
    }
}


template< typename TMetricImage, typename TDisplacementImage >
void
BayesianRegularizationDisplacementCalculator< TMetricImage, TDisplacementImage >
::ImpartLikelihood( MetricImagePointerType& postImage,
                    MetricImagePointerType& priorImage,
                    const unsigned int direction,
                    const VectorType& shift )
{
  // We apply Neumann padding to the prior.
  SizeType lowerBoundPad;
  lowerBoundPad.Fill( 0 );
  SizeType upperBoundPad;
  upperBoundPad.Fill( 0 );
  IndexType priorLeftIndex = priorImage->GetBufferedRegion().GetIndex();
  IndexType priorRightIndex = priorLeftIndex +
    priorImage->GetBufferedRegion().GetSize();
  RegionType postRegion = postImage->GetBufferedRegion();
  postRegion.PadByRadius( m_GaussianKernelRadii[direction] );
  IndexType    postLeftIndex = postRegion.GetIndex();
  IndexType    postRightIndex = postLeftIndex + postRegion.GetSize();
  unsigned int dim;
  unsigned int padding;
  PointType    point;
  IndexType    index;
  postImage->TransformIndexToPhysicalPoint( postLeftIndex, point );
  priorImage->TransformPhysicalPointToIndex( point + shift, index );
  for( dim = 0; dim < ImageDimension; dim++ )
    {
    if( index[dim] < priorLeftIndex[dim] )
      {
      padding = priorLeftIndex[dim] - index[dim];
      if( padding > lowerBoundPad[dim] )
        lowerBoundPad[dim] = padding;
      }
    if( index[dim] > priorRightIndex[dim] )
      {
      padding = index[dim] - priorRightIndex[dim];
      if( padding > upperBoundPad[dim] )
      upperBoundPad[dim] = padding;
      }
    }
  postImage->TransformIndexToPhysicalPoint( postRightIndex, point );
  priorImage->TransformPhysicalPointToIndex( point + shift, index );
  for( dim = 0; dim < ImageDimension; dim++ )
    {
    if( index[dim] < priorLeftIndex[dim] )
      {
      padding = priorLeftIndex[dim] - index[dim];
      if( padding > lowerBoundPad[dim] )
        lowerBoundPad[dim] = padding;
      }
    if( index[dim] > priorRightIndex[dim] )
      {
      padding = index[dim] - priorRightIndex[dim];
      if( padding > upperBoundPad[dim] )
      upperBoundPad[dim] = padding;
      }
    }
  typename PadFilterType::Pointer padFilter;
  padFilter = PadFilterType::New();
  padFilter->SetInput( priorImage );
  padFilter->SetPadLowerBound( lowerBoundPad );
  padFilter->SetPadUpperBound( upperBoundPad );
  padFilter->UpdateLargestPossibleRegion();
  MetricImagePointerType paddedPrior = padFilter->GetOutput();

  MetricImageIteratorType kernelIt( m_GaussianKernels[direction],
                                    m_GaussianKernels[direction]->GetLargestPossibleRegion() );

  postRegion = postImage->GetBufferedRegion();
  MetricImageIteratorType postIt( postImage, postRegion );

  postImage->TransformIndexToPhysicalPoint( postRegion.GetIndex(), point );
  paddedPrior->TransformPhysicalPointToIndex( point + shift, index );
  RegionType priorRegion;
  priorRegion.SetIndex( index );
  priorRegion.SetSize( postRegion.GetSize() );
  typedef NeighborhoodIterator< MetricImageType > NeighborhoodIteratorType;
  NeighborhoodIteratorType priorIt( m_GaussianKernelRadii[direction],
                                    paddedPrior, priorRegion );

  const unsigned int neighborhoodSize = priorIt.Size();

  PixelType       maxPr;
  PixelType       localPr;
  const PixelType nonPositiveMin = NumericTraits< PixelType >::NonpositiveMin();
  unsigned int    neighbor;
  for( priorIt.GoToBegin(), postIt.GoToBegin();
       !postIt.IsAtEnd();
       ++priorIt, ++postIt )
    {
    maxPr = nonPositiveMin;
    for( neighbor = 0, kernelIt.GoToBegin();
         neighbor < neighborhoodSize;
         ++neighbor, ++kernelIt )
      {
      localPr = priorIt.GetPixel( neighbor ) * kernelIt.Get();
      if( localPr > maxPr )
        maxPr = localPr;
      }
    postIt.Value() *= maxPr;
    }
}


template< typename TMetricImage, typename TDisplacementImage >
void
BayesianRegularizationDisplacementCalculator< TMetricImage, TDisplacementImage >
::Compute()
{
  // First we shift the minimum value of the metric image so 0 corresponds to
  // the theoretical lower bound.
  if( !m_MetricLowerBoundDefined )
    {
    itkExceptionMacro( << "The metric lower bound must be specified." );
    }
  this->m_MultiThreader->template ParallelizeImageRegion<ImageDimension>(
    this->m_DisplacementImage->GetRequestedRegion(),
    [this]( const RegionType & outputRegion )
        {
        this->ThreadedSubtractLowerBound( outputRegion );
        },
    nullptr);

  // Scale so every metric image sums to unity like it is a probability image.
  this->ScaleToUnity();

  this->AllocatePriorPrImage();

  this->GenerateGaussianLikeKernels();

  // The radius for regularization.
  typename  MetricImageImageType::SizeType regRadius;
  regRadius.Fill( 1 );

  typedef typename NeighborhoodAlgorithm::
    ImageBoundaryFacesCalculator< MetricImageImageType > FaceCalculatorType;
  FaceCalculatorType faceCalculator;
  typename FaceCalculatorType::FaceListType faceList = faceCalculator(
    this->m_MetricImageImage,
    this->m_MetricImageImage->GetLargestPossibleRegion(), regRadius );
  typename FaceCalculatorType::FaceListType::iterator fit;

  MetricImageImagePointerType tempMetricImageImagePtr;
  MetricImagePointerType      postImage;
  MetricImagePointerType      priorImage;
  m_CurrentIteration = 0;
  // Hack to make the SplitRequestedRegion method operate on the face list
  // regions;
  const RegionType requestedRegion =
    this->m_DisplacementImage->GetRequestedRegion();
  while( m_CurrentIteration < this->m_MaximumIterations &&
                              m_MeanChange > m_MeanChangeThreshold )
    {
    // We evoke iteration events starting from 0,
    // when no regularization has occured yet.
    this->InvokeEvent( IterationEvent() );

    // switcheroo
    tempMetricImageImagePtr = this->m_PriorPr;
    m_PriorPr = this->m_MetricImageImage;
    this->m_MetricImageImage = tempMetricImageImagePtr;

    this->m_MultiThreader->template ParallelizeImageRegion<ImageDimension>(
      this->m_DisplacementImage->GetRequestedRegion(),
      [this]( const RegionType & outputRegion )
          {
          this->ThreadedCopyPriorToPosterior( outputRegion );
          },
      nullptr);

    const ThreadIdType workUnits = this->m_MultiThreader->GetNumberOfWorkUnits();
    for( fit = faceList.begin(); fit != faceList.end(); ++fit )
      {
      // TODO: FIX-ME: hangs if this is not set
      this->m_MultiThreader->SetNumberOfWorkUnits( 2 );
      this->m_MultiThreader->template ParallelizeImageRegion<ImageDimension>(
        *fit,
        [this]( const RegionType & outputRegion )
            {
            this->ThreadedImpartLikelihood( outputRegion );
            },
        nullptr);
      }
    this->m_MultiThreader->SetNumberOfWorkUnits( workUnits );

    ++m_CurrentIteration;

    this->ScaleToUnity();
    }

  // Calculate the displacements from the regularized probablity images.
  m_DisplacementCalculator->SetDisplacementImage( this->m_DisplacementImage );
  itk::ImageRegionConstIteratorWithIndex< MetricImageImageType >
  metricImageImageConstIt( this->m_MetricImageImage,
                           this->m_MetricImageImage->GetLargestPossibleRegion() );

  itk::ImageRegionConstIterator< CenterPointsImageType >
    centerPointsConstIt( this->m_CenterPointsImage,
      this->m_CenterPointsImage->GetLargestPossibleRegion() );
  for( metricImageImageConstIt.GoToBegin(), centerPointsConstIt.GoToBegin();
       !metricImageImageConstIt.IsAtEnd();
       ++metricImageImageConstIt, ++centerPointsConstIt )
    {
    this->m_DisplacementCalculator->SetMetricImagePixel(
      centerPointsConstIt.Get(), metricImageImageConstIt.GetIndex(),
      metricImageImageConstIt.Get() );
    }
  this->m_DisplacementCalculator->Compute();
}


template< typename TMetricImage, typename TDisplacementImage >
void
BayesianRegularizationDisplacementCalculator< TMetricImage, TDisplacementImage >
::ThreadedImpartLikelihood( const RegionType &splitRegion )
{
  MetricImageImageIteratorType imageImageIt( this->m_MetricImageImage,
    splitRegion );
  // The radius for regularization.
  typename  MetricImageImageType::SizeType regRadius;
  regRadius.Fill( 1 );
  MetricImageImageNeighborhoodIteratorType priorImageImageIt( regRadius, this->m_PriorPr, splitRegion );
  priorImageImageIt.OverrideBoundaryCondition( &(this->m_ImageImageBoundaryCondition ) );

  MetricImagePointerType postImage;
  MetricImagePointerType priorImage;

  unsigned int direction;
  VectorType   shift;
  const SpacingType  spacing = this->m_DisplacementImage->GetSpacing();

  for( imageImageIt.GoToBegin(), priorImageImageIt.GoToBegin();
       !imageImageIt.IsAtEnd();
       ++imageImageIt, ++priorImageImageIt )
    {
    postImage = imageImageIt.Get();

    // perform regularization along every direction;
    for( direction = 0; direction < ImageDimension; ++direction )
      {
      priorImage = priorImageImageIt.GetPrevious( direction );
      // If we are inside the boundary.
      if( priorImage.GetPointer() != nullptr )
        {
        shift.Fill( 0.0 );
        shift[direction] = -1 * spacing[direction];
        this->ImpartLikelihood( postImage, priorImage, direction, shift );
        }
      priorImage = priorImageImageIt.GetNext( direction );
      if( priorImage.GetPointer() != nullptr )
        {
        shift.Fill( 0.0 );
        shift[direction] = spacing[direction];
        this->ImpartLikelihood( postImage, priorImage, direction, shift );
        }
      } // for every direction
    }   // for every probability image
}


template< typename TMetricImage, typename TDisplacementImage >
void
BayesianRegularizationDisplacementCalculator< TMetricImage, TDisplacementImage >
::ThreadedSubtractLowerBound( const RegionType& region )
{
  MetricImageImageIteratorType imageImageIt( this->m_MetricImageImage,
    region );
  MetricImagePointerType image;
  // We add eps so we don't take the log of zero, and having a probability of
  // 0.0 is ... pessimistic ;-)
  const PixelType lowerBound = this->m_MetricLowerBound +
    NumericTraits< PixelType >::epsilon();
  for( imageImageIt.GoToBegin(); !imageImageIt.IsAtEnd(); ++imageImageIt )
    {
    image = imageImageIt.Get();
    MetricImageIteratorType it( image,
                                image->GetLargestPossibleRegion() );

    for( it.GoToBegin(); !it.IsAtEnd(); ++it )
      {
      it.Value() -= lowerBound;
      }
    }
}


template< typename TMetricImage, typename TDisplacementImage >
void
BayesianRegularizationDisplacementCalculator< TMetricImage, TDisplacementImage >
::ThreadedScaleToUnity( const RegionType& region )
{
  MetricImageImageIteratorType imageImageIt( this->m_MetricImageImage,
    region );
  MetricImagePointerType image;
  for( imageImageIt.GoToBegin(); !imageImageIt.IsAtEnd(); ++imageImageIt )
    {
    image = imageImageIt.Get();
    PixelType sum = NumericTraits< PixelType >::Zero;
    MetricImageIteratorType it( image,
                                image->GetLargestPossibleRegion() );

    for( it.GoToBegin(); !it.IsAtEnd(); ++it )
      {
      sum += it.Get();
      }

    for( it.GoToBegin(); !it.IsAtEnd(); ++it )
      {
      it.Value() /= sum;
      }
    }
}


template< typename TMetricImage, typename TDisplacementImage >
void
BayesianRegularizationDisplacementCalculator< TMetricImage, TDisplacementImage >
::ThreadedMeanChange( const RegionType& region )
{
  MetricImageImageIteratorType imageImageIt( this->m_MetricImageImage,
                                             region );

  MetricImageImageIteratorType priorImageImageIt( this->m_PriorPr,
    region );
  MetricImagePointerType image;
  MetricImagePointerType priorPtr;
  double                 changeSum = 0.0;
  unsigned long long     changeCount = 0;
  for( imageImageIt.GoToBegin(), priorImageImageIt.GoToBegin();
       !imageImageIt.IsAtEnd();
       ++imageImageIt, ++priorImageImageIt )
    {
    image = imageImageIt.Get();
    priorPtr = priorImageImageIt.Get();

    MetricImageConstIteratorType it( image,
                                     image->GetLargestPossibleRegion() );

    MetricImageConstIteratorType priorIt( priorPtr,
                                          priorPtr->GetLargestPossibleRegion() );
    for( it.GoToBegin(), priorIt.GoToBegin();
         !it.IsAtEnd();
         ++it, ++priorIt )
      {
      changeSum += std::abs( it.Get() - std::exp( priorIt.Get() ) );
      ++changeCount;
      }
    }
  AtomicCompensatedSummation< double >( this->m_ChangeSum, changeSum );
  this->m_ChangeCount.fetch_add( changeCount );
}

template< typename TMetricImage, typename TDisplacementImage >
void
BayesianRegularizationDisplacementCalculator< TMetricImage, TDisplacementImage >
::ThreadedCopyPriorToPosterior( const RegionType& region )
{
  MetricImageImageIteratorType imageImageIt( this->m_MetricImageImage,
                                             region );

  MetricImageImageIteratorType priorImageImageIt( this->m_PriorPr,
    region );
  MetricImagePointerType image;
  MetricImagePointerType priorPtr;
  for( imageImageIt.GoToBegin(), priorImageImageIt.GoToBegin();
       !imageImageIt.IsAtEnd();
       ++imageImageIt, ++priorImageImageIt )
    {
    image = imageImageIt.Get();
    priorPtr = priorImageImageIt.Get();

    MetricImageIteratorType it( image,
                                     image->GetLargestPossibleRegion() );

    MetricImageConstIteratorType priorIt( priorPtr,
                                          priorPtr->GetLargestPossibleRegion() );
    for( it.GoToBegin(), priorIt.GoToBegin();
         !it.IsAtEnd();
         ++it, ++priorIt )
      {
      it.Set( priorIt.Get() );
      }
    }
}

} // end namespace BlockMatching
} // end namespace itk

#endif
