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
#ifndef itkBlockMatchingNormalizedCrossCorrelationFFTMetricImageFilter_hxx
#define itkBlockMatchingNormalizedCrossCorrelationFFTMetricImageFilter_hxx

#include "itkBlockMatchingNormalizedCrossCorrelationFFTMetricImageFilter.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

namespace itk
{
namespace BlockMatching
{

template< typename TFixedImage, typename TMovingImage, typename TMetricImage >
NormalizedCrossCorrelationFFTMetricImageFilter< TFixedImage,
  TMovingImage, TMetricImage >
::NormalizedCrossCorrelationFFTMetricImageFilter()
{
  // Zero pad.
  m_KernelPadFilter = PadFilterType::New();
  m_MovingPadFilter = PadFilterType::New();

  m_FFTShiftFilter  = FFTShiftFilterType::New();
  m_FFTShiftFilter->SetInput( m_KernelPadFilter->GetOutput() );
  m_FFTShiftFilter->SetInverse( true );

  m_KernelFFTFilter = FFTFilterType::New();
  m_KernelFFTFilter->SetInput( m_FFTShiftFilter->GetOutput() );
  m_MovingFFTFilter = FFTFilterType::New();
  m_MovingFFTFilter->SetInput( m_MovingPadFilter->GetOutput() );

  m_SizeGreatestPrimeFactor = m_MovingFFTFilter->GetSizeGreatestPrimeFactor();

  m_ComplexConjugateImageFilter = ComplexConjugateFilterType::New();
  m_ComplexConjugateImageFilter->SetInput( m_KernelFFTFilter->GetOutput() );

  m_MultiplyFilter  = MultiplyFilterType::New();
  m_MultiplyFilter->SetInput1( m_ComplexConjugateImageFilter->GetOutput() );
  m_MultiplyFilter->SetInput2( m_MovingFFTFilter->GetOutput() );
  m_MultiplyFilter->SetInPlace( true );

  m_IFFTFilter = IFFTFilterType::New();
  m_IFFTFilter->SetInput( m_MultiplyFilter->GetOutput() );

  m_CropFilter = CropFilterType::New();
  m_CropFilter->SetInput( m_IFFTFilter->GetOutput() );
}


template< typename TFixedImage, typename TMovingImage, typename TMetricImage >
void
NormalizedCrossCorrelationFFTMetricImageFilter< TFixedImage,
  TMovingImage, TMetricImage >
::GenerateData()
{
  this->AllocateOutputs();
  this->GenerateHelperImages();

  FixedImageConstPointerType  fixedPtr  = this->GetInput( 0 );
  MovingImageConstPointerType movingPtr = this->GetInput( 1 );
  if( !fixedPtr || !movingPtr )
    {
    return;
    }

  // Our output and helper images that have been pre-computed.
  MetricImagePointerType      metricPtr       = this->GetOutput( 0 );
  MetricImageConstPointerType denom           = this->GetOutput( 1 );
  MetricImagePointerType      fixedMinusMean  = this->GetOutput( 2 );
  MetricImageConstPointerType movingMinusMean = this->GetOutput( 3 );

  // The moving search region for this thread.
  m_MovingPadFilter->SetInput( movingMinusMean );
  m_KernelPadFilter->SetInput( fixedMinusMean );

  const MetricImageRegionType & fixedMinusMeanRegion = fixedMinusMean->GetLargestPossibleRegion();
  const typename MetricImageRegionType::SizeType & fixedMinusMeanSize = fixedMinusMeanRegion.GetSize();
  const typename MetricImageRegionType::IndexType & fixedMinusMeanIndex = fixedMinusMeanRegion.GetIndex();
  const MetricImageRegionType & movingMinusMeanRegion = movingMinusMean->GetLargestPossibleRegion();
  const typename MetricImageRegionType::SizeType & movingMinusMeanSize = movingMinusMeanRegion.GetSize();
  const typename MetricImageRegionType::IndexType & movingMinusMeanIndex = movingMinusMeanRegion.GetIndex();
  typename MetricImageRegionType::IndexType paddedIndex;
  typename MetricImageRegionType::SizeType paddedSize;
  for( unsigned int ii = 0; ii < ImageDimension; ++ii )
    {
    SizeValueType padSize = std::max( static_cast< SizeValueType >( 0 ), fixedMinusMeanSize[ii] - 1 );
    if( m_SizeGreatestPrimeFactor > 1 )
      {
      while( Math::GreatestPrimeFactor( movingMinusMeanSize[ii] + padSize ) > m_SizeGreatestPrimeFactor )
        {
        ++padSize;
        }
      }
    else if( m_SizeGreatestPrimeFactor == 1 )
      {
      // make sure the total size is even
      padSize += ( movingMinusMeanSize[ii] + padSize ) % 2;
      }
    paddedIndex[ii] = movingMinusMeanIndex.GetIndex()[ii] - padSize/2;
    paddedSize[ii] = movingMinusMeanSize[ii] + padSize;
    }

  typename MetricImageRegionType::SizeType padding;
  for( unsigned int ii = 0; ii < ImageDimension; ++ii )
    {
    padding[ii] = movingMinusMeanIndex[ii] - paddedIndex[ii];
    }
  m_MovingPadFilter->SetPadLowerBound( padding );
  for( unsigned int ii = 0; ii < ImageDimension; ++ii )
    {
    padding[ii] = paddedSize[ii] - ( movingMinusMeanIndex[ii] - paddedIndex[ii] + movingMinusMeanSize[ii] );
    }
  m_MovingPadFilter->SetPadUpperBound( padding );
  for( unsigned int ii = 0; ii < ImageDimension; ++ii )
    {
    padding[ii] = fixedMinusMeanIndex[ii] - paddedIndex[ii];
    }
  m_KernelPadFilter->SetPadLowerBound( padding );
  for( unsigned int ii = 0; ii < ImageDimension; ++ii )
    {
    padding[ii] = paddedSize[ii] - ( fixedMinusMeanIndex[ii] - paddedIndex[ii] + fixedMinusMeanSize[ii] );
    }
  m_KernelPadFilter->SetPadUpperBound( padding );

  m_MovingPadFilter->SetNumberOfWorkUnits( this->GetNumberOfWorkUnits() );
  m_KernelPadFilter->SetNumberOfWorkUnits( this->GetNumberOfWorkUnits() );
  m_FFTShiftFilter->SetNumberOfWorkUnits( this->GetNumberOfWorkUnits() );
  m_KernelFFTFilter->SetNumberOfWorkUnits( this->GetNumberOfWorkUnits() );
  m_MovingFFTFilter->SetNumberOfWorkUnits( this->GetNumberOfWorkUnits() );
  m_ComplexConjugateImageFilter->SetNumberOfWorkUnits( this->GetNumberOfWorkUnits() );
  m_MultiplyFilter->SetNumberOfWorkUnits( this->GetNumberOfWorkUnits() );
  m_IFFTFilter->SetNumberOfWorkUnits( this->GetNumberOfWorkUnits() );

  m_CropFilter->SetReferenceImage( denom );
  m_CropFilter->SetNumberOfWorkUnits( this->GetNumberOfWorkUnits() );

  m_CropFilter->UpdateLargestPossibleRegion();

  typedef ImageRegionConstIterator< MetricImageType > ConstIteratorType;
  typedef ImageRegionIterator< MetricImageType >      IteratorType;

  ConstIteratorType denomIt( denom, this->m_MovingImageRegion );

  ConstIteratorType corrIt( m_CropFilter->GetOutput(), this->m_MovingImageRegion );
  IteratorType      metricIt( metricPtr, metricPtr->GetLargestPossibleRegion() );

  const MetricImagePixelType negativeOne = -1 * NumericTraits< MetricImagePixelType >::One;
  const MetricImagePixelType positiveOne = NumericTraits< MetricImagePixelType >::One;
  MetricImagePixelType       normXcorr;
  for( denomIt.GoToBegin(), corrIt.GoToBegin(), metricIt.GoToBegin();
       !metricIt.IsAtEnd();
       ++denomIt, ++corrIt, ++metricIt )
    {
    if( denomIt.Get() == NumericTraits< MetricImagePixelType >::Zero )
      {
      metricIt.Set( NumericTraits< MetricImagePixelType >::Zero );
      }
    else
      {
      normXcorr = corrIt.Get() / denomIt.Get();
      // Why does this happen?  Bug?  Funky floating point behavior?
      if( normXcorr < negativeOne || normXcorr > positiveOne )
        {
        metricIt.Set( NumericTraits< MetricImagePixelType >::Zero );
        }
      else
        {
        metricIt.Set( normXcorr );
        }
      }
    }
}

} // end namespace BlockMatching
} //end namespace itk

#endif
