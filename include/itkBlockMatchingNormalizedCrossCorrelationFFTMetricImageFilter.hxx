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

template <class TFixedImage, class TMovingImage, class TMetricImage >
NormalizedCrossCorrelationFFTMetricImageFilter< TFixedImage,
  TMovingImage, TMetricImage >
::NormalizedCrossCorrelationFFTMetricImageFilter():
  m_GreatestPrimeFactor( 13 )
{
  m_PadFilter = PadFilterType::New();
  // Zero pad.
  m_PadFilter->SetPadMethod( 2 );

  m_FFTShiftFilter  = FFTShiftFilterType::New();
  m_FFTShiftFilter->SetInput( m_PadFilter->GetOutputKernel() );
  m_FFTShiftFilter->SetInverse( true );

  m_KernelFFTFilter = FFTFilterType::New();
  m_KernelFFTFilter->SetInput( m_FFTShiftFilter->GetOutput() );
  m_MovingFFTFilter = FFTFilterType::New();
  m_MovingFFTFilter->SetInput( m_PadFilter->GetOutput() );

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


template <class TFixedImage, class TMovingImage, class TMetricImage >
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
  m_PadFilter->SetInput( movingMinusMean );
  m_PadFilter->SetInputKernel( fixedMinusMean );
  m_PadFilter->SetNumberOfThreads( this->GetNumberOfThreads() );
  // vnl filters need a size which is a power of 2
  if( std::string( m_KernelFFTFilter->GetNameOfClass()).find("Vnl") == 0 )
    {
    m_PadFilter->SetPadToPowerOfTwo( true );
    }
  else
    {
    m_PadFilter->SetGreatestPrimeFactor( m_GreatestPrimeFactor );
    }

  m_FFTShiftFilter->SetNumberOfThreads( this->GetNumberOfThreads() );

  m_KernelFFTFilter->SetNumberOfThreads( this->GetNumberOfThreads() );
  m_MovingFFTFilter->SetNumberOfThreads( this->GetNumberOfThreads() );
  m_ComplexConjugateImageFilter->SetNumberOfThreads( this->GetNumberOfThreads() );

  m_MultiplyFilter->SetNumberOfThreads( this->GetNumberOfThreads() );

  m_IFFTFilter->SetNumberOfThreads( this->GetNumberOfThreads() );

  m_CropFilter->SetReferenceImage( denom );
  m_CropFilter->SetNumberOfThreads( this->GetNumberOfThreads() );

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
