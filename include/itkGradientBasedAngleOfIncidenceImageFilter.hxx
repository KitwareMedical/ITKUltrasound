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
#ifndef itkGradientBasedAngleOfIncidenceImageFilter_hxx
#define itkGradientBasedAngleOfIncidenceImageFilter_hxx

#include "itkGradientBasedAngleOfIncidenceImageFilter.h"

#include "itkGradientImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkProgressReporter.h"

namespace itk
{

template< typename TInputImage, typename TOutputImage, typename TOperatorValue >
GradientBasedAngleOfIncidenceImageFilter< TInputImage,
  TOutputImage,
  TOperatorValue >
::GradientBasedAngleOfIncidenceImageFilter()
{
  typedef GradientImageFilter<
      OperatorImageType, OperatorValueType, OperatorValueType >
    DefaultGradientFilterType;
  this->m_GradientFilter = DefaultGradientFilterType::New().GetPointer();

  this->m_CastImageFilter = CastImageFilterType::New();

  this->m_UltrasoundProbeType = CURVILINEAR;
  this->m_UltrasoundProbeOrigin.Fill( 0.0 );
  this->m_UltrasoundProbeBeamDirection.Fill( 0.0 );

  this->m_GradientMagnitudeTolerance = 1.0e-7;
  //Use the ITKv4 Threading Model (call ThreadedGenerateData instead of DynamicThreadedGenerateData)
  this->DynamicMultiThreadingOff();
}


template< typename TInputImage, typename TOutputImage, typename TOperatorValue >
void
GradientBasedAngleOfIncidenceImageFilter< TInputImage,
  TOutputImage,
  TOperatorValue >
::SetUltrasoundProbeBeamDirection( const BeamDirectionType & beamDirection )
{
  if( beamDirection != this->m_UltrasoundProbeBeamDirection )
    {
    this->m_UltrasoundProbeBeamDirection = beamDirection;
    this->m_UltrasoundProbeBeamDirection.Normalize();
    this->Modified();
    }
}


template< typename TInputImage, typename TOutputImage, typename TOperatorValue >
void
GradientBasedAngleOfIncidenceImageFilter< TInputImage,
  TOutputImage,
  TOperatorValue >
::BeforeThreadedGenerateData()
{
  this->m_CastImageFilter->SetInput( this->GetInput() );
  this->m_GradientFilter->SetNumberOfWorkUnits( this->GetNumberOfWorkUnits() );
  this->m_GradientFilter->SetInput( this->m_CastImageFilter->GetOutput() );
  this->m_GradientFilter->Update();

  if( this->m_UltrasoundProbeType == LINEAR &&
      this->m_UltrasoundProbeBeamDirection.GetNorm()
        == static_cast< OperatorValueType >( 0.0 ) )
    {
    itkExceptionMacro(
      << "The BeamDirection must be specified with a linear probe." );
    }
}


template< typename TInputImage, typename TOutputImage, typename TOperatorValue >
void
GradientBasedAngleOfIncidenceImageFilter< TInputImage,
  TOutputImage,
  TOperatorValue >
::ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread,
  ThreadIdType threadId )
{
  const InputImageType * input = this->GetInput();
  const GradientOutputImageType * gradient =
    this->m_GradientFilter->GetOutput();
  OutputImageType * output = this->GetOutput();

  const OriginType origin = this->m_UltrasoundProbeOrigin;
  const double gradientMagnitudeTolerance =
    this->m_GradientMagnitudeTolerance;

  typedef ImageRegionConstIteratorWithIndex< InputImageType >
    InputIteratorType;
  typedef ImageRegionConstIterator< GradientOutputImageType >
    GradientIteratorType;
  typedef ImageRegionIterator< OutputImageType >
    OutputIteratorType;

  InputIteratorType inputIt( input, outputRegionForThread );
  GradientIteratorType gradientIt( gradient, outputRegionForThread );
  OutputIteratorType outputIt( output, outputRegionForThread );

  ProgressReporter progress( this, threadId,
    outputRegionForThread.GetNumberOfPixels() );

  BeamDirectionType beamDirection = this->m_UltrasoundProbeBeamDirection;

  for( inputIt.GoToBegin(), gradientIt.GoToBegin(), outputIt.GoToBegin();
       !outputIt.IsAtEnd();
       ++inputIt, ++gradientIt, ++outputIt )
    {
    typename InputImageType::IndexType index = inputIt.GetIndex();
    typename InputImageType::PointType point;
    input->TransformIndexToPhysicalPoint( index, point );

    if( this->m_UltrasoundProbeType == CURVILINEAR ||
        this->m_UltrasoundProbeType == PHASED )
      {
      beamDirection = point - origin;
      beamDirection.Normalize();
      }

    GradientOutputPixelType gradientPixel = gradientIt.Get();
    const typename GradientOutputPixelType::RealValueType gradientNorm =
      gradientPixel.GetNorm();
    gradientPixel /= gradientNorm;

    // output  scalar product of the two normalized vectors
    typedef typename OutputImageType::PixelType OutputPixelType;
    const OutputPixelType outputPixel = gradientPixel * beamDirection;
    if( vnl_math_isnan( outputPixel )
      || gradientNorm < gradientMagnitudeTolerance )
      {
      outputIt.Set( NumericTraits< OutputPixelType >::Zero );
      }
    else
      {
      outputIt.Set( gradientPixel * beamDirection );
      }
    progress.CompletedPixel();
    }
}


template< typename TInputImage, typename TOutputImage, typename TOperatorValue >
void
GradientBasedAngleOfIncidenceImageFilter< TInputImage,
  TOutputImage,
  TOperatorValue >
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "GradientMagnitudeTolerance: "
     << this->m_GradientMagnitudeTolerance
     << indent << "UltrasoundProbeType: ";
  switch( this->m_UltrasoundProbeType )
    {
    case CURVILINEAR:
      os << "CURVILINEAR";
      break;
    case PHASED:
      os << "PHASED";
      break;
    case LINEAR:
      os << "LINEAR";
      break;
    default:
      os << "INVALID";
    }
  os << indent << "UltrasoundProbeOrigin: "
     << this->m_UltrasoundProbeOrigin
     << indent << "BeamDirection: "
     << this->m_UltrasoundProbeBeamDirection
     << std::endl;
}

} // End namespace itk

#endif // End !defined( __itkGradientBasedAngleOfIncidenceImageFilter_hxx )
