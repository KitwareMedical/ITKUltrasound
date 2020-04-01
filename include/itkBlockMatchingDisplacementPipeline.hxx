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
#ifndef itkBlockMatchingDisplacementPipeline_hxx
#define itkBlockMatchingDisplacementPipeline_hxx

#include "itkBlockMatchingDisplacementPipeline.h"

namespace itk
{
namespace BlockMatching
{

template< typename TFixedPixel, typename TMovingPixel,
          typename TMetricPixel, typename TCoordRep,
          unsigned int VImageDimension >
DisplacementPipeline< TFixedPixel, TMovingPixel, TMetricPixel, TCoordRep, VImageDimension >
::DisplacementPipeline():
  m_LevelRegistrationMethodTextProgressBar( false ),
  m_Direction( 0 ),
  m_MaximumAbsStrainAllowed( 0.075 ),
  m_BlockOverlap( 0.75 ),
  m_ScaleBlockByStrain( true ),
  m_RegularizationMaximumNumberOfIterations( 2 )
{
  this->SetNumberOfRequiredInputs( 2 );

  m_FixedResampler = FixedResamplerType::New();
  m_MovingResampler = MovingResamplerType::New();
  m_FixedResamplerInterpolator = FixedResamplerInterpolatorType::New();
  m_MovingResamplerInterpolator = MovingResamplerInterpolatorType::New();
  m_FixedResampler->SetInterpolator( m_FixedResamplerInterpolator );
  m_MovingResampler->SetInterpolator( m_MovingResamplerInterpolator );

  m_BlockRadiusCalculator = BlockRadiusCalculatorType::New();

  m_SearchRegionImageSource = SearchRegionImageSourceType::New();

  m_LevelRegistrationMethod = LevelRegistrationMethodType::New();
  m_TextProgressBar    = TextProgressBarCommand::New();

  m_ParabolicInterpolator = ParabolicInterpolatorType::New();
  m_MaximumPixelInterpolator = MaximumPixelDisplacementCalculatorType::New();
  m_FinalInterpolator = FinalInterpolatorType::New();
  // Optimizing interpolator specific stuff
  m_SubsampleInterpolator = SubsampleInterpolatorType::New();
  m_FinalInterpolator->SetInterpolator( m_SubsampleInterpolator );
  m_SubsampleOptimizer = SubsampleOptimizerType::New();
  typename SubsampleOptimizerType::ParametersType simplexDelta( ImageDimension );
  simplexDelta.Fill( 0.3 );
  m_SubsampleOptimizer->AutomaticInitialSimplexOff();
  m_SubsampleOptimizer->SetInitialSimplexDelta( simplexDelta );
  m_SubsampleOptimizer->SetMaximumNumberOfIterations( 250 );
  m_SubsampleOptimizer->SetParametersConvergenceTolerance( 1.0e-5 );
  m_SubsampleOptimizer->SetFunctionConvergenceTolerance( 10.0 );
  m_FinalInterpolator->SetOptimizer( m_SubsampleOptimizer );

  m_StrainWindower = StrainWindowDisplacementCalculatorType::New();
  m_StrainWindower->SetMaximumIterations( 2 );
  m_StrainWindower->SetDisplacementCalculator( m_ParabolicInterpolator );
  m_StrainWindowStrainFilter = StrainWindowStrainFilterType::New();
  m_HigherOrderAccurateGradientFilter = HigherOrderAccurateGradientFilterType::New();
  m_HigherOrderAccurateGradientFilter->SetOrderOfAccuracy( 2 );
  m_LinearLeastSquaresGradientFilter = LinearLeastSquaresGradientFilterType::New();
  m_LinearLeastSquaresGradientFilter->SetRadius( 2 );
  //m_StrainWindowStrainFilter->SetGradientFilter( m_HigherOrderAccurateGradientFilter );
  m_StrainWindowStrainFilter->SetGradientFilter( m_LinearLeastSquaresGradientFilter );
  m_StrainWindower->SetStrainImageFilter( m_StrainWindowStrainFilter.GetPointer() );

  m_MetricImageFilter = MetricImageFilterType::New();

  m_BlockTransformMetricImageFilter = BlockTransformMetricImageFilterType::New();
  m_BlockTransformMetricImageFilter->SetMetricImageFilter( m_MetricImageFilter );
  m_BlockTransformCommand = BlockTransformCommandType::New();
  m_BlockTransformCommand->SetBlockAffineTransformMetricImageFilter( m_BlockTransformMetricImageFilter );
  m_StrainWindower->AddObserver( itk::EndEvent(), m_BlockTransformCommand );

  m_Regularizer = DisplacmentRegularizerType::New();
  m_Regularizer->SetMetricLowerBound( -1.0 );
  m_Regularizer->SetDisplacementCalculator( m_StrainWindower );

  m_MultiResolutionRegistrationMethod = RegistrationMethodType::New();
  m_MultiResolutionRegistrationMethod->SetBlockRadiusCalculator( m_BlockRadiusCalculator );
  m_MultiResolutionRegistrationMethod->SetSearchRegionImageSource( m_SearchRegionImageSource );
  m_MultiResolutionRegistrationMethod->SetImageRegistrationMethod( m_LevelRegistrationMethod );
  m_DisplacementCalculatorCommand = DisplacementCalculatorCommandType::New();
  m_DisplacementCalculatorCommand->SetLevel0ToNMinus1DisplacementCalculator( m_StrainWindower );
  m_DisplacementCalculatorCommand->SetLevelNDisplacementCalculator( m_FinalInterpolator );
  m_DisplacementCalculatorCommand->SetRegularizer( m_Regularizer );

  m_UpsamplingRatio[0] = 2.0;
  m_UpsamplingRatio[1] = 2.0;

  m_TopBlockRadius[0] = 15;
  m_TopBlockRadius[1] = 10;

  m_BottomBlockRadius[0] = 12;
  m_BottomBlockRadius[1] = 7;

  m_SearchRegionTopFactor[0] = 2.2;
  m_SearchRegionTopFactor[1] = 1.4;

  m_SearchRegionBottomFactor[0] = 1.1;
  m_SearchRegionBottomFactor[1] = 1.1;

  m_RegularizationStrainSigma[0] = 0.075;
  m_RegularizationStrainSigma[1] = 0.15;
}


template< typename TFixedPixel, class TMovingPixel,
          typename TMetricPixel, class TCoordRep,
          unsigned int VImageDimension >
void
DisplacementPipeline< TFixedPixel, TMovingPixel, TMetricPixel, TCoordRep, VImageDimension >
::GenerateOutputInformation()
{
  this->SetupPipeline();

  m_MultiResolutionRegistrationMethod->UpdateOutputInformation();

  DisplacementImageType* output = this->GetOutput( 0 );
  output->CopyInformation( m_MultiResolutionRegistrationMethod->GetOutput( 0 ) );
}


template< typename TFixedPixel, typename TMovingPixel,
          typename TMetricPixel, typename TCoordRep,
          unsigned int VImageDimension >
void
DisplacementPipeline< TFixedPixel, TMovingPixel, TMetricPixel, TCoordRep, VImageDimension >
::SetupPipeline()
{
  typename FixedImageType::Pointer  fixed  = const_cast< FixedImageType * >(
     static_cast< const FixedImageType * >( this->GetInput( 0 )));
  typename MovingImageType::Pointer moving = const_cast< MovingImageType * >(
      static_cast< const MovingImageType * >( this->GetInput( 1 )));

  if( fixed.GetPointer() == nullptr )
    {
    itkExceptionMacro(<< "Fixed image image pointer is nullptr." );
    }
  if( moving.GetPointer() == nullptr )
    {
    itkExceptionMacro(<< "Moving image image pointer is nullptr." );
    }

  fixed->UpdateOutputInformation();
  moving->UpdateOutputInformation();

  // Upsampling.
  m_FixedResampler->SetInput( fixed );
  m_FixedResampler->SetOutputOrigin( fixed->GetOrigin() );
  m_FixedResampler->SetOutputDirection( fixed->GetDirection() );
  m_FixedResampler->SetOutputStartIndex( fixed->GetLargestPossibleRegion().GetIndex() );

  typename FixedImageType::SizeType    size;
  typename FixedImageType::SpacingType spacing;
  size[0] = static_cast< typename FixedImageType::SizeType::SizeValueType >( fixed->GetLargestPossibleRegion().GetSize()[0] * m_UpsamplingRatio[0] );
  size[1] = static_cast< typename FixedImageType::SizeType::SizeValueType >( fixed->GetLargestPossibleRegion().GetSize()[1] * m_UpsamplingRatio[1] );
  spacing[0] = fixed->GetSpacing()[0] / m_UpsamplingRatio[0];
  spacing[1] = fixed->GetSpacing()[1] / m_UpsamplingRatio[1];
  m_FixedResampler->SetOutputSpacing( spacing );
  m_FixedResampler->SetSize( size );
  m_MovingResampler->SetInput( moving );
  m_MovingResampler->SetOutputOrigin( moving->GetOrigin() );
  m_MovingResampler->SetOutputDirection( moving->GetDirection() );
  m_MovingResampler->SetOutputStartIndex( moving->GetLargestPossibleRegion().GetIndex() );
  size[0] = static_cast< typename MovingImageType::SizeType::SizeValueType >( moving->GetLargestPossibleRegion().GetSize()[0] * m_UpsamplingRatio[0] );
  size[1] = static_cast< typename MovingImageType::SizeType::SizeValueType >( moving->GetLargestPossibleRegion().GetSize()[1] * m_UpsamplingRatio[1] );
  spacing[0] = moving->GetSpacing()[0] / m_UpsamplingRatio[0];
  spacing[1] = moving->GetSpacing()[1] / m_UpsamplingRatio[1];
  m_MovingResampler->SetOutputSpacing( spacing );
  m_MovingResampler->SetSize( size );

  // Block Radius Calculator
  RadiusType                         minBlockRadius;
  RadiusType                         maxBlockRadius;
  minBlockRadius[0] = m_BottomBlockRadius[0];
  minBlockRadius[1] = m_BottomBlockRadius[1];
  maxBlockRadius[0] = m_TopBlockRadius[0];
  maxBlockRadius[1] = m_TopBlockRadius[1];
  m_BlockRadiusCalculator->SetMinRadius( minBlockRadius );
  m_BlockRadiusCalculator->SetMaxRadius( maxBlockRadius );

  // Search Region Image Source
  m_SearchRegionImageSource->SetMaxFactor( m_SearchRegionTopFactor );
  m_SearchRegionImageSource->SetMinFactor( m_SearchRegionBottomFactor );

  typename SearchRegionImageSourceType::PyramidScheduleType pyramidSchedule( 3, ImageDimension );
  if( m_Direction == 1 )
    {
    pyramidSchedule( 0, 0 ) = 2;
    pyramidSchedule( 0, 1 ) = 3;
    pyramidSchedule( 1, 0 ) = 1;
    pyramidSchedule( 1, 1 ) = 2;
    pyramidSchedule( 2, 0 ) = 1;
    pyramidSchedule( 2, 1 ) = 1;
    }
  else
    {
    pyramidSchedule( 0, 0 ) = 3;
    pyramidSchedule( 0, 1 ) = 2;
    pyramidSchedule( 1, 0 ) = 2;
    pyramidSchedule( 1, 1 ) = 1;
    pyramidSchedule( 2, 0 ) = 1;
    pyramidSchedule( 2, 1 ) = 1;
    }
  m_SearchRegionImageSource->SetPyramidSchedule( pyramidSchedule );
  m_SearchRegionImageSource->SetOverlapSchedule( m_BlockOverlap );

  // The registration method.
  m_LevelRegistrationMethod->RemoveAllObservers();
  if( m_LevelRegistrationMethodTextProgressBar )
    {
    m_LevelRegistrationMethod->AddObserver( itk::ProgressEvent(), m_TextProgressBar );
    }

  // Filter out peak hopping.
  typedef typename StrainWindowDisplacementCalculatorType::StrainTensorType StrainTensorType;
  StrainTensorType maxStrain;
  maxStrain.Fill( m_MaximumAbsStrainAllowed );
  m_StrainWindower->SetMaximumAbsStrain( maxStrain );

  // Scale the fixed block by the strain at higher levels.
  // Initialize to nullptr because there is initially no previous strain at the top level of the pyramid.
  m_BlockTransformMetricImageFilter->SetStrainImage( nullptr );
  if( m_ScaleBlockByStrain )
    {
    m_LevelRegistrationMethod->SetMetricImageFilter( m_BlockTransformMetricImageFilter );
    }
  else
    {
    m_LevelRegistrationMethod->SetMetricImageFilter( m_MetricImageFilter );
    }

  // Perform regularization.
  m_Regularizer->SetStrainSigma( m_RegularizationStrainSigma );
  m_Regularizer->SetMaximumIterations( m_RegularizationMaximumNumberOfIterations );
// @todo re-enable the ability to use this point examination code.
  // typedef itk::DisplacementRegularizationIterationCommand<
  // DisplacmentRegularizerType >
  // RegularizerCommandType;
  // RegularizerCommandType::Pointer regularizerObserver =
  // RegularizerCommandType::New();
  // regularizerObserver->SetOutputFilePrefix( args.outputPrefix );
  // MetricImageType::PointType targetPoint;
  // targetPoint[0] = args.targetPointAxial;
  // targetPoint[1] = args.targetPointLateral;
  // regularizerObserver->SetTargetPoint( targetPoint );
  // typedef itk::DisplacementRegularizationIterationCommand<
  // DisplacmentRegularizerType >
  // RegularizerCommandType;
  // RegularizerCommandType::Pointer regularizerObserver =
  // RegularizerCommandType::New();
  // regularizerObserver->SetOutputFilePrefix( args.outputPrefix );
  // MetricImageType::PointType targetPoint;
  // targetPoint[0] = args.targetPointAxial;
  // targetPoint[1] = args.targetPointLateral;
  // regularizerObserver->SetTargetPoint( targetPoint );
  // regularizerObserver->SetRegularizer( regularizer );
  // regularizerObserver->SetTruthFile( args.truthImage );
  // regularizer->AddObserver( itk::IterationEvent(), regularizerObserver );
  // regularizer->SetMeanChangeThreshold( 1.0e-25 );
  // regularizer->SetDisplacementCalculator( interpolator );

  if( m_UpsamplingRatio[0] == 1.0 && m_UpsamplingRatio[1] == 1.0 )
    {
    m_MultiResolutionRegistrationMethod->SetFixedImage( fixed );
    m_MultiResolutionRegistrationMethod->SetMovingImage( moving );
    }
  else
    {
    m_MultiResolutionRegistrationMethod->SetFixedImage( m_FixedResampler->GetOutput() );
    m_MultiResolutionRegistrationMethod->SetMovingImage( m_MovingResampler->GetOutput() );
    }
  m_MultiResolutionRegistrationMethod->SetSchedules( pyramidSchedule, pyramidSchedule );
}


template< typename TFixedPixel, typename TMovingPixel,
          typename TMetricPixel, typename TCoordRep,
          unsigned int VImageDimension >
void
DisplacementPipeline< TFixedPixel, TMovingPixel, TMetricPixel, TCoordRep, VImageDimension >
::GenerateData()
{
  this->AllocateOutputs();

  // Set the displacement calculator and regularizer iterations at every level.
  m_MultiResolutionRegistrationMethod->GetImageRegistrationMethod()->SetMetricImageToDisplacementCalculator( m_Regularizer );
  m_DisplacementCalculatorCommand->SetLevel0ToNMinus1RegularizerIterations( 0 );
  m_DisplacementCalculatorCommand->SetLevelNRegularizerIterations( m_RegularizationMaximumNumberOfIterations );
  m_DisplacementCalculatorCommand->SetMultiResolutionMethod( m_MultiResolutionRegistrationMethod );
  m_MultiResolutionRegistrationMethod->AddObserver( itk::IterationEvent(), m_DisplacementCalculatorCommand );
  m_MultiResolutionRegistrationMethod->GraftOutput( this->GetOutput() );
  m_MultiResolutionRegistrationMethod->Update();
  this->GraftOutput( m_MultiResolutionRegistrationMethod->GetOutput() );
}

} // end namespace BlockMatching
} // end namespace itk

#endif
