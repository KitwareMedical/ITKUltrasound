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
  m_LevelRegistrationMethodTextProgressBar( false )
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
  m_RegistrationObserver    = TextProgressBarCommand::New();

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
  m_StrainWindowStrainFilter->SetGradientFilter( m_HigherOrderAccurateGradientFilter );
  m_StrainWindower->SetStrainImageFilter( m_StrainWindowStrainFilter );

  m_MetricImageFilter = MetricImageFilterType::New();

  m_BlockTransformMetricImageFilter = BlockTransformMetricImageFilterType::New();
  m_BlockTransformMetricImageFilter->SetMetricImageFilter( m_MetricImageFilter );
  m_BlockTransformCommand = BlockTransformCommandType::New();
  m_BlockTransformCommand->SetBlockAffineTransformMetricImageFilter( m_BlockTransformMetricImageFilter );
  m_StrainWindower->AddObserver( itk::EndEvent(), m_BlockTransformCommand );

  m_Regularizer = DisplacmentRegularizerType::New();
  m_Regularizer->SetMetricLowerBound( -1.0 );
  m_Regularizer->SetDisplacementCalculator( m_StrainWindower );

  m_MultiResRegistrationMethod = RegistrationMethodType::New();
  m_MultiResRegistrationMethod->SetBlockRadiusCalculator( m_BlockRadiusCalculator );
  m_MultiResRegistrationMethod->SetSearchRegionImageSource( m_SearchRegionImageSource );
  m_MultiResRegistrationMethod->SetImageRegistrationMethod( m_LevelRegistrationMethod );
  m_MultiResObserver = MultiResolutionObserverType::New();
  m_MultiResObserver->SetMultiResolutionMethod( m_MultiResRegistrationMethod );
  m_DisplacementCalculatorCommand = DisplacementCalculatorCommandType::New();
  m_DisplacementCalculatorCommand->SetLevel0ToNMinus1DisplacementCalculator( m_StrainWindower );
  m_DisplacementCalculatorCommand->SetLevelNDisplacementCalculator( m_FinalInterpolator );
  m_DisplacementCalculatorCommand->SetRegularizer( m_Regularizer );

  m_SearchRegionWriterCommand = SearchRegionWriterCommandType::New();
  m_DisplacementWriter = DisplacementWriterType::New();
  m_StrainFilter = StrainFilterType::New();
  m_StrainFilter->SetStrainForm( StrainFilterType::EULERIANALMANSI );
  m_FinalGradientFilter = FinalGradientFilterType::New();
  m_FinalGradientFilter->SetNumberOfLevels( 1 );

  m_UpsamplingRatio[0] = 1.0;
  m_UpsamplingRatio[1] = 2.0;

  m_TopBlockRadius[0] = 15;
  m_TopBlockRadius[1] = 10;

  m_BottomBlockRadius[0] = 12;
  m_BottomBlockRadius[1] = 7;
}


template< typename TFixedPixel, class TMovingPixel,
          typename TMetricPixel, class TCoordRep,
          unsigned int VImageDimension >
void
DisplacementPipeline< TFixedPixel, TMovingPixel, TMetricPixel, TCoordRep, VImageDimension >
::GenerateOutputInformation()
{
  this->SetupPipeline();

  m_MultiResRegistrationMethod->UpdateOutputInformation();

  DisplacementImageType* output = this->GetOutput( 0 );
  if( !output )
    {
    return;
    }
  output->CopyInformation( m_MultiResRegistrationMethod->GetOutput( 0 ) );
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

  if( fixed.GetPointer() == NULL )
    {
    itkExceptionMacro(<< "Fixed image image pointer is NULL." );
    }
  if( moving.GetPointer() == NULL )
    {
    itkExceptionMacro(<< "Moving image image pointer is NULL." );
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
  size[0] = static_cast< typename FixedImageType::SizeType::SizeValueType >( fixed->GetLargestPossibleRegion().GetSize()[0] * m_UpsamplingRatio[0];
  size[1] = static_cast< typename FixedImageType::SizeType::SizeValueType >( fixed->GetLargestPossibleRegion().GetSize()[1] * m_UpsamplingRatio[1];
  spacing[0] = fixed->GetSpacing()[0] / m_UpsamplingRatio[0];
  spacing[1] = fixed->GetSpacing()[1] / m_UpsamplingRatio[1];
  m_FixedResampler->SetOutputSpacing( spacing );
  m_FixedResampler->SetSize( size );
  m_MovingResampler->SetInput( moving );
  m_MovingResampler->SetOutputOrigin( moving->GetOrigin() );
  m_MovingResampler->SetOutputDirection( moving->GetDirection() );
  m_MovingResampler->SetOutputStartIndex( moving->GetLargestPossibleRegion().GetIndex() );
  size[0] = static_cast< typename MovingImageType::SizeType::SizeValueType >( moving->GetLargestPossibleRegion().GetSize()[0] * m_UpsamplingRatio[0];
  size[1] = static_cast< typename MovingImageType::SizeType::SizeValueType >( moving->GetLargestPossibleRegion().GetSize()[1] * m_UpsamplingRatio[1];
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
  typename SearchRegionImageSourceType::FactorType topFactor;
  typename SearchRegionImageSourceType::FactorType bottomFactor;
  topFactor[0] = m_MMMStrainOptions->GetParameters()->searchRegion.topFactor[0];
  topFactor[1] = m_MMMStrainOptions->GetParameters()->searchRegion.topFactor[1];
  m_SearchRegionImageSource->SetMaxFactor( topFactor );
  bottomFactor[0] = m_MMMStrainOptions->GetParameters()->searchRegion.bottomFactor[0];
  bottomFactor[1] = m_MMMStrainOptions->GetParameters()->searchRegion.bottomFactor[1];
  m_SearchRegionImageSource->SetMinFactor( bottomFactor );

  typename SearchRegionImageSourceType::PyramidScheduleType pyramidSchedule( 3, ImageDimension );
  if( m_MMMStrainOptions->GetParameters()->axialDirection == 1 )
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
  m_SearchRegionImageSource->SetOverlapSchedule( m_MMMStrainOptions->GetParameters()->block.blockOverlap );

  // The registration method.
  m_LevelRegistrationMethod->RemoveAllObservers();
  if( m_LevelRegistrationMethodTextProgressBar )
    {
    m_LevelRegistrationMethod->AddObserver( itk::ProgressEvent(), m_RegistrationObserver );
    }

  // Filter out peak hopping.
  typedef typename StrainWindowDisplacementCalculatorType::StrainTensorType StrainTensorType;
  StrainTensorType maxStrain;
  maxStrain.Fill( m_MMMStrainOptions->GetParameters()->maximumAbsStrainAllowed );
  m_StrainWindower->SetMaximumAbsStrain( maxStrain );

  // Scale the fixed block by the strain at higher levels.
  // Initialize to NULL because there is initially no previous strain at the top level of the pyramid.
  m_BlockTransformMetricImageFilter->SetStrainImage( NULL );
  if( m_MMMStrainOptions->GetParameters()->block.scaleByStrain )
    {
    m_LevelRegistrationMethod->SetMetricImageFilter( m_BlockTransformMetricImageFilter );
    }
  else
    {
    m_LevelRegistrationMethod->SetMetricImageFilter( m_MetricImageFilter );
    }

  // Perform regularization.
  typename MetricImageType::SpacingType strainSigma;
  strainSigma[0] = m_MMMStrainOptions->GetParameters()->regularization.strainSigma[0];
  strainSigma[1] = m_MMMStrainOptions->GetParameters()->regularization.strainSigma[1];
  m_Regularizer->SetStrainSigma( strainSigma );
  m_Regularizer->SetMaximumIterations( m_MMMStrainOptions->GetParameters()->regularization.maximumIterations );
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

  if( m_MMMStrainOptions->GetParameters()->upsample[0] == 1.0 && m_MMMStrainOptions->GetParameters()->upsample[1] == 1.0 )
    {
    m_MultiResRegistrationMethod->SetFixedImage( fixed );
    m_MultiResRegistrationMethod->SetMovingImage( moving );
    }
  else
    {
    m_MultiResRegistrationMethod->SetFixedImage( m_FixedResampler->GetOutput() );
    m_MultiResRegistrationMethod->SetMovingImage( m_MovingResampler->GetOutput() );
    }
  m_MultiResRegistrationMethod->SetSchedules( pyramidSchedule, pyramidSchedule );
  m_MultiResObserver->SetOutputFilePrefix( m_MMMStrainOptions->GetFiles()->outputPrefix );
  m_MultiResRegistrationMethod->RemoveAllObservers();
}

template< typename TFixedPixel, typename TMovingPixel,
          typename TMetricPixel, typename TCoordRep,
          unsigned int VImageDimension >
void
DisplacementPipeline< TFixedPixel, TMovingPixel, TMetricPixel, TCoordRep, VImageDimension >
::GenerateData()
{
  this->AllocateOutputs();

  if( m_WriteOutputImagesToFile )
    {
    m_MultiResRegistrationMethod->AddObserver( itk::IterationEvent(), m_MultiResObserver );
    }

  // Set the displacement calculator and regularizer iterations at every level.
  m_MultiResRegistrationMethod->GetImageRegistrationMethod()->SetMetricImageToDisplacementCalculator( m_Regularizer );
  //if ( args.maximumIterations == 0 )
    //{
    //displacementCalculatorCommand->SetLevel0ToNMinus1RegularizerIterations( 0 );
    //}
  //else
    //{
    m_DisplacementCalculatorCommand->SetLevel0ToNMinus1RegularizerIterations( 3 );
    //}
  m_DisplacementCalculatorCommand->SetLevelNRegularizerIterations( m_MMMStrainOptions->GetParameters()->regularization.maximumIterations );
  m_DisplacementCalculatorCommand->SetMultiResolutionMethod( m_MultiResRegistrationMethod );
  m_MultiResRegistrationMethod->AddObserver( itk::IterationEvent(), m_DisplacementCalculatorCommand );
  m_MultiResRegistrationMethod->GraftOutput( this->GetOutput() );

  // Write out the search region images at every level.
  if( m_WriteOutputImagesToFile )
    {
    m_SearchRegionWriterCommand->SetOutputFilePrefix( m_MMMStrainOptions->GetFiles()->outputPrefix );
    m_SearchRegionWriterCommand->SetMultiResolutionMethod( m_MultiResRegistrationMethod );
    if( m_MMMStrainOptions->GetVerbosity()->writeSearchRegionImages )
      {
      m_MultiResRegistrationMethod->AddObserver( itk::IterationEvent(), m_SearchRegionWriterCommand );
      }

    // Write displacement vector.
    m_DisplacementWriter->SetFileName( m_MMMStrainOptions->GetFiles()->outputPrefix + "_DisplacementVectors.mha" );
    m_DisplacementWriter->SetInput( m_MultiResRegistrationMethod->GetOutput() );
    m_DisplacementWriter->Update();

    // Calculate strains.
    m_StrainFilter->SetInput( m_MultiResRegistrationMethod->GetOutput() );
    m_FinalGradientFilter->SetControlPointSpacingRatio( m_MMMStrainOptions->GetParameters()->controlPointRatio );
    if( m_MMMStrainOptions->GetParameters()->controlPointRatio == 1.0 )
      {
      //m_StrainFilter->SetGradientFilter( m_HigherOrderAccurateGradientFilter );
      m_StrainFilter->SetGradientFilter( m_LinearLeastSquaresGradientFilter );
      }
    else
      {
      m_StrainFilter->SetVectorGradientFilter( m_FinalGradientFilter );
      }

    // Write out strains.
    m_StrainWriter->SetFileName( m_MMMStrainOptions->GetFiles()->outputPrefix + "_StrainTensors.vtk" );
    m_StrainWriter->Update();

    // Write out strain components.
    std::ostringstream ostr;
    for( unsigned int i = 0; i < 3; i++ )
      {
      m_StrainComponentWriter->SetInput( m_StrainComponentsFilter->GetOutput( i ) );
      ostr.str( "" );
      ostr << "_StrainComponent" << i << ".mha";
      m_StrainComponentWriter->SetFileName( m_MMMStrainOptions->GetFiles()->outputPrefix + ostr.str() );
      m_StrainComponentWriter->Update();
      }

    if( m_MMMStrainOptions->GetVerbosity()->writeDisplacementVectorComponents )
      {
      // Write displacement vector components.
      m_DisplacementComponentsFilter->SetInput( m_MultiResRegistrationMethod->GetOutput() );
      m_DisplacementComponentsFilter->Update();
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        m_DisplacementComponentWriter->SetInput( m_DisplacementComponentsFilter->GetOutput( i ) );
        ostr.str( "" );
        ostr << "_DisplacementComponent" << i << ".mha";
        m_DisplacementComponentWriter->SetFileName( m_MMMStrainOptions->GetFiles()->outputPrefix + ostr.str() );
        m_DisplacementComponentWriter->Update();
        }
      }
    } // end if write output images
  else
    {
    m_MultiResRegistrationMethod->Update();
    }
  this->GraftOutput( m_MultiResRegistrationMethod->GetOutput() );
}

} // end namespace BlockMatching
} // end namespace itk
