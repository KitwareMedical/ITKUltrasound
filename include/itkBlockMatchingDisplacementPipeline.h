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
#ifndef itkBlockMatchingDisplacementPipeline_h
#define itkBlockMatchingDisplacementPipeline_h

#include "itkAmoebaOptimizer.h"
#include "itkExpNegativeImageFilter.h"
#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include "itkVector.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"

#include "itkBlockMatchingBayesianRegularizationDisplacementCalculator.h"
#include "itkBlockMatchingBlockAffineTransformMetricImageFilter.h"
#include "itkBlockMatchingImageRegistrationMethod.h"
#include "itkBlockMatchingMaximumPixelDisplacementCalculator.h"
#include "itkBlockMatchingMultiResolutionImageRegistrationMethod.h"
#include "itkBlockMatchingMultiResolutionIterationDisplacementCalculatorCommand.h"
#include "itkBlockMatchingMultiResolutionMinMaxBlockRadiusCalculator.h"
#include "itkBlockMatchingMultiResolutionMinMaxSearchRegionImageSource.h"
#include "itkMultiResolutionImageRegistrationMethod.h"
#include "itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter.h"
#include "itkBlockMatchingOptimizingInterpolationDisplacementCalculator.h"
#include "itkBlockMatchingParabolicInterpolationDisplacementCalculator.h"
#include "itkBlockMatchingStrainWindowDisplacementCalculator.h"
#include "itkBlockMatchingStrainWindowBlockAffineTransformCommand.h"
#include "itkTextProgressBarCommand.h"

#include "itkBSplineApproximationGradientImageFilter.h"
#include "itkLinearLeastSquaresGradientImageFilter.h"
#include "itkHigherOrderAccurateGradientImageFilter.h"
#include "itkStrainImageFilter.h"


namespace itk
{
namespace BlockMatching
{

/** \class DisplacementPipeline
 *
 * \brief Sets up and runs deformable image registration pipeline with block-matching.
 *
 * \ingroup Ultrasound
 */
template< typename TFixedPixel = signed short, typename TMovingPixel = TFixedPixel,
          typename TMetricPixel = double, typename TCoordRep = double,
          unsigned int VImageDimension = 2 >
class ITK_TEMPLATE_EXPORT DisplacementPipeline : public ImageToImageFilter<
  Image< TFixedPixel, VImageDimension >, Image< Vector< TMetricPixel, VImageDimension>, VImageDimension > >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(DisplacementPipeline);

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int, VImageDimension);

  typedef TFixedPixel                              FixedPixelType;
  typedef Image< TFixedPixel, ImageDimension >     FixedImageType;
  typedef typename FixedImageType::SizeType        RadiusType;

  typedef TMovingPixel                             MovingPixelType;
  typedef Image< MovingPixelType, ImageDimension > MovingImageType;

  typedef TMetricPixel                             MetricPixelType;
  typedef Image< MetricPixelType, ImageDimension > MetricImageType;

  typedef Vector<MetricPixelType, ImageDimension > VectorType;
  typedef Image<VectorType, ImageDimension >       DisplacementImageType;
  typedef DisplacementImageType                    OutputImageType;

  /** Standard class typedefs. */
  typedef DisplacementPipeline                                        Self;
  typedef ImageToImageFilter< FixedImageType, DisplacementImageType > Superclass;
  typedef SmartPointer< Self >                                        Pointer;
  typedef SmartPointer< const Self >                                  ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( DisplacementPipeline, ImageToImageFilter );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  typedef TCoordRep CoordRepType;

  typedef ResampleImageFilter< FixedImageType, FixedImageType, CoordRepType >   FixedResamplerType;
  typedef ResampleImageFilter< MovingImageType, MovingImageType, CoordRepType > MovingResamplerType;

  const static unsigned int RESAMPLE_RADIUS = 4;
  typedef Function::WelchWindowFunction< RESAMPLE_RADIUS >    ResampleWindowType;
  typedef ZeroFluxNeumannBoundaryCondition< FixedImageType >  FixedBoundaryConditionType;
  typedef WindowedSincInterpolateImageFunction< FixedImageType,
                                                    RESAMPLE_RADIUS, ResampleWindowType, FixedBoundaryConditionType,
                                                    CoordRepType> FixedResamplerInterpolatorType;
  typedef ZeroFluxNeumannBoundaryCondition< MovingImageType >  MovingBoundaryConditionType;
  typedef WindowedSincInterpolateImageFunction< MovingImageType,
                                                    RESAMPLE_RADIUS, ResampleWindowType, MovingBoundaryConditionType,
                                                    CoordRepType> MovingResamplerInterpolatorType;

  /** The block radius calculator. */
  typedef BlockMatching::MultiResolutionMinMaxBlockRadiusCalculator< FixedImageType >
    BlockRadiusCalculatorType;

  /** The search region image source. */
  typedef BlockMatching::MultiResolutionMinMaxSearchRegionImageSource< FixedImageType,
          MovingImageType, DisplacementImageType > SearchRegionImageSourceType;

  /** The registration method. */
  typedef BlockMatching::ImageRegistrationMethod<FixedImageType,
                                                      MovingImageType, MetricImageType, DisplacementImageType,
                                                      CoordRepType>
    LevelRegistrationMethodType;

  /** Interpolation classes. */
  typedef BlockMatching::ParabolicInterpolationDisplacementCalculator<
    MetricImageType,
    DisplacementImageType> ParabolicInterpolatorType;
  typedef BlockMatching::MaximumPixelDisplacementCalculator<MetricImageType,
                                                                 DisplacementImageType>
  MaximumPixelDisplacementCalculatorType;
  typedef BlockMatching::OptimizingInterpolationDisplacementCalculator<MetricImageType,
                                                                            DisplacementImageType> FinalInterpolatorType;
  const static unsigned int OPTIMIZING_INTERPOLATOR_RADIUS = 4;
  typedef Function::WelchWindowFunction< OPTIMIZING_INTERPOLATOR_RADIUS >  WindowType;
  typedef ZeroFluxNeumannBoundaryCondition< MetricImageType >              ResampleBoundaryConditionType;
  typedef WindowedSincInterpolateImageFunction<MetricImageType,
                                                    OPTIMIZING_INTERPOLATOR_RADIUS, WindowType, ResampleBoundaryConditionType,
                                                    CoordRepType> SubsampleInterpolatorType;
  typedef AmoebaOptimizer SubsampleOptimizerType;

  /** Filter out peak hopping. */
  typedef BlockMatching::StrainWindowDisplacementCalculator<MetricImageType, DisplacementImageType, MetricPixelType >
    StrainWindowDisplacementCalculatorType;
  typedef StrainImageFilter<DisplacementImageType, MetricPixelType, MetricPixelType>
  StrainWindowStrainFilterType;
  typedef HigherOrderAccurateGradientImageFilter<MetricImageType, MetricPixelType, MetricPixelType>
  HigherOrderAccurateGradientFilterType;
  typedef LinearLeastSquaresGradientImageFilter< MetricImageType, MetricPixelType, MetricPixelType >
    LinearLeastSquaresGradientFilterType;

  /** The image similarity metric. */
  typedef BlockMatching::NormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter<
    FixedImageType, MovingImageType, MetricImageType> MetricImageFilterType;

  /** Scale the fixed block by the strain at higher levels. */
  typedef BlockMatching::BlockAffineTransformMetricImageFilter< FixedImageType,
          MovingImageType, MetricImageType, MetricPixelType > BlockTransformMetricImageFilterType;
  typedef BlockMatching::StrainWindowBlockAffineTransformCommand< StrainWindowDisplacementCalculatorType,
          BlockTransformMetricImageFilterType, StrainWindowStrainFilterType > BlockTransformCommandType;

  /** Perform regularization. */
  typedef BlockMatching::BayesianRegularizationDisplacementCalculator<
    MetricImageType, DisplacementImageType> DisplacmentRegularizerType;

  /** Multi-resolution registration method. */
  typedef BlockMatching::MultiResolutionImageRegistrationMethod< FixedImageType, MovingImageType, MetricImageType,
                                                                     DisplacementImageType,
                                                                     CoordRepType> RegistrationMethodType;

  /** Set the displacement calculator and regularizer iterations at every level. */
  typedef BlockMatching::MultiResolutionIterationDisplacementCalculatorCommand<RegistrationMethodType>
    DisplacementCalculatorCommandType;

  /** Calculate strains. */
  typedef StrainImageFilter<DisplacementImageType, MetricPixelType, MetricPixelType> StrainFilterType;
  typedef typename StrainFilterType::OutputImageType                                 TensorImageType;

  typedef FixedArray< double, ImageDimension > UpsamplingRatioType;

  typedef FixedArray< SizeValueType, ImageDimension > BlockRadiusType;

  typedef typename SearchRegionImageSourceType::FactorType SearchRegionFactorType;

  typedef typename MetricImageType::SpacingType RegularizationStrainSigmaType;

  /** Set the fixed image. */
  void SetFixedImage( FixedImageType * fixed )
    {
    this->SetNthInput( 0, fixed );
    }
  const FixedImageType * GetFixedImage()
    {
    return static_cast< const FixedImageType * >( this->ProcessObject::GetInput( 0 ) );
    }

  /** Set the moving image. */
  void SetMovingImage( MovingImageType * moving )
    {
    this->SetNthInput( 1, moving );
    }
  const MovingImageType * GetMovingImage()
    {
    return static_cast< const MovingImageType * >( this->ProcessObject::GetInput( 1 ) );
    }

  /** Any point with a strain component above the given value in the higher levels
   * will have its displacement interpolated by surrounding areas. */
  itkSetMacro( MaximumAbsStrainAllowed, double );
  itkGetConstMacro( MaximumAbsStrainAllowed, double );

  /** Display a text progress bar for the level registration method's creation
   * of the metric images. */
  itkSetMacro( LevelRegistrationMethodTextProgressBar, bool );
  itkGetConstMacro( LevelRegistrationMethodTextProgressBar, bool );

  /** Upsampling ratio applied to the input images */
  itkSetMacro( UpsamplingRatio, UpsamplingRatioType );
  itkGetConstReferenceMacro( UpsamplingRatio, UpsamplingRatioType );


  /** Set the block radius at the top level */
  itkSetMacro( TopBlockRadius, BlockRadiusType );
  itkGetConstReferenceMacro( TopBlockRadius, BlockRadiusType );

  /** Set the block radius at the bottom level */
  itkSetMacro( BottomBlockRadius, BlockRadiusType );
  itkGetConstReferenceMacro( BottomBlockRadius, BlockRadiusType );

  /** Block overlap. 1.0 is no overlap. 0.5 is 50% overlap. */
  itkSetMacro( BlockOverlap, double );
  itkGetConstMacro( BlockOverlap, double );

  /** In the multiresolution method, scale the matching block by the strain
   *   estimated at higher levels. */
  itkSetMacro( ScaleBlockByStrain, bool );
  itkGetConstMacro( ScaleBlockByStrain, bool );


  /** Search region radius at the top level is the following factor times the block radius.
   *  The factors at intermediate levels between the top level and bottom level
   *  are linearly interpolated.  */
  itkSetMacro( SearchRegionTopFactor, SearchRegionFactorType );
  itkGetConstReferenceMacro( SearchRegionTopFactor, SearchRegionFactorType );

  /** Search region radius at the bottom level is the following factor times the block radius. */
  itkSetMacro( SearchRegionBottomFactor, SearchRegionFactorType );
  itkGetConstReferenceMacro( SearchRegionBottomFactor, SearchRegionFactorType );


  /** RF beamline, i.e. axial, direction. */
  itkSetMacro( Direction, unsigned int );
  itkGetConstMacro( Direction, unsigned int );


  /** Strain regularization parameter. */
  itkSetMacro( RegularizationStrainSigma, RegularizationStrainSigmaType );
  itkGetConstReferenceMacro( RegularizationStrainSigma, RegularizationStrainSigmaType );

  /** Maximum number of iterations during regularization at the bottom level. */
  itkSetMacro( RegularizationMaximumNumberOfIterations, unsigned int );
  itkGetConstMacro( RegularizationMaximumNumberOfIterations, unsigned int );

  /** Get the multi-resolution image registration method. */
  itkGetModifiableObjectMacro( MultiResolutionRegistrationMethod, RegistrationMethodType );

protected:
  DisplacementPipeline();

  void GenerateOutputInformation() override;

  virtual void SetupPipeline();

  /** Execute the pipeline. */
  void GenerateData() override;


private:
  typename FixedResamplerType::Pointer              m_FixedResampler;
  typename MovingResamplerType::Pointer             m_MovingResampler;
  typename FixedResamplerInterpolatorType::Pointer  m_FixedResamplerInterpolator;
  typename MovingResamplerInterpolatorType::Pointer m_MovingResamplerInterpolator;

  typename BlockRadiusCalculatorType::Pointer m_BlockRadiusCalculator;

  typename SearchRegionImageSourceType::Pointer m_SearchRegionImageSource;

  typename LevelRegistrationMethodType::Pointer m_LevelRegistrationMethod;
  TextProgressBarCommand::Pointer               m_TextProgressBar;
  bool                                          m_LevelRegistrationMethodTextProgressBar;

  typename ParabolicInterpolatorType::Pointer              m_ParabolicInterpolator;
  typename MaximumPixelDisplacementCalculatorType::Pointer m_MaximumPixelInterpolator;
  typename FinalInterpolatorType::Pointer                  m_FinalInterpolator;
  typename SubsampleInterpolatorType::Pointer              m_SubsampleInterpolator;
  typename SubsampleOptimizerType::Pointer                 m_SubsampleOptimizer;

  typename StrainWindowDisplacementCalculatorType::Pointer m_StrainWindower;
  typename StrainWindowStrainFilterType::Pointer           m_StrainWindowStrainFilter;
  typename HigherOrderAccurateGradientFilterType::Pointer  m_HigherOrderAccurateGradientFilter;
  typename LinearLeastSquaresGradientFilterType::Pointer   m_LinearLeastSquaresGradientFilter;

  typename MetricImageFilterType::Pointer m_MetricImageFilter;

  typename BlockTransformMetricImageFilterType::Pointer m_BlockTransformMetricImageFilter;
  typename BlockTransformCommandType::Pointer           m_BlockTransformCommand;

  typename DisplacmentRegularizerType::Pointer m_Regularizer;

  typename RegistrationMethodType::Pointer            m_MultiResolutionRegistrationMethod;
  typename DisplacementCalculatorCommandType::Pointer m_DisplacementCalculatorCommand;

  UpsamplingRatioType m_UpsamplingRatio;

  BlockRadiusType m_TopBlockRadius;
  BlockRadiusType m_BottomBlockRadius;

  SearchRegionFactorType m_SearchRegionTopFactor;
  SearchRegionFactorType m_SearchRegionBottomFactor;

  unsigned int m_Direction;
  double       m_BlockOverlap;
  double       m_MaximumAbsStrainAllowed;
  bool         m_ScaleBlockByStrain;

  RegularizationStrainSigmaType m_RegularizationStrainSigma;
  unsigned int                  m_RegularizationMaximumNumberOfIterations;

};

} // end namespace BlockMatching
} // end itk namespace

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlockMatchingDisplacementPipeline.hxx"
#endif

#endif
