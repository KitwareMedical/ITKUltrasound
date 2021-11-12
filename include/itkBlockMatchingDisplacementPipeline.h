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
template <typename TFixedPixel = signed short,
          typename TMovingPixel = TFixedPixel,
          typename TMetricPixel = double,
          typename TCoordRep = double,
          unsigned int VImageDimension = 2>
class ITK_TEMPLATE_EXPORT DisplacementPipeline
  : public ImageToImageFilter<Image<TFixedPixel, VImageDimension>,
                              Image<Vector<TMetricPixel, VImageDimension>, VImageDimension>>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(DisplacementPipeline);

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int, VImageDimension);

  using FixedPixelType = TFixedPixel;
  using FixedImageType = Image<TFixedPixel, ImageDimension>;
  using RadiusType = typename FixedImageType::SizeType;

  using MovingPixelType = TMovingPixel;
  using MovingImageType = Image<MovingPixelType, ImageDimension>;

  using MetricPixelType = TMetricPixel;
  using MetricImageType = Image<MetricPixelType, ImageDimension>;

  using VectorType = Vector<MetricPixelType, ImageDimension>;
  using DisplacementImageType = Image<VectorType, ImageDimension>;
  using OutputImageType = DisplacementImageType;

  /** Standard class type alias. */
  using Self = DisplacementPipeline;
  using Superclass = ImageToImageFilter<FixedImageType, DisplacementImageType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information (and related methods). */
  itkTypeMacro(DisplacementPipeline, ImageToImageFilter);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  using CoordRepType = TCoordRep;

  using FixedResamplerType = ResampleImageFilter<FixedImageType, FixedImageType, CoordRepType>;
  using MovingResamplerType = ResampleImageFilter<MovingImageType, MovingImageType, CoordRepType>;

  const static unsigned int RESAMPLE_RADIUS = 4;
  using ResampleWindowType = Function::WelchWindowFunction<RESAMPLE_RADIUS>;
  using FixedBoundaryConditionType = ZeroFluxNeumannBoundaryCondition<FixedImageType>;
  using FixedResamplerInterpolatorType = WindowedSincInterpolateImageFunction<FixedImageType,
                                                                              RESAMPLE_RADIUS,
                                                                              ResampleWindowType,
                                                                              FixedBoundaryConditionType,
                                                                              CoordRepType>;
  using MovingBoundaryConditionType = ZeroFluxNeumannBoundaryCondition<MovingImageType>;
  using MovingResamplerInterpolatorType = WindowedSincInterpolateImageFunction<MovingImageType,
                                                                               RESAMPLE_RADIUS,
                                                                               ResampleWindowType,
                                                                               MovingBoundaryConditionType,
                                                                               CoordRepType>;

  /** The block radius calculator. */
  using BlockRadiusCalculatorType = BlockMatching::MultiResolutionMinMaxBlockRadiusCalculator<FixedImageType>;

  /** The search region image source. */
  using SearchRegionImageSourceType =
    BlockMatching::MultiResolutionMinMaxSearchRegionImageSource<FixedImageType, MovingImageType, DisplacementImageType>;

  /** The registration method. */
  using LevelRegistrationMethodType = BlockMatching::
    BlockMatchingImageRegistrationMethod<FixedImageType, MovingImageType, MetricImageType, DisplacementImageType, CoordRepType>;

  /** Interpolation classes. */
  using ParabolicInterpolatorType =
    BlockMatching::ParabolicInterpolationDisplacementCalculator<MetricImageType, DisplacementImageType>;
  using MaximumPixelDisplacementCalculatorType =
    BlockMatching::MaximumPixelDisplacementCalculator<MetricImageType, DisplacementImageType>;
  using FinalInterpolatorType =
    BlockMatching::OptimizingInterpolationDisplacementCalculator<MetricImageType, DisplacementImageType>;
  const static unsigned int OPTIMIZING_INTERPOLATOR_RADIUS = 4;
  using WindowType = Function::WelchWindowFunction<OPTIMIZING_INTERPOLATOR_RADIUS>;
  using ResampleBoundaryConditionType = ZeroFluxNeumannBoundaryCondition<MetricImageType>;
  using SubsampleInterpolatorType = WindowedSincInterpolateImageFunction<MetricImageType,
                                                                         OPTIMIZING_INTERPOLATOR_RADIUS,
                                                                         WindowType,
                                                                         ResampleBoundaryConditionType,
                                                                         CoordRepType>;
  using SubsampleOptimizerType = AmoebaOptimizer;

  /** Filter out peak hopping. */
  using StrainWindowDisplacementCalculatorType =
    BlockMatching::StrainWindowDisplacementCalculator<MetricImageType, DisplacementImageType, MetricPixelType>;
  using StrainWindowStrainFilterType = StrainImageFilter<DisplacementImageType, MetricPixelType, MetricPixelType>;
  using HigherOrderAccurateGradientFilterType =
    HigherOrderAccurateGradientImageFilter<MetricImageType, MetricPixelType, MetricPixelType>;
  using LinearLeastSquaresGradientFilterType =
    LinearLeastSquaresGradientImageFilter<MetricImageType, MetricPixelType, MetricPixelType>;

  /** The image similarity metric. */
  using MetricImageFilterType = BlockMatching::
    NormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter<FixedImageType, MovingImageType, MetricImageType>;

  /** Scale the fixed block by the strain at higher levels. */
  using BlockTransformMetricImageFilterType = BlockMatching::
    BlockAffineTransformMetricImageFilter<FixedImageType, MovingImageType, MetricImageType, MetricPixelType>;
  using BlockTransformCommandType =
    BlockMatching::StrainWindowBlockAffineTransformCommand<StrainWindowDisplacementCalculatorType,
                                                           BlockTransformMetricImageFilterType,
                                                           StrainWindowStrainFilterType>;

  /** Perform regularization. */
  using DisplacmentRegularizerType =
    BlockMatching::BayesianRegularizationDisplacementCalculator<MetricImageType, DisplacementImageType>;

  /** Multi-resolution registration method. */
  using RegistrationMethodType = BlockMatching::BlockMatchingMultiResolutionImageRegistrationMethod<FixedImageType,
                                                                                       MovingImageType,
                                                                                       MetricImageType,
                                                                                       DisplacementImageType,
                                                                                       CoordRepType>;

  /** Set the displacement calculator and regularizer iterations at every level. */
  using DisplacementCalculatorCommandType =
    BlockMatching::MultiResolutionIterationDisplacementCalculatorCommand<RegistrationMethodType>;

  /** Calculate strains. */
  using StrainFilterType = StrainImageFilter<DisplacementImageType, MetricPixelType, MetricPixelType>;
  using TensorImageType = typename StrainFilterType::OutputImageType;

  using UpsamplingRatioType = FixedArray<double, ImageDimension>;

  using BlockRadiusType = FixedArray<SizeValueType, ImageDimension>;

  using SearchRegionFactorType = typename SearchRegionImageSourceType::FactorType;

  using RegularizationStrainSigmaType = typename MetricImageType::SpacingType;

  /** Set the fixed image. */
  void
  SetFixedImage(FixedImageType * fixed)
  {
    this->SetNthInput(0, fixed);
  }
  const FixedImageType *
  GetFixedImage()
  {
    return static_cast<const FixedImageType *>(this->ProcessObject::GetInput(0));
  }

  /** Set the moving image. */
  void
  SetMovingImage(MovingImageType * moving)
  {
    this->SetNthInput(1, moving);
  }
  const MovingImageType *
  GetMovingImage()
  {
    return static_cast<const MovingImageType *>(this->ProcessObject::GetInput(1));
  }

  /** Any point with a strain component above the given value in the higher levels
   * will have its displacement interpolated by surrounding areas. */
  itkSetMacro(MaximumAbsStrainAllowed, double);
  itkGetConstMacro(MaximumAbsStrainAllowed, double);

  /** Display a text progress bar for the level registration method's creation
   * of the metric images. */
  itkSetMacro(LevelRegistrationMethodTextProgressBar, bool);
  itkGetConstMacro(LevelRegistrationMethodTextProgressBar, bool);

  /** Upsampling ratio applied to the input images */
  itkSetMacro(UpsamplingRatio, UpsamplingRatioType);
  itkGetConstReferenceMacro(UpsamplingRatio, UpsamplingRatioType);


  /** Set the block radius at the top level */
  itkSetMacro(TopBlockRadius, BlockRadiusType);
  itkGetConstReferenceMacro(TopBlockRadius, BlockRadiusType);

  /** Set the block radius at the bottom level */
  itkSetMacro(BottomBlockRadius, BlockRadiusType);
  itkGetConstReferenceMacro(BottomBlockRadius, BlockRadiusType);

  /** Block overlap. 1.0 is no overlap. 0.5 is 50% overlap. */
  itkSetMacro(BlockOverlap, double);
  itkGetConstMacro(BlockOverlap, double);

  /** In the multiresolution method, scale the matching block by the strain
   *   estimated at higher levels. */
  itkSetMacro(ScaleBlockByStrain, bool);
  itkGetConstMacro(ScaleBlockByStrain, bool);


  /** Search region radius at the top level is the following factor times the block radius.
   *  The factors at intermediate levels between the top level and bottom level
   *  are linearly interpolated.  */
  itkSetMacro(SearchRegionTopFactor, SearchRegionFactorType);
  itkGetConstReferenceMacro(SearchRegionTopFactor, SearchRegionFactorType);

  /** Search region radius at the bottom level is the following factor times the block radius. */
  itkSetMacro(SearchRegionBottomFactor, SearchRegionFactorType);
  itkGetConstReferenceMacro(SearchRegionBottomFactor, SearchRegionFactorType);


  /** RF beamline, i.e. axial, direction. */
  itkSetMacro(Direction, unsigned int);
  itkGetConstMacro(Direction, unsigned int);


  /** Strain regularization parameter. */
  itkSetMacro(RegularizationStrainSigma, RegularizationStrainSigmaType);
  itkGetConstReferenceMacro(RegularizationStrainSigma, RegularizationStrainSigmaType);

  /** Maximum number of iterations during regularization at the bottom level. */
  itkSetMacro(RegularizationMaximumNumberOfIterations, unsigned int);
  itkGetConstMacro(RegularizationMaximumNumberOfIterations, unsigned int);

  /** Get the multi-resolution image registration method. */
  itkGetModifiableObjectMacro(MultiResolutionRegistrationMethod, RegistrationMethodType);

protected:
  DisplacementPipeline();

  void
  GenerateOutputInformation() override;

  virtual void
  SetupPipeline();

  /** Execute the pipeline. */
  void
  GenerateData() override;


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
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkBlockMatchingDisplacementPipeline.hxx"
#endif

#endif
