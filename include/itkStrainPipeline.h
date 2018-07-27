#ifndef __StrainPipeline_h
#define __StrainPipeline_h

#include "itkAmoebaOptimizer.h"
#include "itkExpNegativeImageFilter.h"
#include "itkGradientImageFilter.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageToImageFilter.h"
#include "itkVector.h"
#include "itkVTKImageIO.h"
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
#include "itkBlockMatchingMultiResolutionIterationObserver.h"
#include "itkBlockMatchingMultiResolutionSearchRegionWriterCommand.h"
#include "itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter.h"
#include "itkBlockMatchingOptimizingInterpolationDisplacementCalculator.h"
#include "itkBlockMatchingParabolicInterpolationDisplacementCalculator.h"
#include "itkBlockMatchingSearchRegionImageInitializer.h"
#include "itkBlockMatchingStrainWindowDisplacementCalculator.h"
#include "itkBlockMatchingStrainWindowBlockAffineTransformCommand.h"
#include "itkTextProgressBarCommand.h"

#include "itkBSplineApproximationGradientImageFilter.h"
#include "itkLinearLeastSquaresGradientImageFilter.h"
#include "itkHigherOrderAccurateGradientImageFilter.h"
#include "itkStrainImageFilter.h"
#include "itkSplitComponentsImageFilter.h"

#include "MMMStrainOptions.h"

namespace itk
{

/** \class StrainPipeline
 *
 * \brief Sets up and runs deformable image registration pipeline with block-matching.
 * Optionally also saves output images.
 *
 */
template< class TFixedPixel = signed short, class TMovingPixel= TFixedPixel,
          class TMetricPixel = double, class TCoordRep = double,
          unsigned int VImageDimension = 2 >
class ITK_EXPORT StrainPipeline : public ::itk::ImageToImageFilter<
  ::itk::Image< TFixedPixel, VImageDimension >,
  ::itk::Image< ::itk::Vector< TMetricPixel, VImageDimension>, VImageDimension > >
{
public:
  /** Standard class typedefs. */
  typedef StrainPipeline                 Self;
  typedef SmartPointer< Self >           Pointer;
  typedef SmartPointer< const Self >     ConstPointer;

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      VImageDimension);

  typedef TFixedPixel FixedPixelType;
  typedef Image< TFixedPixel, ImageDimension > FixedImageType;
  typedef typename FixedImageType::SizeType RadiusType;

  typedef TMovingPixel MovingPixelType;
  typedef Image< MovingPixelType, ImageDimension > MovingImageType;

  typedef TMetricPixel MetricPixelType;
  typedef Image< MetricPixelType, ImageDimension > MetricImageType;

  typedef Vector<MetricPixelType, ImageDimension > VectorType;
  typedef Image<VectorType, ImageDimension >       DisplacementImageType;
  typedef DisplacementImageType                    OutputImageType;

  typedef ImageToImageFilter< FixedImageType, DisplacementImageType > Superclass;

  /** Run-time type information (and related methods). */
  itkTypeMacro( StrainPipeline, Superclass );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  typedef TCoordRep CoordRepType;

  typedef ResampleImageFilter< FixedImageType, FixedImageType, CoordRepType > FixedResamplerType;
  typedef ResampleImageFilter< MovingImageType, MovingImageType, CoordRepType > MovingResamplerType;

  const static unsigned int RESAMPLE_RADIUS = 4;
  typedef itk::Function::WelchWindowFunction< RESAMPLE_RADIUS >    ResampleWindowType;
  typedef itk::ZeroFluxNeumannBoundaryCondition< FixedImageType >  FixedBoundaryConditionType;
  typedef itk::WindowedSincInterpolateImageFunction< FixedImageType,
                                                    RESAMPLE_RADIUS, ResampleWindowType, FixedBoundaryConditionType,
                                                    CoordRepType> FixedResamplerInterpolatorType;
  typedef itk::ZeroFluxNeumannBoundaryCondition< MovingImageType >  MovingBoundaryConditionType;
  typedef itk::WindowedSincInterpolateImageFunction< MovingImageType,
                                                    RESAMPLE_RADIUS, ResampleWindowType, MovingBoundaryConditionType,
                                                    CoordRepType> MovingResamplerInterpolatorType;

  /** The block radius calculator. */
  typedef itk::BlockMatching::MultiResolutionMinMaxBlockRadiusCalculator< FixedImageType >
    BlockRadiusCalculatorType;

  /** The search region image source. */
  typedef itk::BlockMatching::MultiResolutionMinMaxSearchRegionImageSource< FixedImageType,
          MovingImageType, DisplacementImageType > SearchRegionImageSourceType;

  /** The registration method. */
  typedef itk::BlockMatching::ImageRegistrationMethod<FixedImageType,
                                                      MovingImageType, MetricImageType, DisplacementImageType,
                                                      CoordRepType>
    LevelRegistrationMethodType;

  /** Interpolation classes. */
  typedef itk::BlockMatching::ParabolicInterpolationDisplacementCalculator<
    MetricImageType,
    DisplacementImageType> ParabolicInterpolatorType;
  typedef itk::BlockMatching::MaximumPixelDisplacementCalculator<MetricImageType,
                                                                 DisplacementImageType>
  MaximumPixelDisplacementCalculatorType;
  typedef itk::BlockMatching::OptimizingInterpolationDisplacementCalculator<MetricImageType,
                                                                            DisplacementImageType> FinalInterpolatorType;
  const static unsigned int OPTIMIZING_INTERPOLATOR_RADIUS = 4;
  typedef itk::Function::WelchWindowFunction< OPTIMIZING_INTERPOLATOR_RADIUS >  WindowType;
  typedef itk::ZeroFluxNeumannBoundaryCondition< MetricImageType >              ResampleBoundaryConditionType;
  typedef itk::WindowedSincInterpolateImageFunction<MetricImageType,
                                                    OPTIMIZING_INTERPOLATOR_RADIUS, WindowType, ResampleBoundaryConditionType,
                                                    CoordRepType> SubsampleInterpolatorType;
  typedef itk::AmoebaOptimizer SubsampleOptimizerType;

  /** Filter out peak hopping. */
  typedef itk::BlockMatching::StrainWindowDisplacementCalculator<MetricImageType, DisplacementImageType, double>
    StrainWindowDisplacementCalculatorType;
  typedef itk::StrainImageFilter<DisplacementImageType, MetricPixelType, MetricPixelType>
  StrainWindowStrainFilterType;
  typedef itk::HigherOrderAccurateGradientImageFilter<MetricImageType, MetricPixelType, MetricPixelType>
  HigherOrderAccurateGradientFilterType;
  typedef itk::LinearLeastSquaresGradientImageFilter< MetricImageType, MetricPixelType, MetricPixelType >
    LinearLeastSquaresGradientFilterType;

  /** The image similarity metric. */
  typedef itk::BlockMatching::NormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter<
    FixedImageType, MovingImageType, MetricImageType> MetricImageFilterType;

  /** Scale the fixed block by the strain at higher levels. */
  typedef itk::BlockMatching::BlockAffineTransformMetricImageFilter< FixedImageType,
          MovingImageType, MetricImageType, MetricPixelType > BlockTransformMetricImageFilterType;
  typedef itk::BlockMatching::StrainWindowBlockAffineTransformCommand< StrainWindowDisplacementCalculatorType,
          BlockTransformMetricImageFilterType, StrainWindowStrainFilterType > BlockTransformCommandType;

  /** Perform regularization. */
  typedef itk::BlockMatching::BayesianRegularizationDisplacementCalculator<
    MetricImageType, DisplacementImageType> DisplacmentRegularizerType;

  /** Multi-resolution registration method. */
  typedef itk::BlockMatching::MultiResolutionImageRegistrationMethod< FixedImageType, MovingImageType, MetricImageType,
                                                                     DisplacementImageType,
                                                                     CoordRepType> RegistrationMethodType;
  typedef itk::BlockMatching::MultiResolutionIterationObserver<RegistrationMethodType>
  MultiResolutionObserverType;

  /** Set the displacement calculator and regularizer iterations at every level. */
  typedef itk::BlockMatching::MultiResolutionIterationDisplacementCalculatorCommand<RegistrationMethodType>
  DisplacementCalculatorCommandType;

  /** Write out the search region images at every level. */
  typedef itk::BlockMatching::MultiResolutionSearchRegionWriterCommand<RegistrationMethodType>
  SearchRegionWriterCommandType;

  /** Write displacement vector image. */
  typedef itk::ImageFileWriter< DisplacementImageType > DisplacementWriterType;

  /** Calculate strains. */
  typedef itk::StrainImageFilter<DisplacementImageType, MetricPixelType, MetricPixelType>
                                            StrainFilterType;
  typedef typename StrainFilterType::OutputImageType TensorImageType;
  typedef itk::BSplineApproximationGradientImageFilter< DisplacementImageType, MetricPixelType >
    FinalGradientFilterType;
  typedef itk::ImageFileWriter< TensorImageType > StrainWriterType;
  typedef itk::VTKImageIO StrainIOType;
  typedef itk::SplitComponentsImageFilter<TensorImageType, MetricImageType, 3>
  StrainComponentsFilterType;
  typedef itk::ImageFileWriter<MetricImageType> ComponentWriterType;
  typedef itk::SplitComponentsImageFilter<DisplacementImageType,
                                          MetricImageType, ImageDimension>
  DisplacementComponentsFilterType;

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

  /** Set the options structure extracted from the configuration file. */
  void SetMMMStrainOptions( const MMMStrainOptions * options );
  const MMMStrainOptions * GetMMMStrainOptions() const
    {
    return m_MMMStrainOptions;
    }

  /** Display a text progress bar for the level registration method's creation
   * of the metric images. */
  itkSetMacro( LevelRegistrationMethodTextProgressBar, bool );
  itkGetConstMacro( LevelRegistrationMethodTextProgressBar, bool );

  /** Write output images to file. */
  itkSetMacro( WriteOutputImagesToFile, bool );
  itkGetConstMacro( WriteOutputImagesToFile, bool );

protected:
  StrainPipeline();

  virtual void GenerateOutputInformation();

  virtual void SetupPipeline();

  /** Execute the pipeline. */
  virtual void GenerateData();

  const MMMStrainOptions * m_MMMStrainOptions;
  bool                     m_WriteOutputImagesToFile;

  typename FixedResamplerType::Pointer  m_FixedResampler;
  typename MovingResamplerType::Pointer m_MovingResampler;
  typename FixedResamplerInterpolatorType::Pointer  m_FixedResamplerInterpolator;
  typename MovingResamplerInterpolatorType::Pointer m_MovingResamplerInterpolator;

  typename BlockRadiusCalculatorType::Pointer m_BlockRadiusCalculator;

  typename SearchRegionImageSourceType::Pointer m_SearchRegionImageSource;

  typename LevelRegistrationMethodType::Pointer m_LevelRegistrationMethod;
  TextProgressBarCommand::Pointer m_RegistrationObserver;
  bool                            m_LevelRegistrationMethodTextProgressBar;

  typename ParabolicInterpolatorType::Pointer m_ParabolicInterpolator;
  typename MaximumPixelDisplacementCalculatorType::Pointer m_MaximumPixelInterpolator;
  typename FinalInterpolatorType::Pointer m_FinalInterpolator;
  typename SubsampleInterpolatorType::Pointer m_SubsampleInterpolator;
  typename SubsampleOptimizerType::Pointer m_SubsampleOptimizer;

  typename StrainWindowDisplacementCalculatorType::Pointer m_StrainWindower;
  typename StrainWindowStrainFilterType::Pointer m_StrainWindowStrainFilter;
  typename HigherOrderAccurateGradientFilterType::Pointer m_HigherOrderAccurateGradientFilter;
  typename LinearLeastSquaresGradientFilterType::Pointer m_LinearLeastSquaresGradientFilter;

  typename MetricImageFilterType::Pointer m_MetricImageFilter;

  typename BlockTransformMetricImageFilterType::Pointer m_BlockTransformMetricImageFilter;
  typename BlockTransformCommandType::Pointer m_BlockTransformCommand;

  typename DisplacmentRegularizerType::Pointer m_Regularizer;

  typename RegistrationMethodType::Pointer m_MultiResRegistrationMethod;
  typename MultiResolutionObserverType::Pointer m_MultiResObserver;
  typename DisplacementCalculatorCommandType::Pointer m_DisplacementCalculatorCommand;

  typename SearchRegionWriterCommandType::Pointer m_SearchRegionWriterCommand;
  typename DisplacementWriterType::Pointer m_DisplacementWriter;
  typename StrainFilterType::Pointer m_StrainFilter;
  typename FinalGradientFilterType::Pointer m_FinalGradientFilter;
  typename StrainWriterType::Pointer m_StrainWriter;
  typename StrainIOType::Pointer m_StrainIO;
  typename StrainComponentsFilterType::Pointer m_StrainComponentsFilter;
  typename ComponentWriterType::Pointer m_StrainComponentWriter;
  typename DisplacementComponentsFilterType::Pointer m_DisplacementComponentsFilter;
  typename ComponentWriterType::Pointer m_DisplacementComponentWriter;

private:
  StrainPipeline( const Self & );
  void operator=( const Self & );
};

} // end itk namespace

#endif
