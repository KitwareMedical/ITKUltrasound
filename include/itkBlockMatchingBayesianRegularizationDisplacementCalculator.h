#ifndef __itkBlockMatchingBayesianRegularizationDisplacementCalculator_h
#define __itkBlockMatchingBayesianRegularizationDisplacementCalculator_h

#include "itkConstantBoundaryCondition.h"
#include "itkConstantImagePointerBoundaryCondition.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkNeighborhoodIterator.h"

#include "itkBlockMatchingMetricImageToDisplacementCalculator.h"
#include "itkZeroFluxNeumannPadImageFilter.h"

namespace itk
{
namespace BlockMatching
{

/** \class BayesianRegularizationDisplacementCalculator
 *
 * \brief The metric image is regularized by neighboring metric images according
 * to a Bayesian model before displacments are calculated.
 *
 * A metric image is viewed as a probability density function by shifting so the
 * theoretical lower bound of the metric has the value 0 and scaled so the sum
 * image values equals unity.
 *
 * A local probability image is then iteratively modified by neighboring
 * probability images according to
 *
 * @todo insert equation.
 *
 * This was originally described in Hayton et al.  A non-rigid registration
 * algorithm for dynamic breast MR images.  Artificial Intelligence.  1999.
 *
 * Regularization is performed, then a second displacment calculator is used to
 * calculate the displacement from modified metric images.  This can be set with
 * SetDisplacementCalculator(), and it defaults to
 * MaximumPixelDisplacementCalculator.
 *
 * Assumes that all metric images have th same spacing.
 *
 * \ingroup Ultrasound
 */
template < class TMetricImage, class TDisplacementImage >
class ITK_EXPORT BayesianRegularizationDisplacementCalculator:
  public itk::BlockMatching::MetricImageToDisplacementCalculator<
    TMetricImage, TDisplacementImage >
{
public:
  /** Standard class typedefs. */
  typedef BayesianRegularizationDisplacementCalculator Self;
  typedef MetricImageToDisplacementCalculator< TMetricImage, TDisplacementImage >
  Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TDisplacementImage::ImageDimension);

  /** Run-time type information (and related methods). */
  itkTypeMacro( BayesianRegularizationDisplacementCalculator,
                MetricImageToDisplacementCalculator );

  typedef typename Superclass::MetricImageType        MetricImageType;
  typedef typename Superclass::MetricImagePointerType MetricImagePointerType;
  typedef typename MetricImageType::ConstPointer
    MetricImageConstPointerType;

  typedef typename MetricImageType::PixelType   PixelType;
  typedef typename MetricImageType::SpacingType SpacingType;
  typedef typename MetricImageType::SizeType    SizeType;
  typedef typename MetricImageType::RegionType  RegionType;
  typedef typename itk::ImageRegionIterator< MetricImageType >
    MetricImageIteratorType;
  typedef typename itk::ImageRegionConstIterator< MetricImageType >
    MetricImageConstIteratorType;

  typedef typename Superclass::IndexType IndexType;

  typedef typename Superclass::CenterPointsImageType CenterPointsImageType;

  typedef typename Superclass::MetricImageImageType  MetricImageImageType;
  typedef typename Superclass::MetricImageImagePointerType
    MetricImageImagePointerType;
  typedef typename itk::ImageRegionIterator< MetricImageImageType >
    MetricImageImageIteratorType;
  typedef NeighborhoodIterator< MetricImageImageType >
    MetricImageImageNeighborhoodIteratorType;

  typedef typename Superclass::DisplacementImageType DisplacementImageType;
  typedef typename Superclass::PointType             PointType;
  typedef typename PointType::VectorType             VectorType;

  virtual void Compute();

  /** Maximum number of iterations before regularization stops. */
  itkSetMacro( MaximumIterations, unsigned int );
  itkGetConstMacro( MaximumIterations, unsigned int );

  /** Threshold for minimum mean change in probability image when regularization
   * stops. */
  void SetMeanChangeThreshold( const  double& threshold )
    {
    m_MeanChangeThresholdDefined = true;
    m_MeanChangeThreshold = threshold;
    this->Modified();
    }
  itkGetConstMacro( MeanChangeThreshold, double );

  /** Get the mean absolute value of the probability images' delta at the
   * current iteration. This will only get updated if the MeanChangeThreshold
   * has been set.  */
  itkGetConstMacro( MeanChange, double );

  /** Theoretical lower bound of the metric.  This must be specified. */
  void SetMetricLowerBound( const double& bound )
    {
    this->m_MetricLowerBound = static_cast< PixelType >( bound );
    this->m_MetricLowerBoundDefined = true;
    this->Modified();
    }
  itkGetConstMacro( MetricLowerBound, PixelType );

  /** @todo document */
  itkSetMacro( StrainSigma, SpacingType );
  itkGetConstMacro( StrainSigma, SpacingType );

  /** @todo document Defaults ot 3*StrainSigma.  The maximum change in
   * displacement when regularizing along direction i is the MaximumStrain[i] *
   * DisplacementSpacing[i].  */
  itkSetMacro( MaximumStrain, SpacingType );
  itkGetConstMacro( MaximumStrain, SpacingType );

  /** The current iteration. */
  itkGetConstMacro( CurrentIteration, unsigned int );

  virtual void ModifyGenerateInputRequestedRegion( RegionType & region )
    {
    region = this->m_DisplacementImage->GetLargestPossibleRegion();
    }

  virtual void ModifyEnlargeOutputRequestedRegion( DataObject* data )
    {
    data->SetRequestedRegionToLargestPossibleRegion();
    }

  /** Set/Get the internal displacement calculator that is used to calculate the
   * displacements after regularization.  Defaults to a
   * MaximumPixelDisplacementCalcultor. */
  itkSetObjectMacro( DisplacementCalculator, Superclass );
  itkGetObjectMacro( DisplacementCalculator, Superclass );

protected:
  BayesianRegularizationDisplacementCalculator();

  void AllocatePriorPrImage();

  /** Scale the metric image to unity, and calculate the mean change in
   * probability if needed. */
  virtual void ScaleToUnity();

  /** Generate the gaussian like kernel for conditional probability calculation.
   * */
  void GenerateGaussianLikeKernels();

  /** Calculate apply the conditional probability of neighbors */
  void ImpartLikelihood( MetricImagePointerType& postImage,
    MetricImagePointerType& priorImage,
    const unsigned int direction,
    const VectorType& shift
  );

  typename Superclass::Pointer m_DisplacementCalculator;

  unsigned int  m_MaximumIterations;
  double        m_MeanChangeThreshold;
  bool          m_MeanChangeThresholdDefined;
  double        m_MeanChange;
  PixelType     m_MetricLowerBound;
  bool          m_MetricLowerBoundDefined;

  MetricImageImagePointerType m_PriorPr;

  typedef typename itk::ConstantImagePointerBoundaryCondition <
    MetricImageImageType > ImageImageBoundaryConditionType;
  ImageImageBoundaryConditionType m_ImageImageBoundaryCondition;

  typedef MetricImageType GaussianKernelType;
  typedef typename std::vector< typename GaussianKernelType::Pointer >
    GaussianKernelArrayType;
  GaussianKernelArrayType m_GaussianKernels;
  typedef typename GaussianKernelType::SizeType GaussianKernelRadiusType;
  typedef typename std::vector< GaussianKernelRadiusType >
    GaussianKernelRadiusArrayType;
  GaussianKernelRadiusArrayType m_GaussianKernelRadii;

  SpacingType m_StrainSigma;
  SpacingType m_MaximumStrain;

  unsigned int m_CurrentIteration;

  typedef ZeroFluxNeumannPadImageFilter< MetricImageType,
        MetricImageType > PadFilterType;

  /** Function used as a "callback" by the MultiThreader.  The threading
   * library will call this routine for each thread, which will delegate the
   * control to ThreadedGenerateData(). */
  static ITK_THREAD_RETURN_TYPE ImpartLikelihoodThreaderCallback( void *arg );

  /** Internal structure used for passing values to the threading library. */
  struct ImpartLikelihoodThreadStruct
    {
    Self * self;
    };

  typedef typename Superclass::ThreadFunctor  ThreadFunctor;
  typedef typename Superclass::ThreadStruct   ThreadStruct;

  /** We shift the minimum value of the metric image so 0 corresponds to
   * the theoretical lower bound. */
  class SubtractLowerBoundThreadFunctor : public ThreadFunctor
    {
  public:
    virtual ITK_THREAD_RETURN_TYPE operator() ( Superclass *superclass,
      RegionType& region, int threadId );
    };
  SubtractLowerBoundThreadFunctor m_SubtractLowerBoundThreadFunctor;

  /** Scale the metric images so all values sum to unity. */
  class ScaleToUnityThreadFunctor : public ThreadFunctor
    {
  public:
    virtual ITK_THREAD_RETURN_TYPE operator() ( Superclass *superclass,
      RegionType& region, int threadId );
    };
  ScaleToUnityThreadFunctor m_ScaleToUnityThreadFunctor;

  /** Copy the contents of the prior metric images to the posterior metric
   * images. */
  class CopyPriorToPosteriorThreadFunctor : public ThreadFunctor
    {
  public:
    virtual ITK_THREAD_RETURN_TYPE operator() ( Superclass *superclass,
      RegionType& region, int threadId );
    };
  CopyPriorToPosteriorThreadFunctor m_CopyPriorToPosteriorThreadFunctor;


  /** Calculate the mean change in probability. */
  class MeanChangeThreadFunctor : public ThreadFunctor
    {
  public:
    virtual ITK_THREAD_RETURN_TYPE operator() ( Superclass *superclass,
      RegionType& region, int threadId );
    };
  MeanChangeThreadFunctor m_MeanChangeThreadFunctor;
  /** Some helper quanitites for the mean change calculator. */
  double             m_ChangeSum;
  unsigned long long m_ChangeCount;

private:
  BayesianRegularizationDisplacementCalculator( const Self & );
  void operator=( const Self & );
};

} // end namespace itk
} // end namespace BlockMatching

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlockMatchingBayesianRegularizationDisplacementCalculator.txx"
#endif

#endif
