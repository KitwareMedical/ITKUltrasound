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
#ifndef itkBlockMatchingBayesianRegularizationDisplacementCalculator_h
#define itkBlockMatchingBayesianRegularizationDisplacementCalculator_h

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
template <typename TMetricImage, typename TDisplacementImage>
class ITK_TEMPLATE_EXPORT BayesianRegularizationDisplacementCalculator
  : public MetricImageToDisplacementCalculator<TMetricImage, TDisplacementImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(BayesianRegularizationDisplacementCalculator);

  /** Standard class type alias. */
  using Self = BayesianRegularizationDisplacementCalculator;
  using Superclass = MetricImageToDisplacementCalculator<TMetricImage, TDisplacementImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int, TDisplacementImage::ImageDimension);

  /** Run-time type information (and related methods). */
  itkTypeMacro(BayesianRegularizationDisplacementCalculator, MetricImageToDisplacementCalculator);

  using MetricImageType = typename Superclass::MetricImageType;
  using MetricImagePointerType = typename Superclass::MetricImagePointerType;
  typedef typename MetricImageType::ConstPointer MetricImageConstPointerType;

  using PixelType = typename MetricImageType::PixelType;
  using SpacingType = typename MetricImageType::SpacingType;
  using SizeType = typename MetricImageType::SizeType;
  using RegionType = typename MetricImageType::RegionType;
  using MetricImageIteratorType = typename itk::ImageRegionIterator<MetricImageType>;
  using MetricImageConstIteratorType = typename itk::ImageRegionConstIterator<MetricImageType>;

  using IndexType = typename Superclass::IndexType;

  using CenterPointsImageType = typename Superclass::CenterPointsImageType;

  using MetricImageImageType = typename Superclass::MetricImageImageType;
  typedef typename Superclass::MetricImageImagePointerType MetricImageImagePointerType;
  using MetricImageImageIteratorType = typename itk::ImageRegionIterator<MetricImageImageType>;
  using MetricImageImageNeighborhoodIteratorType = NeighborhoodIterator<MetricImageImageType>;

  using DisplacementImageType = typename Superclass::DisplacementImageType;
  using PointType = typename Superclass::PointType;
  using VectorType = typename PointType::VectorType;

  void
  Compute() override;

  /** Maximum number of iterations before regularization stops. */
  itkSetMacro(MaximumIterations, unsigned int);
  itkGetConstMacro(MaximumIterations, unsigned int);

  /** Threshold for minimum mean change in probability image when regularization
   * stops. */
  void
  SetMeanChangeThreshold(const double & threshold)
  {
    m_MeanChangeThresholdDefined = true;
    m_MeanChangeThreshold = threshold;
    this->Modified();
  }
  itkGetConstMacro(MeanChangeThreshold, double);

  /** Get the mean absolute value of the probability images' delta at the
   * current iteration. This will only get updated if the MeanChangeThreshold
   * has been set.  */
  itkGetConstMacro(MeanChange, double);

  /** Theoretical lower bound of the metric.  This must be specified. */
  void
  SetMetricLowerBound(const double & bound)
  {
    this->m_MetricLowerBound = static_cast<PixelType>(bound);
    this->m_MetricLowerBoundDefined = true;
    this->Modified();
  }
  itkGetConstMacro(MetricLowerBound, PixelType);

  /** @todo document */
  itkSetMacro(StrainSigma, SpacingType);
  itkGetConstMacro(StrainSigma, SpacingType);

  /** @todo document Defaults ot 3*StrainSigma.  The maximum change in
   * displacement when regularizing along direction i is the MaximumStrain[i] *
   * DisplacementSpacing[i].  */
  itkSetMacro(MaximumStrain, SpacingType);
  itkGetConstMacro(MaximumStrain, SpacingType);

  /** The current iteration. */
  itkGetConstMacro(CurrentIteration, unsigned int);

  void
  ModifyGenerateInputRequestedRegion(RegionType & region) override
  {
    region = this->m_DisplacementImage->GetLargestPossibleRegion();
  }

  void
  ModifyEnlargeOutputRequestedRegion(DataObject * data) override
  {
    data->SetRequestedRegionToLargestPossibleRegion();
  }

  /** Set/Get the internal displacement calculator that is used to calculate the
   * displacements after regularization.  Defaults to a
   * MaximumPixelDisplacementCalcultor. */
  itkSetObjectMacro(DisplacementCalculator, Superclass);
  itkGetConstObjectMacro(DisplacementCalculator, Superclass);

protected:
  BayesianRegularizationDisplacementCalculator();

  void
  AllocatePriorPrImage();

  /** Scale the metric image to unity, and calculate the mean change in
   * probability if needed. */
  virtual void
  ScaleToUnity();

  /** Generate the gaussian like kernel for conditional probability calculation.
   * */
  void
  GenerateGaussianLikeKernels();

  /** Calculate apply the conditional probability of neighbors */
  void
  ImpartLikelihood(MetricImagePointerType & postImage,
                   MetricImagePointerType & priorImage,
                   const unsigned int       direction,
                   const VectorType &       shift);
  void
  ThreadedImpartLikelihood(const RegionType & region);

  typename Superclass::Pointer m_DisplacementCalculator;

  unsigned int m_MaximumIterations;
  double       m_MeanChangeThreshold;
  bool         m_MeanChangeThresholdDefined;
  double       m_MeanChange;
  PixelType    m_MetricLowerBound;
  bool         m_MetricLowerBoundDefined;

  MetricImageImagePointerType m_PriorPr;

  using ImageImageBoundaryConditionType = typename itk::ConstantImagePointerBoundaryCondition<MetricImageImageType>;
  ImageImageBoundaryConditionType m_ImageImageBoundaryCondition;

  using GaussianKernelType = MetricImageType;
  using GaussianKernelArrayType = typename std::vector<typename GaussianKernelType::Pointer>;
  GaussianKernelArrayType m_GaussianKernels;
  using GaussianKernelRadiusType = typename GaussianKernelType::SizeType;
  using GaussianKernelRadiusArrayType = typename std::vector<GaussianKernelRadiusType>;
  GaussianKernelRadiusArrayType m_GaussianKernelRadii;

  SpacingType m_StrainSigma;
  SpacingType m_MaximumStrain;

  unsigned int m_CurrentIteration;

  using PadFilterType = ZeroFluxNeumannPadImageFilter<MetricImageType, MetricImageType>;

  /** We shift the minimum value of the metric image so 0 corresponds to
   * the theoretical lower bound. */
  void
  ThreadedSubtractLowerBound(const RegionType & region);

  /** Scale the metric images so all values sum to unity. */
  void
  ThreadedScaleToUnity(const RegionType & region);

  /** Copy the contents of the prior metric images to the posterior metric
   * images. */
  void
  ThreadedCopyPriorToPosterior(const RegionType & region);

  /** Calculate the mean change in probability. */
  void
  ThreadedMeanChange(const RegionType & region);

private:
  /** Some helper quanitites for the mean change calculator. */
  std::atomic<double>   m_ChangeSum;
  std::atomic<uint64_t> m_ChangeCount;
};

} // namespace BlockMatching
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkBlockMatchingBayesianRegularizationDisplacementCalculator.hxx"
#endif

#endif
