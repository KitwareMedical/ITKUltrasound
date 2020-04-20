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
#ifndef itkBlockMatchingMultiResolutionSimilarityFunctionSearchRegionImageSource_h
#define itkBlockMatchingMultiResolutionSimilarityFunctionSearchRegionImageSource_h

#include "itkArray.h"
#include "itkInterpolateImageFunction.h"
#include "itkVector.h"

#include "itkBlockMatchingMultiResolutionSearchRegionImageSource.h"
#include "itkBlockMatchingMetricImageToDisplacementCalculator.h"

#include "itkResampleIdentityNeumannImageFilter.h"

namespace itk
{
namespace BlockMatching
{

/** \class MultiResolutionSimilarityFunctionSearchRegionImageSource
 *
 * \brief The search region at subsequent levels is a function of the similarity
 * metric at the displacement.
 *
 * This is a SearchRegionImageSource that defines the search region size
 * dynamically according the metric image performance at higher levels.  The
 * search region given the metric images at the current level is determined by
 * the value of the similarity at the determined displacement.
 * This size is interpolated is interpolated at the following levels.
 *
 * In order to examine the metric image values, this class also inherits from
 * MetricImageToDisplacementCalculator.
 *
 * A second displacment calculator is used to
 * calculate the displacement from similarity metric images.  This can be set with
 * SetDisplacementCalculator(), and it defaults to
 * MaximumPixelDisplacementCalculator.
 *
 * The specified functor is defined with SetFunctor() and it should evaluate to
 * 1.0 to keep the search region the same size at the next level or less than
 * one to decrease the size at the next level.  Candidate functions are
 * available in the itk::Function:: namespace.
 *
 * @todo more with this.  top level / 2 ?
 * The search region radius will be truncated to the top level radius after
 * evaluation with the functor.
 *
 * \ingroup Ultrasound
 * */
template <class TFixedImage,
          class TMovingImage,
          class TMetricImage,
          class TDisplacementImage,
          class TFunctor,
          class TInterpolatorPrecisionType = double>
class ITK_TEMPLATE_EXPORT MultiResolutionSimilarityFunctionSearchRegionImageSource
  : public MultiResolutionSearchRegionImageSource<TFixedImage, TMovingImage, TDisplacementImage>
  , public MetricImageToDisplacementCalculator<TMetricImage, TDisplacementImage>
{
public:
  /** Standard class type alias. */
  using Self = MultiResolutionSimilarityFunctionSearchRegionImageSource;
  using Superclass = MultiResolutionSearchRegionImageSource<TFixedImage, TMovingImage, TDisplacementImage>;
  using DisplacementCalculatorSuperclass = MetricImageToDisplacementCalculator<TMetricImage, TDisplacementImage>;

  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int, TMovingImage::ImageDimension);

  /** Type of the fixed image. */
  using FixedImageType = TFixedImage;
  using FixedRegionType = typename FixedImageType::RegionType;

  /** Type of the radius used to characterized the fixed image block. */
  using RadiusType = typename FixedImageType::SizeType;

  /** Type of the moving image. */
  using MovingImageType = TMovingImage;
  using MovingRegionType = typename MovingImageType::RegionType;

  /** Type of the metric image. */
  using MetricImageType = TMetricImage;

  /** Type of the image of metric images. */
  using MetricImageImageType = typename DisplacementCalculatorSuperclass::MetricImageImageType;
  using CenterPointsImageType = typename DisplacementCalculatorSuperclass::CenterPointsImageType;

  /** Interpolator type alias. */
  using InterpolatorType = InterpolateImageFunction<MetricImageType, TInterpolatorPrecisionType>;
  using InterpolatorPointerType = typename InterpolatorType::Pointer;

  /** Type of the search region image. */
  using OutputImageType = typename Superclass::OutputImageType;
  using OutputRegionType = typename OutputImageType::RegionType;
  using SearchRegionRadiusImageType =
    Image<typename itk::Vector<typename RadiusType::SizeValueType, ImageDimension>, ImageDimension>;
  using SearchRegionRadiusImagePointer = typename SearchRegionRadiusImageType::Pointer;
  using SearchRegionRadiusResamplerType =
    VectorResampleIdentityNeumannImageFilter<SearchRegionRadiusImageType, SearchRegionRadiusImageType>;

  /** ScheduleType type alias support. */
  using PyramidScheduleType = typename Superclass::PyramidScheduleType;

  /** OverlapScheduleType type alias support. */
  using OverlapScheduleType = typename Superclass::OverlapScheduleType;

  /** Type of the displacement image from the previous level. */
  using DisplacementImageType = TDisplacementImage;

  /** Type of the function object used to calculate the search region size. */
  using FunctorType = TFunctor;

  /** Type of the filter used to resample the deformations. */
  using DisplacementResamplerType = typename Superclass::DisplacementResamplerType;
  using DisplacementResamplerPointer = typename DisplacementResamplerType::Pointer;

  /** Types inherited from the DisplacementCalculator superclass. */
  using PointType = typename DisplacementCalculatorSuperclass::PointType;
  using IndexType = typename DisplacementCalculatorSuperclass::IndexType;

  /** Type of the minimum search region radius factor. */
  using RadiusFactorType = FixedArray<double, ImageDimension>;

  /** Run-time type information (and related methods). */
  itkTypeMacro(MultiResolutionSimilarityFunctionSearchRegionImageSource, MultiResolutionFixedSearchRegionImageSource);

  /** New macro for creation of through a Smart Pointer is not used because of
   * ambiguities with LightObject. */
  static Pointer
  New(void)
  {
    Pointer smartPtr = ObjectFactory<Self>::Create();

    if (smartPtr.GetPointer() == nullptr)
    {
      smartPtr = new Self;
    }
    smartPtr->UnRegister();
    return smartPtr;
  }

  // virtual LightObject::Pointer CreateAnother(void) const
  // {
  // LightObject::Pointer smartPtr;
  // // @todo fix me -- itk::LightObject ambiguity somewhere in the following
  // line
  // smartPtr = Self::New().GetPointer();
  // return smartPtr;
  // }

  /** Set/Get the search region radius at the top, highest level ( level 0 ). */
  virtual void
  SetTopLevelRadius(const RadiusType & radius)
  {
    m_TopLevelRadius = radius;
    m_TopLevelRadiusSpecified = true;
    this->Modified();
  }

  /** Set/Get the minimum search region radius factor.  This is multiplied by
   * the fixed block size to obtain the minimum search region radius. */
  virtual void
  SetMinimumSearchRegionRadiusFactor(const RadiusFactorType & factor)
  {
    m_MinimumSearchRegionRadiusFactor = factor;
    this->Modified();
  }
  virtual void
  SetMinimumSearchRegionRadiusFactor(const double & factor)
  {
    RadiusFactorType rf(factor);
    this->SetMinimumSearchRegionRadiusFactor(rf);
  }
  itkGetConstReferenceMacro(MinimumSearchRegionRadiusFactor, RadiusFactorType);

  /** We allocate the previous search region image. */
  virtual void
  SetDisplacementImage(DisplacementImageType * image);

  /** Calculates the search region radius based on the metric image. */
  virtual void
  SetMetricImagePixel(const PointType & point, const IndexType & index, MetricImageType * image);

  virtual void
  Compute();

  /** This is needed to avoid resolution ambiguities that occur with multiple
   * inheritance. */
  virtual void
  Register() const
  {
    Superclass::Register();
  }

  virtual void
  UnRegister() const
  {
    Superclass::UnRegister();
  }

  virtual void
  Modified() const
  {
    Superclass::Modified();
  }

  bool
  GetDebug() const
  {
    return Superclass::GetDebug();
  }

  /** Set/Get the internal displacement calculator that is used to calculate the
   * displacements after regularization.  Defaults to a
   * MaximumPixelDisplacementCalcultor. */
  itkSetObjectMacro(DisplacementCalculator, DisplacementCalculatorSuperclass);
  itkGetConstObjectMacro(DisplacementCalculator, DisplacementCalculatorSuperclass);

  /** Set the interpolator function.  The default is
   * itk::LinearInterpolateImageFunction<InputImageType, TInterpolatorPrecisionType>. Some
   * other options are itk::NearestNeighborInterpolateImageFunction
   * (useful for binary masks and other images with a small number of
   * possible pixel values), and itk::BSplineInterpolateImageFunction
   * (which provides a higher order of interpolation).  */
  itkSetObjectMacro(Interpolator, InterpolatorType);

  /** Get a pointer to the interpolator function. */
  itkGetConstObjectMacro(Interpolator, InterpolatorType);

  /** Get the functor object.  The functor is returned by reference.
   * (Functors do not have to derive from itk::LightObject, so they do
   * not necessarily have a reference count. So we cannot return a
   * SmartPointer.) */
  const FunctorType &
  GetFunctor() const
  {
    return m_Functor;
  }

  /** Set the functor object.  This replaces the current Functor with a
   * copy of the specified Functor. This allows the user to specify a
   * functor that has ivars set differently than the default functor.
   * This method requires an operator!=() be defined on the functor
   * (or the compiler's default implementation of operator!=() being
   * appropriate). */
  void
  SetFunctor(const FunctorType & functor)
  {
    if (m_Functor != functor)
    {
      m_Functor = functor;
      this->Modified();
    }
  }

protected:
  MultiResolutionSimilarityFunctionSearchRegionImageSource();

  virtual void
  BeforeThreadedGenerateData() override;

  virtual void
  ThreadedGenerateData(const OutputRegionType & outputRegion, ThreadIdType threadID);

  RadiusType m_TopLevelRadius;
  bool       m_TopLevelRadiusSpecified;

  RadiusFactorType m_MinimumSearchRegionRadiusFactor;

  typename SearchRegionRadiusImageType::Pointer     m_SearchRegionRadiusImage;
  typename SearchRegionRadiusResamplerType::Pointer m_SearchRegionRadiusResampler;

  typename DisplacementCalculatorSuperclass::Pointer m_DisplacementCalculator;

  InterpolatorPointerType m_Interpolator; // Image function for
                                          // metric interpolation
  FunctorType m_Functor;

private:
  MultiResolutionSimilarityFunctionSearchRegionImageSource(const Self &);
  void
  operator=(const Self &);
};

} // end namespace BlockMatching
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkBlockMatchingMultiResolutionSimilarityFunctionSearchRegionImageSource.hxx"
#endif

#endif
