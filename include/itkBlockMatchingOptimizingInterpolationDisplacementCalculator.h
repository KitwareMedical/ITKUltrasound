/*=========================================================================
 *
 *  Copyright NumFOCUS
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         https://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkBlockMatchingOptimizingInterpolationDisplacementCalculator_h
#define itkBlockMatchingOptimizingInterpolationDisplacementCalculator_h

#include "itkBlockMatchingMetricImageToDisplacementCalculator.h"

#include "itkInterpolateImageFunction.h"

#include "itkSingleValuedNonLinearOptimizer.h"

namespace itk
{
namespace BlockMatching
{

/** \class OptimizingInterpolationDisplacementCalculator
 *
 * \brief The displacement around the maximum pixel is interpolated with an
 * interpolator and the maximizer is search for with a numerical optimization
 * algorithm.
 *
 * A highly accurate interpolator, such as a windowed sinc interpolator, is
 * recommended.
 *
 * \sa WindowedSincInterpolateImageFunction
 *
 * \ingroup Ultrasound
 */
template <typename TMetricImage, typename TDisplacementImage, typename TCoordRep = double>
class ITK_TEMPLATE_EXPORT OptimizingInterpolationDisplacementCalculator
  : public MetricImageToDisplacementCalculator<TMetricImage, TDisplacementImage>
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(OptimizingInterpolationDisplacementCalculator);

  /** Standard class type alias. */
  using Self = OptimizingInterpolationDisplacementCalculator;
  using Superclass = MetricImageToDisplacementCalculator<TMetricImage, TDisplacementImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(OptimizingInterpolationDisplacementCalculator);

  using MetricImageType = typename Superclass::MetricImageType;
  using MetricImagePointerType = typename Superclass::MetricImagePointerType;
  using PixelType = typename MetricImageType::PixelType;
  using SpacingType = typename MetricImageType::SpacingType;
  using PointType = typename Superclass::PointType;
  using IndexType = typename Superclass::IndexType;

  /** Type of the interpolator. */
  using InterpolatorType = typename itk::InterpolateImageFunction<MetricImageType, TCoordRep>;
  using ContinuousIndexType = typename InterpolatorType::ContinuousIndexType;

  /** Type of the optimizer. */
  using OptimizerType = SingleValuedNonLinearOptimizer;

  void
  SetMetricImagePixel(const PointType & point, const IndexType & index, MetricImageType * image) override;

  void
  Compute() override
  {
    // We do this here instead of SetMetricImagePixel so it only has to be done
    // once.
    this->m_DisplacementImage->Modified();
  };

  /** Set the interpolator.  Windowed sinc interpolators are recommended. */
  virtual void
  SetInterpolator(InterpolatorType * interpolator)
  {
    m_CostFunction->SetInterpolator(interpolator);
    this->Modified();
  }
  InterpolatorType *
  GetInterpolator()
  {
    return m_CostFunction->GetInterpolator();
  }

  /** Set the optimizer.  The parameter step size is in pixel indices; it should
   * be less than unity.  Since some optimizer only support minimization, the
   * metric image value is inverted. */
  virtual void
  SetOptimizer(OptimizerType * optimizer)
  {
    m_Optimizer = optimizer;
    m_Optimizer->SetCostFunction(m_CostFunction);
    this->Modified();
  }
  itkGetConstObjectMacro(Optimizer, OptimizerType);

  using CostFunctionType = typename OptimizerType::CostFunctionType;

  /** Makes the interpolator into a cost function. */
  class OptimizingInterpolationCostFunction : public CostFunctionType
  {
  public:
    using Self = OptimizingInterpolationCostFunction;
    using Superclass = CostFunctionType;

    using Pointer = SmartPointer<Self>;

    itkOverrideGetNameOfClassMacro(OptimizingInterpolationCostFunction);

    itkNewMacro(Self);

    /**  MeasureType type alias.
     *  It defines a type used to return the cost function value. */
    using MeasureType = Superclass::MeasureType;

    /**  ParametersType type alias.
     *  It defines a position in the optimization search space. */
    using ParametersType = Superclass::ParametersType;

    /** DerivativeType type alias.
     *  It defines a type used to return the cost function derivative.  */
    using DerivativeType = Superclass::DerivativeType;

    /** This method returns the value of the cost function corresponding
     * to the specified parameters.    */
    MeasureType
    GetValue(const ParametersType & parameters) const override;

    /** This method returns the derivative of the cost function corresponding
     * to the specified parameters.   */
    void
    GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const override;

    /** Return the number of parameters required to compute
     *  this cost function.
     *  This method MUST be overloaded by derived classes. */
    unsigned int
    GetNumberOfParameters() const override
    {
      return ImageDimension;
    }

    /** Set the previous parameters before the start of optimization.  Used for
     * calculating the step length in GetDerivative(). */
    void
    Initialize(const ParametersType & parameters)
    {
      m_PreviousParameters = parameters;
      this->m_PreviousParametersPtr = this->m_PreviousParameters.data_block();
    }

    itkSetObjectMacro(Interpolator, InterpolatorType);
    itkGetConstObjectMacro(Interpolator, InterpolatorType);

    typename InterpolatorType::Pointer m_Interpolator;

  protected:
    OptimizingInterpolationCostFunction();

  private:
    ParametersType              m_PreviousParameters;
    ParametersType::ValueType * m_PreviousParametersPtr;

    ContinuousIndexType                       m_ContinuousIndex;
    typename ContinuousIndexType::ValueType * m_ContinuousIndexPtr;
  };

protected:
  OptimizingInterpolationDisplacementCalculator();

  typename OptimizingInterpolationCostFunction::Pointer m_CostFunction;
  typename OptimizerType::Pointer                       m_Optimizer;

private:
};


} // namespace BlockMatching
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkBlockMatchingOptimizingInterpolationDisplacementCalculator.hxx"
#endif

#endif
