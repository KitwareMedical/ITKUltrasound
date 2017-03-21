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
#ifndef __itkBlockMatchingOptimizingInterpolationDisplacementCalculator_h
#define __itkBlockMatchingOptimizingInterpolationDisplacementCalculator_h

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
template < class TMetricImage, class TDisplacementImage, class TCoordRep=double > 
  class ITK_EXPORT OptimizingInterpolationDisplacementCalculator: public
      itk::BlockMatching::MetricImageToDisplacementCalculator<
      TMetricImage, TDisplacementImage > 
{ 
public:
  /** Standard class typedefs. */
  typedef OptimizingInterpolationDisplacementCalculator     Self;
  typedef MetricImageToDisplacementCalculator< TMetricImage, TDisplacementImage >
    Superclass;

  typedef SmartPointer< Self >           Pointer;
  typedef SmartPointer< const Self >     ConstPointer;

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      Superclass::ImageDimension);

  /** Method for creation through the object factory. */
  itkNewMacro( Self );
  
  /** Run-time type information (and related methods). */
  itkTypeMacro( OptimizingInterpolationDisplacementCalculator, MetricImageToDisplacementCalculator );

  typedef typename Superclass::MetricImageType MetricImageType;
  typedef typename Superclass::MetricImagePointerType MetricImagePointerType;
  typedef typename MetricImageType::PixelType   PixelType;
  typedef typename MetricImageType::SpacingType SpacingType;

  typedef typename Superclass::PointType PointType;
  
  typedef typename Superclass::IndexType IndexType;

  /** Type of the interpolator. */
  typedef typename itk::InterpolateImageFunction< MetricImageType, TCoordRep >
    InterpolatorType;
  typedef typename InterpolatorType::ContinuousIndexType ContinuousIndexType;

  /** Type of the optimizer. */
  typedef typename itk::SingleValuedNonLinearOptimizer OptimizerType;

  virtual void SetMetricImagePixel( const PointType & point, const IndexType& index, MetricImageType* image ); 

  virtual void Compute() {
    // We do this here instead of SetMetricImagePixel so it only has to be done
    // once.
    this->m_DisplacementImage->Modified();
  };

  /** Set the interpolator.  Windowed sinc interpolators are recommended. */
  virtual void SetInterpolator( InterpolatorType* interpolator )
    {
    m_CostFunction->SetInterpolator( interpolator );
    this->Modified();
    }
  InterpolatorType* GetInterpolator()
    {
    return m_CostFunction->GetInterpolator();
    }

  /** Set the optimizer.  The parameter step size is in pixel indices; it should
   * be less than unity.  Since some optimizer only support minimization, the
   * metric image value is inverted. */
  virtual void SetOptimizer( OptimizerType * optimizer )
    {
    m_Optimizer = optimizer;
    m_Optimizer->SetCostFunction( m_CostFunction );
    this->Modified();
    }
  itkGetConstObjectMacro( Optimizer, OptimizerType );

  typedef typename OptimizerType::CostFunctionType CostFunctionType;

   /** Makes the interpolator into a cost function. */
  class OptimizingInterpolationCostFunction : public CostFunctionType
  {
  public:
    typedef OptimizingInterpolationCostFunction Self;
    typedef CostFunctionType                    Superclass;

    typedef SmartPointer< Self >                Pointer;

    itkTypeMacro( OptimizingInterpolationCostFunction, SingleValuedCostFunction );

    itkNewMacro( Self );

    /**  MeasureType typedef.
     *  It defines a type used to return the cost function value. */
    typedef Superclass::MeasureType    MeasureType;

    /**  ParametersType typedef.
     *  It defines a position in the optimization search space. */
    typedef Superclass::ParametersType ParametersType;

    /** DerivativeType typedef.
     *  It defines a type used to return the cost function derivative.  */
    typedef Superclass::DerivativeType DerivativeType;

    /** This method returns the value of the cost function corresponding
      * to the specified parameters.    */ 
    virtual MeasureType GetValue( const ParametersType & parameters ) const;

    /** This method returns the derivative of the cost function corresponding
      * to the specified parameters.   */ 
    virtual void GetDerivative( const ParametersType & parameters,
                                DerivativeType & derivative ) const;

    /** Return the number of parameters required to compute 
     *  this cost function.
     *  This method MUST be overloaded by derived classes. */
    virtual unsigned int GetNumberOfParameters(void) const
      { return ImageDimension; }

    /** Set the previous parameters before the start of optimization.  Used for
     * calculating the step length in GetDerivative(). */
    void Initialize( const ParametersType& parameters )
      {
      m_PreviousParameters = parameters;
      this->m_PreviousParametersPtr = this->m_PreviousParameters.data_block();
      }

    itkSetObjectMacro( Interpolator, InterpolatorType );
    itkGetObjectMacro( Interpolator, InterpolatorType );

    typename InterpolatorType::Pointer m_Interpolator;

  protected:
    OptimizingInterpolationCostFunction();

  private:
    ParametersType m_PreviousParameters;
    ParametersType::ValueType* m_PreviousParametersPtr;

    ContinuousIndexType m_ContinuousIndex;
    typename ContinuousIndexType::ValueType* m_ContinuousIndexPtr;
  };

protected:
  OptimizingInterpolationDisplacementCalculator();

  typename OptimizingInterpolationCostFunction::Pointer m_CostFunction;
  typename OptimizerType::Pointer m_Optimizer;

private:
  OptimizingInterpolationDisplacementCalculator( const Self & );
  void operator=( const Self & );
};


} // end namespace itk
} // end namespace BlockMatching

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlockMatchingOptimizingInterpolationDisplacementCalculator.hxx"
#endif

#endif
