#ifndef __itkBlockMatchingMultiResolutionIterationDisplacementCalculatorCommand_h
#define __itkBlockMatchingMultiResolutionIterationDisplacementCalculatorCommand_h

#include "itkBlockMatchingMultiResolutionIterationCommand.h"

#include "itkBlockMatchingBayesianRegularizationDisplacementCalculator.h"

namespace itk
{
namespace BlockMatching
{

/** \class MultiResolutionIterationDisplacementCalculatorCommand
 * \brief Change the displacement calculator and regulator iterations depending on the level in the pyramid. */
template< class TMultiResolutionMethod >
class MultiResolutionIterationDisplacementCalculatorCommand :
  public MultiResolutionIterationCommand< TMultiResolutionMethod >
{
public:
  typedef MultiResolutionIterationDisplacementCalculatorCommand   Self;
  typedef MultiResolutionIterationCommand<TMultiResolutionMethod> Superclass;
  typedef SmartPointer<Self>                                      Pointer;

  itkNewMacro( Self );

  virtual void Execute(itk::Object *caller, const itk::EventObject & event)
    {
    Execute( (const itk::Object *)caller, event);
    }

  virtual void Execute(const itk::Object * object, const itk::EventObject & event);

  typedef TMultiResolutionMethod                                MultiResolutionMethod;
  typedef typename MultiResolutionMethod::MetricImageType       MetricImageType;
  typedef typename MultiResolutionMethod::DisplacementImageType DisplacementImageType;

  typedef typename MultiResolutionMethod::ImageRegistrationMethodType ImageRegistrationMethodType;
  typedef typename ImageRegistrationMethodType::MetricImageToDisplacementCalculatorType
    MetricImageToDisplacementCalculatorType;
  typedef typename MetricImageToDisplacementCalculatorType::Pointer MetricImageToDisplacementCalculatorPointer;

  typedef BayesianRegularizationDisplacementCalculator< MetricImageType, DisplacementImageType >
                                            RegularizerType;
  typedef typename RegularizerType::Pointer RegularizerPointer;

  itkSetObjectMacro( Level0ToNMinus1DisplacementCalculator, MetricImageToDisplacementCalculatorType );
  itkGetObjectMacro( Level0ToNMinus1DisplacementCalculator, MetricImageToDisplacementCalculatorType );

  itkSetObjectMacro( LevelNDisplacementCalculator, MetricImageToDisplacementCalculatorType );
  itkGetObjectMacro( LevelNDisplacementCalculator, MetricImageToDisplacementCalculatorType );

  /** This is only used if it is set. */
  itkSetObjectMacro( Regularizer, RegularizerType );
  itkGetObjectMacro( Regularizer, RegularizerType );

  itkSetMacro( Level0ToNMinus1RegularizerIterations, unsigned int );
  itkGetConstMacro( Level0ToNMinus1RegularizerIterations, unsigned int );

  itkSetMacro( LevelNRegularizerIterations, unsigned int );
  itkGetConstMacro( LevelNRegularizerIterations, unsigned int );

protected:
  MultiResolutionIterationDisplacementCalculatorCommand():
    m_Level0ToNMinus1DisplacementCalculator( NULL ),
    m_LevelNDisplacementCalculator( NULL ),
    m_Regularizer( NULL ),
    m_Level0ToNMinus1RegularizerIterations( 2 ),
    m_LevelNRegularizerIterations( 1 )
  {}

  MetricImageToDisplacementCalculatorPointer m_Level0ToNMinus1DisplacementCalculator;
  MetricImageToDisplacementCalculatorPointer m_LevelNDisplacementCalculator;

  RegularizerPointer m_Regularizer;

  unsigned int m_Level0ToNMinus1RegularizerIterations;
  unsigned int m_LevelNRegularizerIterations;

private:
  MultiResolutionIterationDisplacementCalculatorCommand( const Self& );
  void operator=( const Self & );
};

} // end namespace BlockMatching
} // end namespace itk

#include "itkBlockMatchingMultiResolutionIterationDisplacementCalculatorCommand.txx"

#endif
