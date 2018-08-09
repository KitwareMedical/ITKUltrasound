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
#ifndef itkBlockMatchingMultiResolutionIterationDisplacementCalculatorCommand_h
#define itkBlockMatchingMultiResolutionIterationDisplacementCalculatorCommand_h

#include "itkBlockMatchingMultiResolutionIterationCommand.h"

#include "itkBlockMatchingBayesianRegularizationDisplacementCalculator.h"

namespace itk
{
namespace BlockMatching
{

/** \class MultiResolutionIterationDisplacementCalculatorCommand
 *
 * \brief Change the displacement calculator and regulator iterations depending on the level in the pyramid.
 *
 * \ingroup Ultrasound
 * */
template< typename TMultiResolutionMethod >
class MultiResolutionIterationDisplacementCalculatorCommand :
  public MultiResolutionIterationCommand< TMultiResolutionMethod >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(MultiResolutionIterationDisplacementCalculatorCommand);

  typedef MultiResolutionIterationDisplacementCalculatorCommand   Self;
  typedef MultiResolutionIterationCommand<TMultiResolutionMethod> Superclass;
  typedef SmartPointer<Self>                                      Pointer;

  itkNewMacro( Self );

  void Execute(itk::Object *caller, const itk::EventObject & event) override
    {
    Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event) override;

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
    m_Level0ToNMinus1DisplacementCalculator( nullptr ),
    m_LevelNDisplacementCalculator( nullptr ),
    m_Regularizer( nullptr ),
    m_Level0ToNMinus1RegularizerIterations( 2 ),
    m_LevelNRegularizerIterations( 1 )
  {}

  MetricImageToDisplacementCalculatorPointer m_Level0ToNMinus1DisplacementCalculator;
  MetricImageToDisplacementCalculatorPointer m_LevelNDisplacementCalculator;

  RegularizerPointer m_Regularizer;

  unsigned int m_Level0ToNMinus1RegularizerIterations;
  unsigned int m_LevelNRegularizerIterations;

private:
};

} // end namespace BlockMatching
} // end namespace itk

#include "itkBlockMatchingMultiResolutionIterationDisplacementCalculatorCommand.hxx"

#endif
