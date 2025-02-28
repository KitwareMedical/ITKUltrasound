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
template <typename TMultiResolutionMethod>
class MultiResolutionIterationDisplacementCalculatorCommand
  : public MultiResolutionIterationCommand<TMultiResolutionMethod>
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(MultiResolutionIterationDisplacementCalculatorCommand);

  using Self = MultiResolutionIterationDisplacementCalculatorCommand;
  using Superclass = MultiResolutionIterationCommand<TMultiResolutionMethod>;
  using Pointer = SmartPointer<Self>;

  itkNewMacro(Self);

  void
  Execute(itk::Object * caller, const itk::EventObject & event) override
  {
    Execute((const itk::Object *)caller, event);
  }

  void
  Execute(const itk::Object * object, const itk::EventObject & event) override;

  using MultiResolutionMethod = TMultiResolutionMethod;
  using MetricImageType = typename MultiResolutionMethod::MetricImageType;
  using DisplacementImageType = typename MultiResolutionMethod::DisplacementImageType;

  using ImageRegistrationMethodType = typename MultiResolutionMethod::ImageRegistrationMethodType;
  typedef typename ImageRegistrationMethodType::MetricImageToDisplacementCalculatorType
    MetricImageToDisplacementCalculatorType;
  using MetricImageToDisplacementCalculatorPointer = typename MetricImageToDisplacementCalculatorType::Pointer;

  using RegularizerType = BayesianRegularizationDisplacementCalculator<MetricImageType, DisplacementImageType>;
  using RegularizerPointer = typename RegularizerType::Pointer;

  itkSetObjectMacro(Level0ToNMinus1DisplacementCalculator, MetricImageToDisplacementCalculatorType);
  itkGetConstObjectMacro(Level0ToNMinus1DisplacementCalculator, MetricImageToDisplacementCalculatorType);

  itkSetObjectMacro(LevelNDisplacementCalculator, MetricImageToDisplacementCalculatorType);
  itkGetConstObjectMacro(LevelNDisplacementCalculator, MetricImageToDisplacementCalculatorType);

  /** This is only used if it is set. */
  itkSetObjectMacro(Regularizer, RegularizerType);
  itkGetConstObjectMacro(Regularizer, RegularizerType);

  itkSetMacro(Level0ToNMinus1RegularizerIterations, unsigned int);
  itkGetConstMacro(Level0ToNMinus1RegularizerIterations, unsigned int);

  itkSetMacro(LevelNRegularizerIterations, unsigned int);
  itkGetConstMacro(LevelNRegularizerIterations, unsigned int);

protected:
  MultiResolutionIterationDisplacementCalculatorCommand()
    : m_Level0ToNMinus1DisplacementCalculator(nullptr)
    , m_LevelNDisplacementCalculator(nullptr)
    , m_Regularizer(nullptr)
    , m_Level0ToNMinus1RegularizerIterations(2)
    , m_LevelNRegularizerIterations(1)
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
