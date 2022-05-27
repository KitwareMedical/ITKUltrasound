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
#ifndef itkBlockMatchingMultiResolutionIterationDisplacementCalculatorCommand_hxx
#define itkBlockMatchingMultiResolutionIterationDisplacementCalculatorCommand_hxx


namespace itk
{
namespace BlockMatching
{

template <typename TMultiResolutionMethod>
void
MultiResolutionIterationDisplacementCalculatorCommand<TMultiResolutionMethod>::Execute(const itk::Object *      object,
                                                                                       const itk::EventObject & event)
{
  Superclass::Execute(object, event);

  const unsigned long level = this->m_MultiResolutionMethod->GetCurrentLevel();

  if (m_Level0ToNMinus1DisplacementCalculator.GetPointer() == nullptr)
  {
    itkExceptionMacro(<< "Level0ToNMinus1ToNMinus1DisplacementCalculator is not present.");
  }

  if (m_LevelNDisplacementCalculator.GetPointer() == nullptr)
  {
    itkExceptionMacro(<< "LevelNDisplacementCalculator is not present.");
  }

  if (m_Regularizer.GetPointer() == nullptr)
  {
    if (level == this->m_MultiResolutionMethod->GetNumberOfLevels() - 1)
    {
      this->m_MultiResolutionMethod->GetModifiableImageRegistrationMethod()->SetMetricImageToDisplacementCalculator(
        this->m_LevelNDisplacementCalculator);
    }
    else
    {
      this->m_MultiResolutionMethod->GetImageRegistrationMethod()->SetMetricImageToDisplacementCalculator(
        this->m_Level0ToNMinus1DisplacementCalculator);
    }
  }
  else
  {
    if (level == this->m_MultiResolutionMethod->GetNumberOfLevels() - 1)
    {
      m_Regularizer->SetDisplacementCalculator(this->m_LevelNDisplacementCalculator);
      m_Regularizer->SetMaximumIterations(m_LevelNRegularizerIterations);
    }
    else
    {
      m_Regularizer->SetDisplacementCalculator(this->m_Level0ToNMinus1DisplacementCalculator);
      m_Regularizer->SetMaximumIterations(m_Level0ToNMinus1RegularizerIterations);
    }
    // this->m_MultiResolutionMethod->GetImageRegistrationMethod()->SetMetricImageToDisplacementCalculator(
    // this->m_Regularizer );
  }
}

} // end namespace BlockMatching
} // end namespace itk

#endif
