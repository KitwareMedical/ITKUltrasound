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
#ifndef itkBlockMatchingMultiResolutionIterationCommand_h
#define itkBlockMatchingMultiResolutionIterationCommand_h

#include "itkCommand.h"

#include "itkBlockMatchingMultiResolutionImageRegistrationMethod.h"

namespace itk
{
namespace BlockMatching
{

/**
 * \class MultiResolutionIterationCommand
 *
 * \brief This is a base class for classes that want to observe/adjust a
 * BlockMatching::MultiResolutionImageRegistrationMethod at every iteration.
 *
 * \ingroup Ultrasound
 * */
template <typename TMultiResolutionMethod>
class MultiResolutionIterationCommand : public Command
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(MultiResolutionIterationCommand);

  using Self = MultiResolutionIterationCommand;
  using Superclass = Command;
  using Pointer = SmartPointer<Self>;

  void
  Execute(itk::Object * caller, const itk::EventObject & event) override
  {
    Execute((const itk::Object *)caller, event);
  }

  void
  Execute(const itk::Object * object, const itk::EventObject & event) override;

  using MultiResolutionMethodType = TMultiResolutionMethod;

  void
  SetMultiResolutionMethod(MultiResolutionMethodType * method)
  {
    m_MultiResolutionMethod = method;
  }

  using FixedImagePyramidType = typename MultiResolutionMethodType::FixedImagePyramidType;
  using FixedImagePyramidPointer = typename FixedImagePyramidType::ConstPointer;
  using MovingImagePyramidType = typename MultiResolutionMethodType::MovingImagePyramidType;
  using MovingImagePyramidPointer = typename MovingImagePyramidType::ConstPointer;

  using BlockRadiusCalculatorType = typename MultiResolutionMethodType::BlockRadiusCalculatorType;
  using BlockRadiusCalculatorPointer = typename BlockRadiusCalculatorType::ConstPointer;

  using SearchRegionImageSourceType = typename MultiResolutionMethodType::SearchRegionImageSourceType;
  using SearchRegionImageSourcePointer = typename SearchRegionImageSourceType::ConstPointer;

protected:
  MultiResolutionIterationCommand() { m_MultiResolutionMethod = nullptr; }
  virtual ~MultiResolutionIterationCommand(){};

  typename MultiResolutionMethodType::Pointer m_MultiResolutionMethod;

  FixedImagePyramidPointer  m_FixedImagePyramid;
  MovingImagePyramidPointer m_MovingImagePyramid;

  BlockRadiusCalculatorPointer m_BlockRadiusCalculator;

  SearchRegionImageSourcePointer m_SearchRegionImageSource;

private:
};

} // end namespace BlockMatching
} // end namespace itk

#include "itkBlockMatchingMultiResolutionIterationCommand.hxx"

#endif
