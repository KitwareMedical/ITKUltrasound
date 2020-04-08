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
#ifndef itkBlockMatchingStrainWindowBlockAffineTransformCommand_h
#define itkBlockMatchingStrainWindowBlockAffineTransformCommand_h

#include "itkCommand.h"

#include "itkStrainImageFilter.h"
#include "itkLinearLeastSquaresGradientImageFilter.h"

namespace itk
{
namespace BlockMatching
{

/** \class StrainWindowBlockAffineTransformCommand
 *
 * \brief Use the strain image calculated in the
 * StrainWindowDisplacementCalculator as the deforming strain image for the
 * BlockAffineTransformMetricImageFilter or calculate new strain with a linear
 * least squares gradient filter.
 *
 * To use, set the objects, then add as an observer to the
 * StrainWindowDisplacementCalculator.
 *
 * \ingroup Ultrasound
 */
template <typename TStrainWindowDisplacemenCalculator,
          typename TBlockAffineTransformMetricImageFilter,
          typename TStrainImageFilter>
class ITK_TEMPLATE_EXPORT StrainWindowBlockAffineTransformCommand : public Command
{
public:
  using Self = StrainWindowBlockAffineTransformCommand;
  using Superclass = Command;
  using Pointer = SmartPointer<Self>;

  itkNewMacro(Self);

  using StrainWindowDisplacementCalculatorType = TStrainWindowDisplacemenCalculator;

  using BlockAffineTransformMetricImageFilterType = TBlockAffineTransformMetricImageFilter;
  using StrainImageFilterType = TStrainImageFilter;

  void
  Execute(itk::Object * caller, const itk::EventObject & event) override
  {
    Execute((const itk::Object *)caller, event);
  }

  void
  Execute(const itk::Object * object, const itk::EventObject & event) override;

  itkSetObjectMacro(BlockAffineTransformMetricImageFilter, BlockAffineTransformMetricImageFilterType);

  /** Use the strain windower strain?  Otherwise the least squares base strain
   * filter is use.  Defaults to true. \@todo figure out why if set to true it
   * performs worse or segfaults.  */
  itkSetMacro(UseStrainWindowStrain, bool);
  itkGetConstMacro(UseStrainWindowStrain, bool);

protected:
  StrainWindowBlockAffineTransformCommand();

  using LeastSquaresFilterType =
    LinearLeastSquaresGradientImageFilter<typename StrainWindowDisplacementCalculatorType::MetricImageType,
                                          typename StrainWindowDisplacementCalculatorType::MetricPixelType,
                                          typename StrainWindowDisplacementCalculatorType::MetricPixelType>;

  using BlockAffineTransformMetricImageFilterPointer = typename BlockAffineTransformMetricImageFilterType::Pointer;

  BlockAffineTransformMetricImageFilterPointer m_BlockAffineTransformMetricImageFilter;

  typename LeastSquaresFilterType::Pointer m_LeastSquaresFilter;

  typename StrainImageFilterType::Pointer m_StrainImageFilter;

  bool m_UseStrainWindowStrain;

private:
  StrainWindowBlockAffineTransformCommand(const Self &);
  void
  operator=(const Self &);
};

} // end namespace BlockMatching
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkBlockMatchingStrainWindowBlockAffineTransformCommand.hxx"
#endif

#endif
