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
#ifndef itkBlockMatchingMultiResolutionMinMaxBlockRadiusCalculator_h
#define itkBlockMatchingMultiResolutionMinMaxBlockRadiusCalculator_h

#include "itkObjectFactory.h"
#include "itkBlockMatchingMultiResolutionBlockRadiusCalculator.h"

namespace itk
{
namespace BlockMatching
{
/** \class MultiResolutionMinMaxBlockRadiusCalculator
 *
 * \brief The radius varies linearly from a minimum at the lowest level to a
 * maximum at the highest level.
 *
 * \todo rename this from 'MinMax' to 'TopBottom'.
 *
 * \sa MultiResolutionBlockRadiusCalculator
 *
 * \ingroup Ultrasound
 */
template <class TFixedImage>
class ITK_TEMPLATE_EXPORT MultiResolutionMinMaxBlockRadiusCalculator
  : public MultiResolutionBlockRadiusCalculator<TFixedImage>
{
public:
  /** Standard class type alias. */
  using Self = MultiResolutionMinMaxBlockRadiusCalculator;
  using Superclass = MultiResolutionBlockRadiusCalculator<TFixedImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information ( and related methods ). */
  itkTypeMacro(MultiResolutionMinMaxBlockRadiusCalculator, MultiResolutionBlockRadiusCalculator);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Type of the fixed image. */
  using FixedImageType = typename Superclass::FixedImageType;
  using RadiusType = typename Superclass::RadiusType;

  /** Set the radius used for at the bottom level. */
  virtual void
  SetMinRadius(const RadiusType & radius)
  {
    m_MinRadius = radius;
  }
  itkGetConstReferenceMacro(MinRadius, RadiusType);

  /** Set the radius used at the top level. */
  virtual void
  SetMaxRadius(const RadiusType & radius)
  {
    m_MaxRadius = radius;
  }
  itkGetConstReferenceMacro(Radius, RadiusType);

  const RadiusType &
  Compute(unsigned long level) override
  {
    double slope;
    double distance = static_cast<double>(this->m_PyramidSchedule.rows() - 1);
    for (unsigned int i = 0; i < this->m_PyramidSchedule.cols(); ++i)
    {
      slope = (static_cast<double>(m_MinRadius[i]) - static_cast<double>(m_MaxRadius[i])) / distance;
      m_Radius[i] = static_cast<typename RadiusType::SizeValueType>(slope * level + m_MaxRadius[i]);
    }
    return m_Radius;
  }

protected:
  MultiResolutionMinMaxBlockRadiusCalculator()
  {
    m_MinRadius.Fill(1);
    m_MaxRadius.Fill(1);
  }

  RadiusType m_MinRadius;
  RadiusType m_MaxRadius;

  RadiusType m_Radius;

private:
  MultiResolutionMinMaxBlockRadiusCalculator(const Self &);
  void
  operator=(const Self &);
};

} // namespace BlockMatching
} // namespace itk


#endif
