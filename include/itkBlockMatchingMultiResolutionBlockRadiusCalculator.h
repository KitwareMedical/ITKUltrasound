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
#ifndef itkBlockMatchingMultiResolutionBlockRadiusCalculator_h
#define itkBlockMatchingMultiResolutionBlockRadiusCalculator_h

#include "itkArray2D.h"
#include "itkObject.h"

namespace itk
{
namespace BlockMatching
{

/** \class MultiResolutionBlockRadiusCalculator
 *
 * \brief Base class calculate the fixed image matching kernel radius for each level in a
 * BlockMatching::MultiResolutionImageRegistrationMethod.
 *
 * This must be able to produce the fixed image block radius for every level of
 * the MultiResolutionPyramidImage filter.
 *
 * \ingroup Ultrasound
 */
template <class TFixedImage>
class ITK_TEMPLATE_EXPORT MultiResolutionBlockRadiusCalculator : public Object
{
public:
  /** Standard class type alias. */
  using Self = MultiResolutionBlockRadiusCalculator;
  using Superclass = Object;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information (and related methods). */
  itkTypeMacro(MultiResolutionBlockRadiusCalculator, Object);

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int, TFixedImage::ImageDimension);

  /** Type of the fixed image type. */
  using FixedImageType = TFixedImage;
  using FixedImageConstPointer = typename FixedImageType::ConstPointer;
  using RadiusType = typename FixedImageType::SizeType;

  /** ScheduleType type alias support. */
  using ScheduleType = Array2D<unsigned int>;

  /** SetPyramidSchedule() gets called with the pyramid schedule after the pyramid has
   * been generated.  This information is the available for child classes if
   * they choose to use it.  */
  virtual void
  SetPyramidSchedule(const ScheduleType & schedule)
  {
    m_PyramidSchedule = schedule;
    this->Modified();
  }

  /** Get the multi-resolution schedule. */
  itkGetConstReferenceMacro(PyramidSchedule, ScheduleType);

  /** Set/Get the Fixed image. The fixed image is set during the
   * BlockMatching::MultiResolutionImageRegistrationMethod and is available for
   * child classes if the choose to use it.  */
  itkSetConstObjectMacro(FixedImage, FixedImageType);
  itkGetConstObjectMacro(FixedImage, FixedImageType);

  /** This abstract method must be implemented by child classes.  Given the
   * current level in the image pyramid, return the fixed block radius. */
  virtual const RadiusType &
  Compute(unsigned long current_level) = 0;

protected:
  MultiResolutionBlockRadiusCalculator() {}

  ScheduleType m_PyramidSchedule;

  FixedImageConstPointer m_FixedImage;

private:
  MultiResolutionBlockRadiusCalculator(const Self &);
  void
  operator=(const Self &);
};

} // end namespace BlockMatching
} // end namespace itk

#endif
