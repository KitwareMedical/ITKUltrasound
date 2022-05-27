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
#ifndef itkBlockMatchingMultiResolutionFixedSearchRegionImageSource_h
#define itkBlockMatchingMultiResolutionFixedSearchRegionImageSource_h

#include "itkBlockMatchingMultiResolutionSearchRegionImageSource.h"

namespace itk
{
namespace BlockMatching
{

/** \class MultiResolutionFixedSearchRegionImageSource
 *
 * \brief The search region has a fixed radius across levels.
 *
 * This class generates the search region with a fixed radius per level.
 *
 * The center of the search region is determined by the center of the tracking
 * kernel offset by the displacement from the previous level.
 *
 * \sa MultiResolutionSearchRegionImageSource
 *
 * \ingroup Ultrasound
 * */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementImage>
class ITK_TEMPLATE_EXPORT MultiResolutionFixedSearchRegionImageSource
  : public MultiResolutionSearchRegionImageSource<TFixedImage, TMovingImage, TDisplacementImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(MultiResolutionFixedSearchRegionImageSource);

  /** Standard class type alias. */
  using Self = MultiResolutionFixedSearchRegionImageSource;
  using Superclass = MultiResolutionSearchRegionImageSource<TFixedImage, TMovingImage, TDisplacementImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int, TMovingImage::ImageDimension);

  /** Type of the fixed image. */
  using FixedImageType = TFixedImage;
  using FixedRegionType = typename FixedImageType::RegionType;

  /** Type of the radius used to characterized the fixed image block. */
  using RadiusType = typename FixedImageType::SizeType;

  /** Type of the moving image. */
  using MovingImageType = TMovingImage;
  using MovingRegionType = typename MovingImageType::RegionType;

  /** Type of the search region image. */
  using OutputImageType = Image<typename MovingImageType::RegionType, ImageDimension>;
  using OutputRegionType = typename OutputImageType::RegionType;

  /** ScheduleType type alias support. */
  using PyramidScheduleType = typename Superclass::PyramidScheduleType;

  /** Type of the radius schedule. */
  using RadiusScheduleType = PyramidScheduleType;

  /** OverlapScheduleType type alias support. */
  using OverlapScheduleType = typename Superclass::OverlapScheduleType;

  /** Type of the displacement image from the previous level. */
  using DisplacementImageType = typename Superclass::DisplacementImageType;
  using DisplacementImagePointer = typename Superclass::DisplacementImagePointer;

  /** Type of the filter used to resample the deformations. */
  using DisplacementResamplerType = typename Superclass::DisplacementResamplerType;
  using DisplacementResamplerPointer = typename DisplacementResamplerType::Pointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(MultiResolutionFixedSearchRegionImageSource, MultiResolutionSearchRegionImageSource);

  /** New macro for creation of through a Smart Pointer. */
  itkNewMacro(Self);

  /** Set the search region, i.e. the radius of the search region in the
   * moving image. This is used across levels.
   * Therefore, the physical size of the search region scales with pyramid
   * schedule.
   * */
  virtual void
  SetSearchRegionRadiusSchedule(const RadiusType & radius);

  /** Set the search region radius the same in all directions. */
  virtual void
  SetSearchRegionRadiusSchedule(const unsigned int rad);

  virtual void
  SetSearchRegionRadiusSchedule(const RadiusScheduleType & schedule)
  {
    m_SearchRegionRadiusSchedule = schedule;
    m_SearchRegionRadiusSet = true;
    this->Modified();
  }

  itkGetConstReferenceMacro(SearchRegionRadiusSchedule, RadiusScheduleType);


protected:
  MultiResolutionFixedSearchRegionImageSource();

  void
  BeforeThreadedGenerateData() override;

  void
  DynamicThreadedGenerateData(const OutputRegionType & outputRegion) override;

  RadiusScheduleType m_SearchRegionRadiusSchedule;

  bool m_SearchRegionRadiusSet;

private:
};

} // namespace BlockMatching
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkBlockMatchingMultiResolutionFixedSearchRegionImageSource.hxx"
#endif

#endif
