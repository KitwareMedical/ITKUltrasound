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
#ifndef itkBlockMatchingMultiResolutionMinMaxSearchRegionImageSource_h
#define itkBlockMatchingMultiResolutionMinMaxSearchRegionImageSource_h

#include "itkBlockMatchingMultiResolutionSearchRegionImageSource.h"

namespace itk
{
namespace BlockMatching
{

/** \class MultiResolutionMinMaxSearchRegionImageSource
 *
 * \brief The search region is a factor of the matching block size.
 *
 * The search region size is a factor of the matching block size ( should be 1.0
 * or higher ).  This factor is set at the bottom level at top level and
 * linearly interpolated between levels.
 *
 * \todo rename this from 'MinMax' to 'TopBottomFactor'.
 *
 * \ingroup Ultrasound
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementImage>
class ITK_TEMPLATE_EXPORT MultiResolutionMinMaxSearchRegionImageSource
  : public MultiResolutionSearchRegionImageSource<TFixedImage, TMovingImage, TDisplacementImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(MultiResolutionMinMaxSearchRegionImageSource);

  /** Standard class type alias. */
  using Self = MultiResolutionMinMaxSearchRegionImageSource;
  using Superclass = MultiResolutionSearchRegionImageSource<TFixedImage, TMovingImage, TDisplacementImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information (and related methods). */
  itkTypeMacro(MultiResolutionMinMaxSearchRegionImageSource, MultiResolutionSearchRegionImageSource);

  /** New macro for creation of through a Smart Pointer. */
  itkNewMacro(Self);

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int, TMovingImage::ImageDimension);

  /** Type of the fixed image. */
  using FixedImageType = typename Superclass::FixedImageType;
  using FixedRegionType = typename Superclass::FixedRegionType;

  /** Type of the radius used to characterized the fixed image block. */
  using RadiusType = typename Superclass::RadiusType;

  /** Type of the moving image. */
  using MovingImageType = typename Superclass::MovingImageType;
  using MovingRegionType = typename Superclass::MovingRegionType;

  /** Type of the search region image. */
  using OutputImageType = typename Superclass::OutputImageType;
  using OutputRegionType = typename Superclass::OutputRegionType;

  /** Type of the displacement image. */
  using DisplacementImageType = typename Superclass::DisplacementImageType;

  /** ScheduleType type alias support. */
  using PyramidScheduleType = typename Superclass::PyramidScheduleType;

  /** Type of the search region to block radius ratio. */
  using FactorType = FixedArray<double, ImageDimension>;

  /** Set the ration of the search region radius to the matching block radius at
   * the bottom level. */
  void
  SetMinFactor(const FactorType & factor)
  {
    m_MinFactor = factor;
    this->Modified();
  }
  void
  SetMinFactor(const double & factor)
  {
    FactorType f;
    f.Fill(factor);
    this->SetMinFactor(f);
  }
  itkGetConstReferenceMacro(MinFactor, FactorType);

  /** Set the ration of the search region radius to the matching block radius at
   * the top level. */
  void
  SetMaxFactor(const FactorType & factor)
  {
    m_MaxFactor = factor;
    this->Modified();
  }
  void
  SetMaxFactor(const double & factor)
  {
    FactorType f;
    f.Fill(factor);
    this->SetMaxFactor(f);
  }
  itkGetConstReferenceMacro(MaxFactor, FactorType);

protected:
  MultiResolutionMinMaxSearchRegionImageSource()
    : m_MinFactor(1.1)
    , m_MaxFactor(3.0)
  {}

  void
  DynamicThreadedGenerateData(const OutputRegionType & outputRegion) override;

  FactorType m_MinFactor;
  FactorType m_MaxFactor;

private:
};

} // end namespace BlockMatching
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkBlockMatchingMultiResolutionMinMaxSearchRegionImageSource.hxx"
#endif

#endif
