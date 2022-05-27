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
#ifndef itkBlockMatchingMultiResolutionThresholdBoundingBoxSearchRegionImageSource_h
#define itkBlockMatchingMultiResolutionThresholdBoundingBoxSearchRegionImageSource_h

#include "itkArray.h"
#include "itkResampleImageFilter.h"
#include "itkVector.h"

#include "itkBlockMatchingMultiResolutionSearchRegionImageSource.h"
#include "itkBlockMatchingMetricImageToDisplacementCalculator.h"

namespace itk
{
namespace BlockMatching
{

/** \class MultiResolutionThresholdBoundingBoxSearchRegionImageSource
 *
 * \brief The search region at subsequent levels is determined by bounding box
 * of the metric image convex hull at the current level.
 *
 * This is a SearchRegionImageSource that defines the search region size
 * dynamically according the metric image performance at higher levels.  The
 * search region given the metric images at the current level is determined by
 * bounding box of the convex hull of metric values above a given threshold.
 * This size is interpolated is interpolated at the following levels.
 *
 * In order to examine the metric image values, this class also inherits from
 * MetricImageToDisplacementCalculator.  The center of the search region is
 * taken to be the displacement.
 *
 * \ingroup Ultrasound
 * */
template <class TFixedImage, class TMovingImage, class TMetricImage, class TDisplacementImage>
class ITK_TEMPLATE_EXPORT MultiResolutionThresholdBoundingBoxSearchRegionImageSource
  : public MultiResolutionSearchRegionImageSource<TFixedImage, TMovingImage, TDisplacementImage>
  , public MetricImageToDisplacementCalculator<TMetricImage, TDisplacementImage>
{
public:
  /** Standard class type alias. */
  using Self = MultiResolutionThresholdBoundingBoxSearchRegionImageSource;
  using Superclass = MultiResolutionSearchRegionImageSource<TFixedImage, TMovingImage, TDisplacementImage>;
  using DisplacementCalculatorSuperclass = MetricImageToDisplacementCalculator<TMetricImage, TDisplacementImage>;

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
  using OutputImageType = typename Superclass::OutputImageType;
  using OutputRegionType = typename OutputImageType::RegionType;

  using SearchRegionRadiusImageType =
    Image<typename itk::Vector<typename RadiusType::SizeValueType, ImageDimension>, ImageDimension>;
  using SearchRegionRadiusImagePointer = typename SearchRegionRadiusImageType::Pointer;

  /** Type of the filter used to resample the search region images. */
  using SearchRegionRadiusResamplerType = ResampleImageFilter<SearchRegionRadiusImageType, SearchRegionRadiusImageType>;
  using SearchRegionRadiusResamplerPointer = typename SearchRegionRadiusResamplerType::Pointer;

  /** ScheduleType type alias support. */
  using PyramidScheduleType = typename Superclass::PyramidScheduleType;

  /** Type of the threshold. */
  using ThresholdType = typename TMetricImage::PixelType;
  using ThresholdScheduleType = Array<ThresholdType>;

  /** OverlapScheduleType type alias support. */
  using OverlapScheduleType = typename Superclass::OverlapScheduleType;

  /** Type of the displacement image from the previous level. */
  using DisplacementImageType = TDisplacementImage;

  /** Type of the filter used to resample the deformations. */
  using DisplacementResamplerType = ResampleImageFilter<DisplacementImageType, DisplacementImageType>;
  using DisplacementResamplerPointer = typename DisplacementResamplerType::Pointer;

  /** Types inherited from the DisplacementCalculator superclass. */
  using MetricImageType = typename DisplacementCalculatorSuperclass::MetricImageType;
  using PointType = typename DisplacementCalculatorSuperclass::PointType;
  using IndexType = typename DisplacementCalculatorSuperclass::IndexType;

  /** Run-time type information (and related methods). */
  itkTypeMacro(MultiResolutionThresholdBoundingBoxSearchRegionImageSource, MultiResolutionFixedSearchRegionImageSource);

  /** New macro for creation of through a Smart Pointer is not used because of
   * ambiguities with LightObject. */
  static Pointer
  New(void)
  {
    Pointer smartPtr = ObjectFactory<Self>::Create();
    if (smartPtr.GetPointer() == nullptr)
    {
      smartPtr = new Self;
    }
    smartPtr->UnRegister();
    return smartPtr;
  }
  // virtual LightObject::Pointer CreateAnother(void) const
  //{
  // LightObject::Pointer smartPtr;
  //// @todo fix me -- itk::LightObject ambiguity somewhere in the following
  // line
  // smartPtr = Self::New().GetPointer();
  // return smartPtr;
  //}

  /** Set/get the threshold schedule.  The search region size in the next level will be
   * interpolated from the size of the area above this threshold.  The threshold
   * is an Array that should have the same size as the number of pyramid levels less one.  */
  virtual void
  SetThresholdSchedule(const ThresholdType threshold);

  virtual void
  SetThresholdSchedule(const ThresholdScheduleType & schedule)
  {
    m_ThresholdSchedule = schedule;
    m_ThresholdScheduleSpecified = true;
    this->Modified();
  }

  itkGetConstReferenceMacro(ThresholdSchedule, ThresholdScheduleType);

  /** Set/Get the search region radius at the top, highest level ( level 0 ). */
  virtual void
  SetTopLevelRadius(const RadiusType & radius)
  {
    m_TopLevelRadius = radius;
    m_TopLevelRadiusSpecified = true;
    this->Modified();
  }

  /** Set/Get the minimum search region radius.  The search region radius is not
   * allowed to go smaller that the given value. */
  virtual void
  SetMinimumSearchRegionRadius(const RadiusType & radius)
  {
    m_MinimumSearchRegionRadius = radius;
    this->Modified();
  }
  /** Set the radius to the given value in all directions. */
  virtual void
  SetMinimumSearchRegionRadius(const unsigned int rad)
  {
    RadiusType radius;
    radius.Fill(rad);
    this->SetMinimumSearchRegionRadius(radius);
  }
  itkGetConstReferenceMacro(MinimumSearchRegionRadius, RadiusType);

  /** We allocate the previous search region image. */
  virtual void
  SetDisplacementImage(DisplacementImageType * image);

  /** Calculates the search region radius based on the metric image. */
  virtual void
  SetMetricImagePixel(const PointType & point, const IndexType & index, MetricImageType * image);

  virtual void
  Compute()
  {
    // We do this here instead of SetMetricImagePixel so it only has to be done
    // once.
    this->m_DisplacementImage->Modified();
  }

  /** This is needed to avoid resolution ambiguities that occur with multiple
   * inheritance. */
  virtual void
  Register() const
  {
    Superclass::Register();
  }

  virtual void
  UnRegister() const
  {
    Superclass::UnRegister();
  }

  virtual void
  Modified() const
  {
    Superclass::Modified();
  }

  bool
  GetDebug() const
  {
    return Superclass::GetDebug();
  }

  itkGetConstObjectMacro(SearchRegionRadiusImage, SearchRegionRadiusImageType);

protected:
  MultiResolutionThresholdBoundingBoxSearchRegionImageSource();

  virtual void
  BeforeThreadedGenerateData() override;

  virtual void
  ThreadedGenerateData(const OutputRegionType & outputRegion, ThreadIdType threadID) override;

  ThresholdScheduleType m_ThresholdSchedule;
  bool                  m_ThresholdScheduleSpecified;
  RadiusType            m_TopLevelRadius;
  bool                  m_TopLevelRadiusSpecified;

  RadiusType m_MinimumSearchRegionRadius;

  DisplacementResamplerPointer m_DisplacementResampler;

  // The search region radius from the previous level.
  SearchRegionRadiusImagePointer     m_SearchRegionRadiusImage;
  SearchRegionRadiusResamplerPointer m_SearchRegionRadiusResampler;

private:
  MultiResolutionThresholdBoundingBoxSearchRegionImageSource(const Self &);
  void
  operator=(const Self &);
};

} // namespace BlockMatching
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkBlockMatchingMultiResolutionThresholdBoundingBoxSearchRegionImageSource.hxx"
#endif

#endif
