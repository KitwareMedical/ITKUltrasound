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
#ifndef itkBlockMatchingMetricImageFilter_h
#define itkBlockMatchingMetricImageFilter_h

#include "itkImageToImageFilter.h"

namespace itk
{
namespace BlockMatching
{

/** \class MetricImageFilter
 *
 * \brief Create a metric image by applying a metric to a kernel over a search
 * region.
 *
 * This class is intended to be used internally by
 * itk::BlockMatching::ImageRegistrationMethod and its ilk.
 *
 * There are two inputs to this image, the Fixed Image, and the Moving Image.
 * A block or kernel from the FixedImage is specified with SetFixedRegion().
 * A metric is evaluted by comparing the kernel to the search area specified by
 * SetMovingRegion().  The fixed image region should always be smaller than the
 * moving image region.  The metric is evaluated multiple times with different
 * translations applied such that a metric image is created.  The output of this
 * filter is the metric image.
 *
 * By default, the metric image spacing is the same as the moving
 * image's spacing.  The origin of the metric image is the same as the fixed image.
 *
 * Only derived classes of this class are intended to be instantiated.
 *
 *
 * \sa ImageRegistrationMethod
 *
 * \ingroup Ultrasound
 */
template <typename TFixedImage, typename TMovingImage, typename TMetricImage>
class ITK_TEMPLATE_EXPORT MetricImageFilter : public ImageToImageFilter<TFixedImage, TMetricImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(MetricImageFilter);

  /** Standard class type alias. */
  using Self = MetricImageFilter;
  using Superclass = ImageToImageFilter<TFixedImage, TMetricImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information (and related methods). */
  itkTypeMacro(MetricImageFilter, ImageToImageFilter);

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int, TFixedImage::ImageDimension);

  /**  Type of the Fixed image. */
  using FixedImageType = TFixedImage;
  using FixedImageRegionType = typename FixedImageType::RegionType;

  /**  Type of the Moving image. */
  using MovingImageType = TMovingImage;
  using MovingImageRegionType = typename MovingImageType::RegionType;

  /** Type of the Metric image. */
  using MetricImageType = TMetricImage;
  using MetricImageSpacingType = typename MetricImageType::SpacingType;
  using MetricImageRegionType = typename MetricImageType::RegionType;

  /** Type of the block radius. */
  using RadiusType = typename FixedImageType::SizeType;

  /** Set the fixed image. */
  void
  SetFixedImage(FixedImageType * fixedImage);

  /** Set the moving image. */
  void
  SetMovingImage(MovingImageType * movingImage);

  /** Set/Get the fixed image region.  This will defined the kernel in the fixed
   * image to be compared against the moving.  Every pixel in this region will
   * be compared against the moving image region to generate a pixel in the
   * output metric image. */
  void
  SetFixedImageRegion(const FixedImageRegionType & region);
  itkGetConstReferenceMacro(FixedImageRegion, FixedImageRegionType);

  /** Set/Get the moving image region.  This will defined the search region over
   * which the FixedImageRegion will be iteratively translated over to create
   * the metric image. */
  void
  SetMovingImageRegion(const MovingImageRegionType & region);
  itkGetConstReferenceMacro(MovingImageRegion, MovingImageRegionType);

protected:
  MetricImageFilter();
  ~MetricImageFilter() override = default;

  void
  GenerateOutputInformation() override;

  /** We set the fixed and moving image's requested region.  The fixed image's
   * requested region is equal to FixedImageRegion and the moving image's requested region is equal to the
   * MovingImageRegion dilated by a radius related to the FixedImageRegion's size. */
  void
  GenerateInputRequestedRegion() override;

  /** This filter produces the LargestPossibleRegion. */
  void
  EnlargeOutputRequestedRegion(DataObject * data) override;

  FixedImageRegionType  m_FixedImageRegion;
  MovingImageRegionType m_MovingImageRegion;

  bool m_FixedImageRegionDefined;
  bool m_MovingImageRegionDefined;

  // This is closely tied to m_FixedImageRegion.  It gets set when setting
  // FixedImageRegion.
  RadiusType m_FixedRadius;
  // This is again determined by m_FixedImageRegion, but it is the radius of the
  // block in the moving image, which may be different if the spacing is
  // different.
  RadiusType m_MovingRadius;

private:
};

} // end namespace BlockMatching
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkBlockMatchingMetricImageFilter.hxx"
#endif

#endif
