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
#ifndef itkBlockMatchingImageRegistrationMethod_h
#define itkBlockMatchingImageRegistrationMethod_h

#include "itkImageToImageFilter.h"

#include "itkBlockMatchingMetricImageFilter.h"
#include "itkBlockMatchingMetricImageToDisplacementCalculator.h"

namespace itk
{
namespace BlockMatching
{

/** \class ImageRegistrationMethod
 *
 * This class defines the interface to perform deformable image registration by
 * block matching.
 *
 * Displacements are calculated at every block from the FixedImage to the Moving
 * Image.
 *
 * Blocks are neighborhoods with a fixed radius, and they are located on a grid
 * in the fixed image.
 *
 * An Image of search Regions in the moving image specifies each block's search
 * area.  The information from the search region image (origin, spacing, region,
 * etc) determines the information in the output displacement image.
 *
 * \sa ImageRegistrationMethod
 *
 * \ingroup RegistrationFilters
 * \ingroup Ultrasound
 */
template <typename TFixedImage,
          typename TMovingImage,
          typename TMetricImage,
          typename TDisplacementImage,
          typename TCoordRep>
class ITK_TEMPLATE_EXPORT ImageRegistrationMethod
  : public ImageToImageFilter<itk::Image<typename TMovingImage::RegionType, TDisplacementImage::ImageDimension>,
                              TDisplacementImage>
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(ImageRegistrationMethod);

  /** ImageDimension enumeration. */
  static constexpr unsigned int ImageDimension = TDisplacementImage::ImageDimension;

  /** Type of the fixed image. */
  using FixedImageType = TFixedImage;
  using FixedRegionType = typename FixedImageType::RegionType;

  /** Type of the radius used to characterized the fixed image block. */
  using RadiusType = typename FixedImageType::SizeType;

  /** Type of the moving image. */
  using MovingImageType = TMovingImage;
  using MovingRegionType = typename MovingImageType::RegionType;

  /** Type of the metric image. */
  using MetricImageType = TMetricImage;

  /** Type of the displacement image. */
  using ImageType = TDisplacementImage;

  using RegionType = typename ImageType::RegionType;
  using IndexType = typename RegionType::IndexType;
  using SizeType = typename RegionType::SizeType;

  using SpacingType = typename ImageType::SpacingType;
  using DirectionType = typename ImageType::DirectionType;
  using OriginType = typename ImageType::PointType;

  /** Type of the search region image. */
  using SearchRegionImageType = Image<typename MovingImageType::RegionType, ImageDimension>;

  /** Standard class type alias. */
  using Self = ImageRegistrationMethod;
  using Superclass = ImageToImageFilter<SearchRegionImageType, TDisplacementImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(ImageRegistrationMethod);

  /** Type of the point use for determing the location in the fixed image of a
   * block's center. */
  using CoordRepType = typename itk::Point<TCoordRep, ImageDimension>;

  /** Type of the MetricImageFilter. */
  using MetricImageFilterType = MetricImageFilter<FixedImageType, MovingImageType, MetricImageType>;

  /** Type of the MetricImageToDisplacementCalculator. */
  using MetricImageToDisplacementCalculatorType = MetricImageToDisplacementCalculator<TMetricImage, TDisplacementImage>;

  /** Set the fixed image. */
  void
  SetFixedImage(FixedImageType * fixedImage);
  const FixedImageType *
  GetFixedImage() const
  {
    return this->m_FixedImage.GetPointer();
  }

  /** Set the moving image. */
  void
  SetMovingImage(MovingImageType * movingImage);
  const MovingImageType *
  GetMovingImage() const
  {
    return this->m_MovingImage.GetPointer();
  }

  /** Set the MetricImageFilter. */
  itkSetObjectMacro(MetricImageFilter, MetricImageFilterType);
  itkGetConstObjectMacro(MetricImageFilter, MetricImageFilterType);

  /** Set/Get the MetricImageToDisplacementCalculator.  This defaults to a
   * MaximumPixelDisplacementCalculator. */
  itkSetObjectMacro(MetricImageToDisplacementCalculator, MetricImageToDisplacementCalculatorType);
  itkGetConstObjectMacro(MetricImageToDisplacementCalculator, MetricImageToDisplacementCalculatorType);

  /** Whether or not to use streaming.  Streaming is achieved by streaming each
   * block matching.  This allows for deformable registration of very large
   * images, but it comes at a performance penalty.  By default it is OFF. */
  itkSetMacro(UseStreaming, bool);
  itkGetConstMacro(UseStreaming, bool);
  itkBooleanMacro(UseStreaming);

  /** Set the radius for blocks in the fixed image to be matched against the
   * moving image.  This is a radius defined similarly to an itk::Neighborhood
   * radius, i.e., the size of the block in the i'th direction is 2*radius[i] +
   * 1.  Every fixed image block to be registered uses the same radius.
   */
  virtual void
  SetRadius(const RadiusType & radius)
  {
    m_Radius = radius;
    this->Modified();
  }
  /** Set the radius to the given value in all directions. */
  virtual void
  SetRadius(const unsigned int rad)
  {
    RadiusType radius;
    radius.Fill(rad);
    this->SetRadius(radius);
  }
  itkGetConstReferenceMacro(Radius, RadiusType);

  /** Set/Get the search region image.  The SearchRegionImage has the same
   * LargestPossibleRegion as the output displacement image.  It contains
   * ImageRegions in the moving image that define the search region for each
   * block in the fixed image.
   *
   * The metric image is created by calculating the
   * value of the metric between the fixed image block and the corresponding
   * area/volume in the moving image.  The center of the fixed image block is
   * translated between the corners of the given search region, evaluating the
   * metric at the spacing of the metric image.  Therefore, contributing region
   * in the moving image to the search is actually the given search region
   * dilated by the radius of the kernel block. */
  virtual void
  SetSearchRegionImage(SearchRegionImageType * searchRegionImage)
  {
    this->SetInput(searchRegionImage);
  }


protected:
  ImageRegistrationMethod();
  ~ImageRegistrationMethod() override = default;

  void
  GenerateOutputInformation() override;

  void
  GenerateInputRequestedRegion() override;

  void
  EnlargeOutputRequestedRegion(DataObject * data) override;

  /** Initialize by setting the interconnects between the components. */
  virtual void
  Initialize();

  void
  GenerateData() override;

  typename FixedImageType::Pointer  m_FixedImage;
  typename MovingImageType::Pointer m_MovingImage;

  typename MetricImageFilterType::Pointer                   m_MetricImageFilter;
  typename MetricImageToDisplacementCalculatorType::Pointer m_MetricImageToDisplacementCalculator;

  bool       m_UseStreaming;
  RadiusType m_Radius;

private:
};

} // end namespace BlockMatching
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkBlockMatchingImageRegistrationMethod.hxx"
#endif

#endif
