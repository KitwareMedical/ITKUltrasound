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
#ifndef itkBlockMatchingMultiResolutionSearchRegionImageSource_h
#define itkBlockMatchingMultiResolutionSearchRegionImageSource_h

#include "itkResampleImageFilter.h"

#include "itkImageDuplicator.h"
#include "itkImageSource.h"
#include "itkArray2D.h"

namespace itk
{
namespace BlockMatching
{

// Forward declaration.
template <typename TFixedImageF,
          typename TMovingImageF,
          typename TMetricImage,
          typename TDisplacementImageF,
          typename TCoordRep>
class MultiResolutionImageRegistrationMethod;


/** \class MultiResolutionSearchRegionImageSource
 *
 * \brief Generates the search region image during a multiresolution block
 * matching deformable image registration.
 *
 * This is a base class that is not intended to be instantiated.
 *
 * Subclases must be able to generate the search region for all levels create by the
 * MultiResolutionPyramidImageFilter in addition to the original image.
 *
 * \ingroup Ultrasound
 */
template <typename TFixedImage, typename TMovingImage, typename TDisplacementImage>
class ITK_TEMPLATE_EXPORT MultiResolutionSearchRegionImageSource
  : public ImageSource<Image<typename TMovingImage::RegionType, TMovingImage::ImageDimension>>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(MultiResolutionSearchRegionImageSource);

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

  /** Standard class type alias. */
  using Self = MultiResolutionSearchRegionImageSource;
  using Superclass = ImageSource<OutputImageType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** ScheduleType type alias support. */
  using PyramidScheduleType = Array2D<unsigned int>;

  /** OverlapScheduleType type alias support. */
  using OverlapScheduleType = Array2D<double>;

  /** Type of the displacement image from the previous level. */
  using DisplacementImageType = TDisplacementImage;
  using DisplacementImagePointer = typename DisplacementImageType::Pointer;

  /** Type of the filter used to resample the deformations. */
  using DisplacementResamplerType = ResampleImageFilter<DisplacementImageType, DisplacementImageType>;
  using DisplacementResamplerPointer = typename DisplacementResamplerType::Pointer;

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(MultiResolutionSearchRegionImageSource);

  /** Set the fixed image. */
  void
  SetFixedImage(FixedImageType * fixedImage)
  {
    m_FixedImage = fixedImage;
    this->Modified();
  }
  const FixedImageType *
  GetFixedImage() const
  {
    return this->m_FixedImage.GetPointer();
  }

  /** Set the moving image. */
  void
  SetMovingImage(MovingImageType * movingImage)
  {
    m_MovingImage = movingImage;
    this->Modified();
  }
  const MovingImageType *
  GetMovingImage() const
  {
    return this->m_MovingImage.GetPointer();
  }

  /** Set the fixed block radius, i.e. the radius of the matching kernel from
   * the fixed image. */
  virtual void
  SetFixedBlockRadius(const RadiusType & radius)
  {
    m_FixedBlockRadius = radius;
    this->Modified();
  }
  itkGetConstMacro(FixedBlockRadius, RadiusType);
  /** Set the fixed block radius to the the same in all directions. */
  virtual void
  SetFixedBlockRadius(const unsigned int rad)
  {
    RadiusType radius;
    radius.Fill(rad);
    this->SetFixedBlockRadius(radius);
  }

  /** SetPyramidSchedule() gets called with the pyramid schedule after the pyramid has
   * been generated.  This information is the available for child classes if
   * they choose to use it.  */
  virtual void
  SetPyramidSchedule(const PyramidScheduleType & schedule)
  {
    m_PyramidSchedule = schedule;
    this->Modified();
  }

  /** Get the multi-resolution schedule. */
  itkGetConstReferenceMacro(PyramidSchedule, PyramidScheduleType);

  /** Set/Get the overlap between fixed image blocks.  This value should be
   * greater than zero and defaults to unity.  Values less than unity will have
   * the blocks overlapping, e.g. 0.5 will render 50% overlap.  Values greater
   * than unity will result in spacing between blocks. The size of the schedule
   * should be the same as the size of the PyramidSchedule plus one.  */
  virtual void
  SetOverlapSchedule(const OverlapScheduleType & schedule)
  {
    m_OverlapSchedule = schedule;
    this->Modified();
  }
  itkGetConstReferenceMacro(OverlapSchedule, OverlapScheduleType);

  /** This is a convenience methods that sets the overlap to be the same across
   * a dimensions and levels. */
  virtual void
  SetOverlapSchedule(const double & schedule);

  /** So that it can call SetCurrentLevel(). */
  template <typename TFixedImageF,
            typename TMovingImageF,
            typename TMetricImage,
            typename TDisplacementImageF,
            typename TCoordRep>
  friend class MultiResolutionImageRegistrationMethod;

  itkGetConstObjectMacro(PreviousDisplacements, DisplacementImageType);

protected:
  using DisplacementDuplicatorType = ImageDuplicator<DisplacementImageType>;

  /** This is called by the MultiResolutionImageRegistration method to let this
   * filter know which level it is interested in. */
  itkSetMacro(CurrentLevel, unsigned long);

  /** This is called by the MultiResolutionImageRegistration method to allow
   * following search regions to be centered around the previous displacements.
   * */
  virtual void
  SetPreviousDisplacements(const DisplacementImageType * displacements);

  MultiResolutionSearchRegionImageSource();

  void
  BeforeThreadedGenerateData() override;

  void
  GenerateOutputInformation() override;

  typename FixedImageType::Pointer  m_FixedImage;
  typename MovingImageType::Pointer m_MovingImage;

  RadiusType m_FixedBlockRadius;

  PyramidScheduleType m_PyramidSchedule;
  OverlapScheduleType m_OverlapSchedule;

  unsigned long m_CurrentLevel;

  DisplacementImagePointer                     m_PreviousDisplacements;
  typename DisplacementDuplicatorType::Pointer m_DisplacementDuplicator;

  DisplacementResamplerPointer m_DisplacementResampler;

private:
};

} // namespace BlockMatching
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkBlockMatchingMultiResolutionSearchRegionImageSource.hxx"
#endif

#endif
