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
#ifndef itkBlockMatchingMultiResolutionImageRegistrationMethod_h
#define itkBlockMatchingMultiResolutionImageRegistrationMethod_h

#include "itkImageSource.h"
#include "itkMultiResolutionPyramidImageFilter.h"

#include "itkBlockMatchingImageRegistrationMethod.h"
#include "itkBlockMatchingMultiResolutionBlockRadiusCalculator.h"
#include "itkBlockMatchingMultiResolutionSearchRegionImageSource.h"

namespace itk
{
namespace BlockMatching
{

/** \class MultiResolutionImageRegistrationMethod
 *
 * \brief Base class for multi-resolution image registration methods in the
 * BlockMatching set of tools.
 *
 * Before each resolution level an IterationEvent is invoked providing an
 * opportunity for a user interface to change any of the components,
 * change component parameters, or stop the registration.
 *
 * The SetNumberOfLevels() or SetSchedule() is used to set up the
 * MultiResolutionPyramidImageFilter.
 *
 * \sa ImageRegistrationMethod
 *
 * \ingroup RegistrationFilters
 * \ingroup Ultrasound
 * */
template <typename TFixedImage,
          typename TMovingImage,
          typename TMetricImage,
          typename TDisplacementImage,
          typename TCoordRep>
class ITK_TEMPLATE_EXPORT MultiResolutionImageRegistrationMethod : public ImageSource<TDisplacementImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(MultiResolutionImageRegistrationMethod);

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int, TDisplacementImage::ImageDimension);

  /** Type of the fixed image. */
  using FixedImageType = TFixedImage;
  using FixedRegionType = typename FixedImageType::RegionType;
  using FixedImagePointer = typename FixedImageType::Pointer;

  /** Type of the radius used to characterized the fixed image block. */
  using RadiusType = typename FixedImageType::SizeType;

  /** Type of the moving image. */
  using MovingImageType = TMovingImage;
  using MovingRegionType = typename MovingImageType::RegionType;
  using MovingImagePointer = typename MovingImageType::Pointer;

  /** Type of the metric image. */
  using MetricImageType = TMetricImage;

  /** Type of the displacement image. */
  using DisplacementImageType = TDisplacementImage;

  using RegionType = typename DisplacementImageType::RegionType;
  using IndexType = typename RegionType::IndexType;
  using SizeType = typename RegionType::SizeType;

  using SpacingType = typename DisplacementImageType::SpacingType;
  using DirectionType = typename DisplacementImageType::DirectionType;
  using OriginType = typename DisplacementImageType::PointType;

  /** Type of the search region image. */
  using SearchRegionImageType = Image<typename MovingImageType::RegionType, ImageDimension>;

  /** Standard class type alias. */
  using Self = MultiResolutionImageRegistrationMethod;
  using Superclass = ImageSource<TDisplacementImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MultiResolutionImageRegistrationMethod, ImageSource);

  /** Type of the Fixed image multiresolution pyramid. */
  using FixedImagePyramidType = MultiResolutionPyramidImageFilter<FixedImageType, FixedImageType>;
  using FixedImagePyramidPointer = typename FixedImagePyramidType::Pointer;

  /** Type of pyramid schedule type */
  using ScheduleType = typename FixedImagePyramidType::ScheduleType;

  /** Type of the moving image multiresolution pyramid. */
  using MovingImagePyramidType = MultiResolutionPyramidImageFilter<MovingImageType, MovingImageType>;
  using MovingImagePyramidPointer = typename MovingImagePyramidType::Pointer;

  /** Type of the registration method used at every level. */
  using ImageRegistrationMethodType = typename BlockMatching::
    ImageRegistrationMethod<TFixedImage, TMovingImage, TMetricImage, TDisplacementImage, TCoordRep>;
  using ImageRegistrationMethodPointer = typename ImageRegistrationMethodType::Pointer;

  /** Type of the class to calculate the fixed image matching kernel block
   * radius at every level. */
  using BlockRadiusCalculatorType = MultiResolutionBlockRadiusCalculator<TFixedImage>;
  using BlockRadiusCalculatorPointer = typename BlockRadiusCalculatorType::Pointer;

  using SearchRegionImageSourceType =
    typename BlockMatching::MultiResolutionSearchRegionImageSource<TFixedImage, TMovingImage, TDisplacementImage>;
  using SearchRegionImageSourcePointer = typename SearchRegionImageSourceType::Pointer;

  /** Method to stop the registration after registering a level. */
  void
  StopRegistration();

  /** Set/Get the Fixed image. */
  itkSetObjectMacro(FixedImage, FixedImageType);
  itkGetConstObjectMacro(FixedImage, FixedImageType);

  /** Set/Get the Moving image. */
  itkSetObjectMacro(MovingImage, MovingImageType);
  itkGetConstObjectMacro(MovingImage, MovingImageType);

  /** Set/Get the Fixed image pyramid. */
  itkSetObjectMacro(FixedImagePyramid, FixedImagePyramidType);
  itkGetConstObjectMacro(FixedImagePyramid, FixedImagePyramidType);

  /** Set/Get the Moving image pyramid. */
  itkSetObjectMacro(MovingImagePyramid, MovingImagePyramidType);
  itkGetConstObjectMacro(MovingImagePyramid, MovingImagePyramidType);

  /** Set/Get the schedules . */
  void
  SetSchedules(const ScheduleType & fixedSchedule, const ScheduleType & movingSchedule);
  itkGetConstMacro(FixedImagePyramidSchedule, ScheduleType);
  itkGetConstMacro(MovingImagePyramidSchedule, ScheduleType);

  /** Set/Get the number of multi-resolution levels. */
  void
  SetNumberOfLevels(unsigned long numberOfLevels);
  itkGetConstMacro(NumberOfLevels, unsigned long);

  /** Get the current resolution level being processed. */
  itkGetConstMacro(CurrentLevel, unsigned long);

  /** Method to return the latest modified time of this object or
   * any of its cached ivars */
  ModifiedTimeType
  GetMTime() const override;

  /** BlockMatching::ImageRegistrationMethod used to register each image at
   * every level. */
  itkSetObjectMacro(ImageRegistrationMethod, ImageRegistrationMethodType);
  itkGetModifiableObjectMacro(ImageRegistrationMethod, ImageRegistrationMethodType);

  /** Set the object used to generate the block radii in the fixed image at
   * every level. */
  itkSetObjectMacro(BlockRadiusCalculator, BlockRadiusCalculatorType);
  itkGetConstObjectMacro(BlockRadiusCalculator, BlockRadiusCalculatorType);

  /** Set the object used to generate the search regions. */
  itkSetObjectMacro(SearchRegionImageSource, SearchRegionImageSourceType);
  itkGetConstObjectMacro(SearchRegionImageSource, SearchRegionImageSourceType);

protected:
  MultiResolutionImageRegistrationMethod();
  ~MultiResolutionImageRegistrationMethod() override = default;

  /** The size and spacing of the search region image at the lowest level is
   * used to generate the information for the output image. */
  void
  GenerateOutputInformation() override;

  /** Generates the entire displacement image. */
  void
  EnlargeOutputRequestedRegion(DataObject * data) override
  {
    TDisplacementImage * output = this->GetOutput(0);
    output->SetRequestedRegionToLargestPossibleRegion();
  }

  /** Method invoked by the pipeline in order to trigger the computation of
   * the registration. */
  void
  GenerateData() override;

  /** Initialize by setting the interconnects between the components.
      This method is executed at every level of the pyramid with the
      values corresponding to this resolution
   */
  void
  Initialize();

  /** Create the image pyramids. */
  void
  PreparePyramids();

  /** Set up the fixed block radius calculator. */
  void
  PrepareBlockRadiusCalculator();

  /** Set up the search region image calculator. */
  void
  PrepareSearchRegionImageSource();

  /** Set the current level to be processed */
  itkSetMacro(CurrentLevel, unsigned long);

  FixedImagePointer  m_FixedImage;
  MovingImagePointer m_MovingImage;

  FixedImagePyramidPointer  m_FixedImagePyramid;
  MovingImagePyramidPointer m_MovingImagePyramid;

  unsigned long m_NumberOfLevels;
  unsigned long m_CurrentLevel;

  bool m_Stop;

  ScheduleType m_FixedImagePyramidSchedule;
  ScheduleType m_MovingImagePyramidSchedule;

  bool m_ScheduleSpecified;
  bool m_NumberOfLevelsSpecified;

  ImageRegistrationMethodPointer m_ImageRegistrationMethod;
  BlockRadiusCalculatorPointer   m_BlockRadiusCalculator;
  SearchRegionImageSourcePointer m_SearchRegionImageSource;

private:
};

} // end namespace BlockMatching
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkBlockMatchingMultiResolutionImageRegistrationMethod.hxx"
#endif

#endif
