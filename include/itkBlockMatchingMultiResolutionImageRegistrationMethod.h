/*=========================================================================
 *
 *  Copyright Insight Software Consortium
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
template < typename TFixedImage, typename TMovingImage,
  typename TMetricImage, typename TDisplacementImage, typename TCoordRep >
class ITK_TEMPLATE_EXPORT MultiResolutionImageRegistrationMethod :
  public ImageSource< TDisplacementImage >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(MultiResolutionImageRegistrationMethod);

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int, TDisplacementImage::ImageDimension);

  /** Type of the fixed image. */
  typedef TFixedImage                         FixedImageType;
  typedef typename FixedImageType::RegionType FixedRegionType;
  typedef typename FixedImageType::Pointer    FixedImagePointer;

  /** Type of the radius used to characterized the fixed image block. */
  typedef typename FixedImageType::SizeType RadiusType;

  /** Type of the moving image. */
  typedef TMovingImage                         MovingImageType;
  typedef typename MovingImageType::RegionType MovingRegionType;
  typedef typename MovingImageType::Pointer    MovingImagePointer;

  /** Type of the metric image. */
  typedef TMetricImage  MetricImageType;

  /** Type of the displacement image. */
  typedef TDisplacementImage DisplacementImageType;

  typedef typename DisplacementImageType::RegionType RegionType;
  typedef typename RegionType::IndexType             IndexType;
  typedef typename RegionType::SizeType              SizeType;

  typedef typename DisplacementImageType::SpacingType   SpacingType;
  typedef typename DisplacementImageType::DirectionType DirectionType;
  typedef typename DisplacementImageType::PointType     OriginType;

  /** Type of the search region image. */
  typedef Image< typename MovingImageType::RegionType, ImageDimension > SearchRegionImageType;

  /** Standard class typedefs. */
  typedef MultiResolutionImageRegistrationMethod  Self;
  typedef ImageSource< TDisplacementImage >       Superclass;
  typedef SmartPointer< Self >                    Pointer;
  typedef SmartPointer< const Self >              ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MultiResolutionImageRegistrationMethod, ImageSource);

  /** Type of the Fixed image multiresolution pyramid. */
  typedef MultiResolutionPyramidImageFilter< FixedImageType,
                                             FixedImageType >
                                                   FixedImagePyramidType;
  typedef typename FixedImagePyramidType::Pointer  FixedImagePyramidPointer;

  /** Type of pyramid schedule type */
  typedef typename FixedImagePyramidType::ScheduleType ScheduleType;

  /** Type of the moving image multiresolution pyramid. */
  typedef MultiResolutionPyramidImageFilter< MovingImageType,
                                             MovingImageType >
                                                   MovingImagePyramidType;
  typedef typename MovingImagePyramidType::Pointer MovingImagePyramidPointer;

  /** Type of the registration method used at every level. */
  typedef typename BlockMatching::ImageRegistrationMethod< TFixedImage, TMovingImage,
          TMetricImage, TDisplacementImage, TCoordRep > ImageRegistrationMethodType;
  typedef typename ImageRegistrationMethodType::Pointer ImageRegistrationMethodPointer;

  /** Type of the class to calculate the fixed image matching kernel block
   * radius at every level. */
  typedef MultiResolutionBlockRadiusCalculator< TFixedImage > BlockRadiusCalculatorType;
  typedef typename BlockRadiusCalculatorType::Pointer         BlockRadiusCalculatorPointer;

  typedef typename BlockMatching::MultiResolutionSearchRegionImageSource< TFixedImage,
          TMovingImage, TDisplacementImage > SearchRegionImageSourceType;
  typedef typename SearchRegionImageSourceType::Pointer SearchRegionImageSourcePointer;

  /** Method to stop the registration after registering a level. */
  void StopRegistration();

  /** Set/Get the Fixed image. */
  itkSetObjectMacro( FixedImage, FixedImageType );
  itkGetConstObjectMacro( FixedImage, FixedImageType );

  /** Set/Get the Moving image. */
  itkSetObjectMacro( MovingImage, MovingImageType );
  itkGetConstObjectMacro( MovingImage, MovingImageType );

  /** Set/Get the Fixed image pyramid. */
  itkSetObjectMacro( FixedImagePyramid, FixedImagePyramidType );
  itkGetConstObjectMacro( FixedImagePyramid, FixedImagePyramidType );

  /** Set/Get the Moving image pyramid. */
  itkSetObjectMacro( MovingImagePyramid, MovingImagePyramidType );
  itkGetConstObjectMacro( MovingImagePyramid, MovingImagePyramidType );

  /** Set/Get the schedules . */
  void SetSchedules( const ScheduleType & fixedSchedule,
                    const ScheduleType & movingSchedule );
  itkGetConstMacro( FixedImagePyramidSchedule, ScheduleType );
  itkGetConstMacro( MovingImagePyramidSchedule, ScheduleType );

  /** Set/Get the number of multi-resolution levels. */
  void SetNumberOfLevels( unsigned long numberOfLevels );
  itkGetConstMacro( NumberOfLevels, unsigned long );

  /** Get the current resolution level being processed. */
  itkGetConstMacro( CurrentLevel, unsigned long );

  /** Method to return the latest modified time of this object or
   * any of its cached ivars */
  unsigned long GetMTime() const override;

  /** BlockMatching::ImageRegistrationMethod used to register each image at
   * every level. */
  itkSetObjectMacro( ImageRegistrationMethod, ImageRegistrationMethodType );
  itkGetModifiableObjectMacro( ImageRegistrationMethod, ImageRegistrationMethodType );

  /** Set the object used to generate the block radii in the fixed image at
   * every level. */
  itkSetObjectMacro( BlockRadiusCalculator, BlockRadiusCalculatorType );
  itkGetConstObjectMacro( BlockRadiusCalculator, BlockRadiusCalculatorType );

  /** Set the object used to generate the search regions. */
  itkSetObjectMacro( SearchRegionImageSource, SearchRegionImageSourceType );
  itkGetConstObjectMacro( SearchRegionImageSource, SearchRegionImageSourceType );

protected:
  MultiResolutionImageRegistrationMethod();
  virtual ~MultiResolutionImageRegistrationMethod() {};

  /** The size and spacing of the search region image at the lowest level is
   * used to generate the information for the output image. */
  void GenerateOutputInformation() override;

  /** Generates the entire displacement image. */
  void EnlargeOutputRequestedRegion( DataObject* data ) override
    {
    TDisplacementImage* output = this->GetOutput( 0 );
    output->SetRequestedRegionToLargestPossibleRegion();
    }

  /** Method invoked by the pipeline in order to trigger the computation of
   * the registration. */
  void GenerateData() override;

  /** Initialize by setting the interconnects between the components.
      This method is executed at every level of the pyramid with the
      values corresponding to this resolution
   */
  void Initialize();

  /** Create the image pyramids. */
  void PreparePyramids();

  /** Set up the fixed block radius calculator. */
  void PrepareBlockRadiusCalculator();

  /** Set up the search region image calculator. */
  void PrepareSearchRegionImageSource();

  /** Set the current level to be processed */
  itkSetMacro( CurrentLevel, unsigned long );

  FixedImagePointer                m_FixedImage;
  MovingImagePointer               m_MovingImage;

  FixedImagePyramidPointer         m_FixedImagePyramid;
  MovingImagePyramidPointer        m_MovingImagePyramid;

  unsigned long                    m_NumberOfLevels;
  unsigned long                    m_CurrentLevel;

  bool                             m_Stop;

  ScheduleType                     m_FixedImagePyramidSchedule;
  ScheduleType                     m_MovingImagePyramidSchedule;

  bool                             m_ScheduleSpecified;
  bool                             m_NumberOfLevelsSpecified;

  ImageRegistrationMethodPointer   m_ImageRegistrationMethod;
  BlockRadiusCalculatorPointer     m_BlockRadiusCalculator;
  SearchRegionImageSourcePointer   m_SearchRegionImageSource;

private:
};

} // end namespace BlockMatching
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlockMatchingMultiResolutionImageRegistrationMethod.hxx"
#endif

#endif
