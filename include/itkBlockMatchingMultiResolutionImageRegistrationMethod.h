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
template < class TFixedImage, class TMovingImage,
  class TMetricImage, class TDisplacementImage, class TCoordRep >
class ITK_EXPORT MultiResolutionImageRegistrationMethod :
  public ImageSource< TDisplacementImage >
{
public:
  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TDisplacementImage::ImageDimension);
  /** Type of the fixed image. */
  typedef TFixedImage                             FixedImageType;
  typedef typename FixedImageType::RegionType     FixedRegionType;
  typedef typename FixedImageType::Pointer   FixedImagePointer;

  /** Type of the radius used to characterized the fixed image block. */
  typedef typename FixedImageType::SizeType RadiusType;

  /** Type of the moving image. */
  typedef TMovingImage                            MovingImageType;
  typedef typename MovingImageType::RegionType    MovingRegionType;
  typedef typename MovingImageType::Pointer  MovingImagePointer;

  /** Type of the metric image. */
  typedef TMetricImage  MetricImageType;

  /** Type of the displacement image. */
  typedef TDisplacementImage DisplacementImageType;

  typedef typename DisplacementImageType::RegionType RegionType;
  typedef typename RegionType::IndexType IndexType;
  typedef typename RegionType::SizeType  SizeType;

  typedef typename DisplacementImageType::SpacingType   SpacingType;
  typedef typename DisplacementImageType::DirectionType DirectionType;
  typedef typename DisplacementImageType::PointType     OriginType;

  /** Type of the search region image. */
  typedef typename ::itk::Image< typename MovingImageType::RegionType, 
    ImageDimension > SearchRegionImageType;

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
  typedef typename ::itk::BlockMatching::ImageRegistrationMethod< TFixedImage, TMovingImage,
          TMetricImage, TDisplacementImage, TCoordRep > ImageRegistrationMethodType;
  typedef typename ImageRegistrationMethodType::Pointer ImageRegistrationMethodPointer;

  /** Type of the class to calculate the fixed image matching kernel block
   * radius at every level. */
  typedef typename ::itk::BlockMatching::MultiResolutionBlockRadiusCalculator< TFixedImage >
    BlockRadiusCalculatorType;
  typedef typename BlockRadiusCalculatorType::Pointer   BlockRadiusCalculatorPointer;

  typedef typename ::itk::BlockMatching::MultiResolutionSearchRegionImageSource< TFixedImage,
          TMovingImage, TDisplacementImage > SearchRegionImageSourceType;
  typedef typename SearchRegionImageSourceType::Pointer SearchRegionImageSourcePointer;

  /** Method to stop the registration after registering a level. */
  void StopRegistration();

  /** Set/Get the Fixed image. */
  itkSetObjectMacro( FixedImage, FixedImageType );
  itkGetObjectMacro( FixedImage, FixedImageType ); 

  /** Set/Get the Moving image. */
  itkSetObjectMacro( MovingImage, MovingImageType );
  itkGetObjectMacro( MovingImage, MovingImageType );

  /** Set/Get the Fixed image pyramid. */
  itkSetObjectMacro( FixedImagePyramid, FixedImagePyramidType );
  itkGetObjectMacro( FixedImagePyramid, FixedImagePyramidType ); 

  /** Set/Get the Moving image pyramid. */
  itkSetObjectMacro( MovingImagePyramid, MovingImagePyramidType );
  itkGetObjectMacro( MovingImagePyramid, MovingImagePyramidType );

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
  unsigned long GetMTime() const;  

  /** BlockMatching::ImageRegistrationMethod used to register each image at
   * every level. */
  itkSetObjectMacro( ImageRegistrationMethod, ImageRegistrationMethodType );
  itkGetObjectMacro( ImageRegistrationMethod, ImageRegistrationMethodType );

  /** Set the object used to generate the block radii in the fixed image at
   * every level. */
  itkSetObjectMacro( BlockRadiusCalculator, BlockRadiusCalculatorType );
  itkGetObjectMacro( BlockRadiusCalculator, BlockRadiusCalculatorType );

  /** Set the object used to generate the search regions. */
  itkSetObjectMacro( SearchRegionImageSource, SearchRegionImageSourceType );
  itkGetObjectMacro( SearchRegionImageSource, SearchRegionImageSourceType );

protected:
  MultiResolutionImageRegistrationMethod();
  virtual ~MultiResolutionImageRegistrationMethod() {};

  /** The size and spacing of the search region image at the lowest level is
   * used to generate the information for the output image. */
  virtual void GenerateOutputInformation();

  /** Generates the entire displacement image. */
  virtual void EnlargeOutputRequestedRegion( DataObject* data )
    {
    TDisplacementImage* output = this->GetOutput( 0 );
    output->SetRequestedRegionToLargestPossibleRegion();
    }

  /** Method invoked by the pipeline in order to trigger the computation of 
   * the registration. */
  virtual void GenerateData ();

  /** Initialize by setting the interconnects between the components.
      This method is executed at every level of the pyramid with the
      values corresponding to this resolution
   */
  void Initialize() throw (ExceptionObject);

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
  MultiResolutionImageRegistrationMethod( const Self& ); // purposely not implemented
  void operator=( const Self& ); //purposely not implemented
};

} // end namespace BlockMatching
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlockMatchingMultiResolutionImageRegistrationMethod.hxx"
#endif

#endif
