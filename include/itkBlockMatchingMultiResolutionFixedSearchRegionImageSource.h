#ifndef __itkBlockMatchingMultiResolutionFixedSearchRegionImageSource_h
#define __itkBlockMatchingMultiResolutionFixedSearchRegionImageSource_h

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
 * */
template < class TFixedImage, class TMovingImage, class TDisplacementImage >
class ITK_EXPORT MultiResolutionFixedSearchRegionImageSource :
  public MultiResolutionSearchRegionImageSource < TFixedImage, TMovingImage,
                                                  TDisplacementImage >
{
public:
  /** Standard class typedefs. */
  typedef MultiResolutionFixedSearchRegionImageSource    Self;
  typedef MultiResolutionSearchRegionImageSource< TFixedImage,
          TMovingImage, TDisplacementImage >             Superclass;
  typedef SmartPointer< Self >                           Pointer;
  typedef SmartPointer< const Self >                     ConstPointer;

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TMovingImage::ImageDimension);

  /** Type of the fixed image. */
  typedef TFixedImage                         FixedImageType;
  typedef typename FixedImageType::RegionType FixedRegionType;

  /** Type of the radius used to characterized the fixed image block. */
  typedef typename FixedImageType::SizeType RadiusType;

  /** Type of the moving image. */
  typedef TMovingImage                         MovingImageType;
  typedef typename MovingImageType::RegionType MovingRegionType;

  /** Type of the search region image. */
  typedef typename::itk::Image< typename MovingImageType::RegionType,
                                ImageDimension > OutputImageType;
  typedef typename OutputImageType::RegionType OutputRegionType;

  /** ScheduleType typedef support. */
  typedef typename Superclass::PyramidScheduleType PyramidScheduleType;

  /** Type of the radius schedule. */
  typedef PyramidScheduleType RadiusScheduleType;

  /** OverlapScheduleType typedef support. */
  typedef typename Superclass::OverlapScheduleType OverlapScheduleType;

  /** Type of the displacement image from the previous level. */
  typedef typename Superclass::DisplacementImageType    DisplacementImageType;
  typedef typename Superclass::DisplacementImagePointer DisplacementImagePointer;

  /** Type of the filter used to resample the deformations. */
  typedef typename Superclass::DisplacementResamplerType DisplacementResamplerType;
  typedef typename DisplacementResamplerType::Pointer    DisplacementResamplerPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( MultiResolutionFixedSearchRegionImageSource, MultiResolutionSearchRegionImageSource );

  /** New macro for creation of through a Smart Pointer. */
  itkNewMacro( Self );

  /** Set the search region, i.e. the radius of the search region in the
   * moving image. This is used across levels.
   * Therefore, the physical size of the search region scales with pyramid
   * schedule.
   * */
  virtual void SetSearchRegionRadiusSchedule( const RadiusType& radius );

  /** Set the search region radius the same in all directions. */
  virtual void SetSearchRegionRadiusSchedule( const unsigned int rad );

  virtual void SetSearchRegionRadiusSchedule( const RadiusScheduleType& schedule )
    {
    m_SearchRegionRadiusSchedule = schedule;
    m_SearchRegionRadiusSet      = true;
    this->Modified();
    }

  itkGetConstReferenceMacro( SearchRegionRadiusSchedule, RadiusScheduleType );


protected:
  MultiResolutionFixedSearchRegionImageSource();

  virtual void BeforeThreadedGenerateData();

  virtual void ThreadedGenerateData( const OutputRegionType& outputRegion,
                                     int threadID );

  RadiusScheduleType m_SearchRegionRadiusSchedule;

  bool m_SearchRegionRadiusSet;

private:
  MultiResolutionFixedSearchRegionImageSource( const Self& );
  void operator=( const Self& );
};

} // end namespace itk
} // end namespace BlockMatching

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlockMatchingMultiResolutionFixedSearchRegionImageSource.txx"
#endif

#endif

