#ifndef __itkBlockMatchingSearchRegionImageInitializer_h
#define __itkBlockMatchingSearchRegionImageInitializer_h

#include "itkImageSource.h"

namespace itk
{
namespace BlockMatching
{

/** \class SearchRegionImageInitializer
 *
 * \brief Creates a SearchRegionImage for input into a
 * BlockMatching::ImageRegistrationMethod.
 *
 * This creates a SearchRegionImage for input into the
 * BlockMatching::ImageRegistrationMethod.  The search regions are centered
 * around the fixed image blocks and evenly spaced.  Overlap between blocks may
 * be set with SetOverlap().  The input fixed image, fixed image block radius,
 * moving image, and search region radius must be set.
 *
 * \ingroup Ultrasound
 */
template < class TFixedImage, class TMovingImage >
class ITK_EXPORT SearchRegionImageInitializer:
  public ::itk::ImageSource< 
    ::itk::Image< typename TMovingImage::RegionType,
      TMovingImage::ImageDimension > >
{
public:
  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TMovingImage::ImageDimension);

  /** Type of the fixed image. */
  typedef TFixedImage                         FixedImageType;
  typedef typename FixedImageType::RegionType FixedRegionType;

  /** Type of the radius used to characterized the fixed image block. */
  typedef typename FixedImageType::SizeType RadiusType;

  /** Type of the moving image. */
  typedef TMovingImage  MovingImageType;
  typedef typename MovingImageType::RegionType MovingRegionType;

  /** Type of the search region image. */
  typedef typename ::itk::Image< typename MovingImageType::RegionType, 
    ImageDimension > OutputImageType;
  typedef typename OutputImageType::RegionType OutputRegionType;

  /** Standard class typedefs. */
  typedef SearchRegionImageInitializer              Self;
  typedef ImageSource< OutputImageType >            Superclass;
  typedef SmartPointer< Self >                      Pointer;
  typedef SmartPointer< const Self >                ConstPointer;

  /** New macro for creation of through a Smart Pointer. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( SearchRegionImageInitializer, ImageSource );

  /** Set the fixed image. */
  void SetFixedImage( FixedImageType * fixedImage );
  const FixedImageType * GetFixedImage() const 
    { return this->m_FixedImage.GetPointer(); }

  /** Set the moving image. */
  void SetMovingImage( MovingImageType * movingImage );
  const MovingImageType * GetMovingImage() const
    { return this->m_MovingImage.GetPointer(); }

  /** Set the fixed block radius, i.e. the radius of the matching kernel from
   * the fixed image. */
  virtual void SetFixedBlockRadius( const RadiusType& radius )
    {
    m_FixedBlockRadius = radius;
    this->Modified();
    }
  itkGetConstMacro( FixedBlockRadius, RadiusType );
  /** Set the fixed block radius to the the same in all directions. */
  virtual void SetFixedBlockRadius( const unsigned int rad )
    {
    RadiusType radius;
    radius.Fill( rad );
    this->SetFixedBlockRadius( radius );
    }

  /** Set the search region image, i.e. the radius of the search region in the
   * moving image. */
  virtual void SetSearchRegionRadius( const RadiusType& radius )
    {
    m_SearchRegionRadius = radius;
    this->Modified();
    }
  itkGetConstMacro( SearchRegionRadius, RadiusType );
  /** Set the search region radius the same in all directions. */
  virtual void SetSearchRegionRadius( const unsigned int rad )
    {
    RadiusType radius;
    radius.Fill( rad );
    this->SetSearchRegionRadius( radius );
    }

  /** Set/Get the overlap between fixed image blocks.  This value should be
   * greater than zero and defaults to unity.  Values less than unity will have
   * the blocks overlapping, e.g. 0.5 will render 50% overlap.  Values greater
   * than unity will result in spacing between blocks. */
  itkSetMacro( Overlap, double );
  itkGetConstMacro( Overlap, double );

  // @todo make overlap into a SpacingType so the overlap can be specified
  // differently in every direction.

protected:
  SearchRegionImageInitializer();

  virtual void GenerateOutputInformation();

  virtual void BeforeThreadedGenerateData();

  virtual void ThreadedGenerateData( const OutputRegionType& outputRegion, int threadId );

  typename FixedImageType::Pointer  m_FixedImage;
  typename MovingImageType::Pointer m_MovingImage;

  RadiusType m_FixedBlockRadius;
  RadiusType m_SearchRegionRadius;

  double m_Overlap;

private:
  SearchRegionImageInitializer( const Self& );
  void operator=( const Self& );
};

} // end namespace BlockMatching
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlockMatchingSearchRegionImageInitializer.hxx"
#endif

#endif

