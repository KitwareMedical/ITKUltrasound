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
template< class TFixedImage, class TMovingImage, class TMetricImage >
class ITK_TEMPLATE_EXPORT MetricImageFilter :
  public ImageToImageFilter< TFixedImage, TMetricImage >
{
public:
  /** Standard class typedefs. */
  typedef MetricImageFilter                               Self;
  typedef ImageToImageFilter< TFixedImage, TMetricImage > Superclass;
  typedef SmartPointer<Self>                              Pointer;
  typedef SmartPointer<const Self>                        ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(MetricImageFilter, ImageToImageFilter);

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TFixedImage::ImageDimension);

  /**  Type of the Fixed image. */
  typedef          TFixedImage                FixedImageType;
  typedef typename FixedImageType::RegionType FixedImageRegionType;

  /**  Type of the Moving image. */
  typedef          TMovingImage                MovingImageType;
  typedef typename MovingImageType::RegionType MovingImageRegionType;

  /** Type of the Metric image. */
  typedef          TMetricImage                 MetricImageType;
  typedef typename MetricImageType::SpacingType MetricImageSpacingType;
  typedef typename MetricImageType::RegionType  MetricImageRegionType;

  /** Type of the block radius. */
  typedef typename FixedImageType::SizeType RadiusType;

  /** Set the fixed image. */
  void SetFixedImage( FixedImageType * fixedImage );

  /** Set the moving image. */
  void SetMovingImage( MovingImageType * movingImage );
  
  /** Set/Get the fixed image region.  This will defined the kernel in the fixed
   * image to be compared against the moving.  Every pixel in this region will
   * be compared against the moving image region to generate a pixel in the
   * output metric image. */
  void SetFixedImageRegion( const FixedImageRegionType& region );
  itkGetConstReferenceMacro( FixedImageRegion, FixedImageRegionType );

  /** Set/Get the moving image region.  This will defined the search region over
   * which the FixedImageRegion will be iteratively translated over to create
   * the metric image. */
  void SetMovingImageRegion( const MovingImageRegionType& region );
  itkGetConstReferenceMacro( MovingImageRegion, MovingImageRegionType );

  /** Set/Get the minimum size in the split requested region's split direction
   * (the last dimension by default).  Set this to a reasonably large number to
   * improve performance and prevent breakage. */
  itkSetMacro( MinimumSplitSize, int );
  itkGetConstMacro( MinimumSplitSize, int );

  /** Overload so we can decrease the number of threads according to the
   * MinimumSplitSize if need be. */
  virtual const int & GetNumberOfThreads();

protected:
  MetricImageFilter();
  virtual ~MetricImageFilter() {};

  virtual void GenerateOutputInformation();

  /** We set the fixed and moving image's requested region.  The fixed image's
   * requested region is equal to FixedImageRegion and the moving image's requested region is equal to the MovingImageRegion dilated by a radius related to the FixedImageRegion's size. */
  virtual void GenerateInputRequestedRegion();

  /** This filter produces the LargestPossibleRegion. */
  virtual void EnlargeOutputRequestedRegion( DataObject * data );

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

  int m_MinimumSplitSize;

  // Because it returns a reference?
  int m_SpecialThreadCount;

private:
  MetricImageFilter( const Self & ); // purposely not implemented
  void operator=( const Self& );     // purposely not implemented

};

} // end namespace BlockMatching
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlockMatchingMetricImageFilter.hxx"
#endif

#endif
