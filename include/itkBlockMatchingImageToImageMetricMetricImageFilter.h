#ifndef __itkBlockMatchingImageToImageMetricMetricImageFilter_h
#define __itkBlockMatchingImageToImageMetricMetricImageFilter_h

#include "itkBlockMatchingImageToImageMetricMetricImageFilter.h"

namespace itk
{
namespace BlockMatching
{

/** \class ImageToImageMetricMetricImageFilter
 *
 * \brief A BlockMatching::MetricImageFilter that uses any of ITK's
 * ImageToImageMetric's to calculate the metric.
 *
 * \sa MetricImageFilter
 * \sa ImageToImageMetric
 *
 * \ingroup Ultrasound
 */
template< class TFixedImage, class TMovingImage, class TMetricImage >
class ITK_EXPORT ImageToImageMetricMetricImageFilter :
  public MetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
{
public:
  /** Standard class typedefs. */
  typedef ImageToImageMetricMetricImageFilter			    Self;
  typedef MetricImageFilter< TFixedImage, TMovingImage, TMetricImage > Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageToImageMetricMetricImageFilter, MetricImageFilter);

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      Superclass::ImageDimension);

  /** Type of the Moving image. */
  typedef typename Superclass::MovingImageType                     MovingImageType;
  typedef typename MovingImageType::SpacingType     MovingImageSpacingType;
  typedef typename MovingImageType::RegionType      MovingImageRegionType;

  /** Type of the Metric image. */
  typedef typename Superclass::MetricImageType                     MetricImageType;
  typedef typename MetricImageType::SpacingType     MetricImageSpacingType;
  typedef typename MetricImageType::RegionType      MetricImageRegionType;

  /** Set the spacing in the output metric image.  This is optional.  If not
   * set, the spacing in the moving image will be used. */
  void SetMetricImageSpacing( const MetricImageSpacingType& spacing );
  itkGetConstReferenceMacro( MetricImageSpacing, MetricImageSpacingType );
protected:
  ImageToImageMetricMetricImageFilter();
  virtual ~ImageToImageMetricMetricImageFilter() {}

  virtual void GenerateOutputInformation();

  MetricImageSpacingType m_MetricImageSpacing;
  bool m_MetricImageSpacingDefined;

private:
  ImageToImageMetricMetricImageFilter( const Self& ); // purposely not implemented
  void operator=( const Self& ); // purposely not implemented

};

}
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlockMatchingImageToImageMetricMetricImageFilter.hxx"
#endif

#endif

