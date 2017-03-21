#ifndef __itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter_h
#define __itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter_h

#include "itkBlockMatchingNormalizedCrossCorrelationMetricImageFilter.h"

namespace itk
{
namespace BlockMatching
{

/** \class NormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter
 *
 * \brief Create an image of the the normalized cross correlation with a kernel
 * calculated with a neighborhood iterator.
 *
 * \sa NormalizedCrossCorrelationMetricImageFilter
 *
 * \ingroup Ultrasound
 */
template< class TFixedImage, class TMovingImage, class TMetricImage >
class ITK_EXPORT NormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter :
  public NormalizedCrossCorrelationMetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
{
public:
  /** Standard class typedefs. */
  typedef NormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter                        Self;
  typedef NormalizedCrossCorrelationMetricImageFilter< TFixedImage, TMovingImage, TMetricImage > Superclass;
  typedef SmartPointer<Self>                                                                     Pointer;
  typedef SmartPointer<const Self>                                                               ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(NormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter,
               NormalizedCrossCorrelationMetricImageFilter);

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TFixedImage::ImageDimension);

  /** Type of the fixed image. */
  typedef typename Superclass::FixedImageType   FixedImageType;
  typedef typename FixedImageType::ConstPointer FixedImageConstPointerType;

  /** Type of the moving image. */
  typedef typename Superclass::MovingImageType   MovingImageType;
  typedef typename MovingImageType::ConstPointer MovingImageConstPointerType;


  /** Type of the metric image. */
  typedef typename Superclass::MetricImageType   MetricImageType;
  typedef typename MetricImageType::Pointer      MetricImagePointerType;
  typedef typename MetricImageType::ConstPointer MetricImageConstPointerType;
  typedef typename MetricImageType::RegionType   MetricImageRegionType;
  typedef typename MetricImageType::PixelType    MetricImagePixelType;

protected:
  NormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter(){}

  virtual void BeforeThreadedGenerateData();

  virtual void ThreadedGenerateData( const MetricImageRegionType& outputRegion, int threadId );

private:
  NormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter( const Self& ); // purposely not implemented
  void operator=( const Self& ); // purposely not implemented
};

} // end namespace BlockMatching
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter.hxx"
#endif

#endif
