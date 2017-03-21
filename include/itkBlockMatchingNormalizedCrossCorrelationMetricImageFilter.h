#ifndef __itkBlockMatchingNormalizedCrossCorrelationMetricImageFilter_h
#define __itkBlockMatchingNormalizedCrossCorrelationMetricImageFilter_h

#include "itkBoxMeanImageFilter.h"
#include "itkBoxSigmaSqrtNMinusOneImageFilter.h"

#include "itkBlockMatchingMetricImageFilter.h"

namespace itk
{
namespace BlockMatching
{

/** \class NormalizedCrossCorrelationMetricImageFilter
 *
 * \brief A MetricImageFilter that calculates a metric image of
 * normalized cross correlation.
 *
 * r_{xy} = \frac{1}{n-1} \sum_{x,y}\frac{(x) - \overline{f})(y -
 * \overline{y})}{\sigma_f \sigma_y}
 *
 * This is an abstract base class that does the mean and standard deviation
 * calculation.  The cross correlation is left to inherited classes.
 *
 * \sa MetricImageFilter
 *
 * \ingroup Ultrasound
 */
template< class TFixedImage, class TMovingImage, class TMetricImage >
class ITK_EXPORT NormalizedCrossCorrelationMetricImageFilter :
  public MetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
{
public:
  /** Standard class typedefs. */
  typedef NormalizedCrossCorrelationMetricImageFilter                  Self;
  typedef MetricImageFilter< TFixedImage, TMovingImage, TMetricImage > Superclass;
  typedef SmartPointer<Self>                                           Pointer;
  typedef SmartPointer<const Self>                                     ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(NormalizedCrossCorrelationMetricImageFilter, MetricImageFilter);

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      Superclass::ImageDimension);

  /** Type of the fixed image. */
  typedef typename Superclass::FixedImageType   FixedImageType;
  typedef typename FixedImageType::ConstPointer FixedImageConstPointerType;

  /** Type of the moving image. */
  typedef typename Superclass::MovingImageType   MovingImageType;
  typedef typename MovingImageType::RegionType   MovingImageRegionType;
  typedef typename MovingImageType::ConstPointer MovingImageConstPointerType;

  /** Type of the metric image. */
  typedef typename Superclass::MetricImageType MetricImageType;
  typedef typename MetricImageType::Pointer    MetricImagePointerType;
  typedef typename MetricImageType::PixelType  MetricImagePixelType;
protected:
  NormalizedCrossCorrelationMetricImageFilter();

  /** The mean and pseudo-standarddeviation images are stored in the outputs so
    they fix in with the pipline architecture. */
  virtual void GenerateOutputInformation();

  /** All outputs generate the largest possible region. */
  virtual void EnlargeOutputRequestedRegion( DataObject * data );

  /** Don't let the default mess with our output requested regions. */
  virtual void GenerateOutputRequestedRegion( DataObject * data ) {};

  /** Generates helper images for the calculation.  These are only needed for
   * internal calculation, but they are put on the
   * outputs for use by subclasses and to use the pipeline memory management
   * system.
   *
   * Calculates an image of means for each block neighborhood.  Also calculates
   * an image of standard deviations (times sqrt(N-1)) for each block neighborhood. */
  virtual void GenerateHelperImages();

  typedef typename ::itk::BoxMeanImageFilter< MovingImageType, MetricImageType >               BoxMeanFilterType;
  typedef typename ::itk::BoxSigmaSqrtNMinusOneImageFilter< MovingImageType, MetricImageType > BoxPseudoSigmaFilterType;

  typename BoxMeanFilterType::Pointer   m_BoxMeanFilter;
  typename BoxPseudoSigmaFilterType::Pointer m_BoxPseudoSigmaFilter;

private:
  NormalizedCrossCorrelationMetricImageFilter( const Self& ); // purposely not implemented
  void operator=( const Self& ); // purposely not implemented

  typedef typename ::itk::ConstantBoundaryCondition< MetricImageType > BoundaryConditionType;
  BoundaryConditionType m_BoundaryCondition;
};

} // end namespace BlockMatching
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlockMatchingNormalizedCrossCorrelationMetricImageFilter.txx"
#endif

#endif

