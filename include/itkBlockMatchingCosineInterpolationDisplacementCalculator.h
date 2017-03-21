#ifndef __itkBlockMatchingCosineInterpolationDisplacementCalculator_h
#define __itkBlockMatchingCosineInterpolationDisplacementCalculator_h

#include "itkBlockMatchingMetricImageToDisplacementCalculator.h"

namespace itk
{
namespace BlockMatching
{

/** \class CosineInterpolationDisplacementCalculator
 *
 * \brief The displacement around the maximum pixel is interpolated by a cosine
 * in each direction and the peak of the cosine is used.
 *
 * Cespedes et. al. Methods for estimation of subsample time delays of digitized
 * echo signals.  Ultrasonic Imaging 17. 142-171.  1995.
 *
 * \ingroup Ultrasound
 */
template < class TMetricImage, class TDisplacementImage, class TCoordRep=double >
class ITK_EXPORT CosineInterpolationDisplacementCalculator:
  public itk::BlockMatching::MetricImageToDisplacementCalculator<
    TMetricImage, TDisplacementImage >
{
public:
  /** Standard class typedefs. */
  typedef CosineInterpolationDisplacementCalculator     Self;
  typedef MetricImageToDisplacementCalculator< TMetricImage, TDisplacementImage >
    Superclass;

  typedef SmartPointer< Self >           Pointer;
  typedef SmartPointer< const Self >     ConstPointer;

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      Superclass::ImageDimension);

  /** Method for creation through the object factory. */
  itkNewMacro( Self );
  
  /** Run-time type information (and related methods). */
  itkTypeMacro( CosineInterpolationDisplacementCalculator, MetricImageToDisplacementCalculator );

  typedef typename Superclass::MetricImageType MetricImageType;
  typedef typename Superclass::MetricImagePointerType MetricImagePointerType;
  typedef typename MetricImageType::PixelType  PixelType;
  typedef typename MetricImageType::SpacingType SpacingType;

  typedef typename Superclass::PointType PointType;
  
  typedef typename Superclass::IndexType IndexType;

  virtual void SetMetricImagePixel( const PointType & point, const IndexType& index, MetricImageType* image ); 

  virtual void Compute() {
    // We do this here instead of SetMetricImagePixel so it only has to be done
    // once.
    this->m_DisplacementImage->Modified();
  };

protected:
  CosineInterpolationDisplacementCalculator();

private:
  CosineInterpolationDisplacementCalculator( const Self & );
  void operator=( const Self & );
};

} // end namespace itk
} // end namespace BlockMatching

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlockMatchingCosineInterpolationDisplacementCalculator.txx"
#endif

#endif
