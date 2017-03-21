#ifndef __itkBlockMatchingMaximumPixelDisplacementCalculator_h
#define __itkBlockMatchingMaximumPixelDisplacementCalculator_h

#include "itkBlockMatchingMetricImageToDisplacementCalculator.h"

namespace itk
{
namespace BlockMatching
{

/** \class MaximumPixelDisplacementCalculator
 *
 * \brief The displacement is taken to be the pixel with the maximum metric
 * value.
 *
 * This is the simplest and fastest of the MetricImageToDisplacementCalculator's.
 *
 * \ingroup Ultrasound
 */
template < class TMetricImage, class TDisplacementImage >
class ITK_EXPORT MaximumPixelDisplacementCalculator:
  public itk::BlockMatching::MetricImageToDisplacementCalculator<
    TMetricImage, TDisplacementImage >
{
public:
  /** Standard class typedefs. */
  typedef MaximumPixelDisplacementCalculator     Self;
  typedef MetricImageToDisplacementCalculator< TMetricImage, TDisplacementImage >
    Superclass;

  typedef SmartPointer< Self >           Pointer;
  typedef SmartPointer< const Self >     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );
  
  /** Run-time type information (and related methods). */
  itkTypeMacro( MaximumPixelDisplacementCalculator, MetricImageToDisplacementCalculator );

  typedef typename Superclass::MetricImageType MetricImageType;
  typedef typename Superclass::MetricImagePointerType MetricImagePointerType;
  typedef typename MetricImageType::PixelType  PixelType;

  typedef typename Superclass::PointType PointType;
  
  typedef typename Superclass::IndexType IndexType;

  virtual void SetMetricImagePixel( const PointType & point, const IndexType& index, MetricImageType* image ); 

  virtual void Compute() {
    // We do this here instead of SetMetricImagePixel so it only has to be done
    // once.
    this->m_DisplacementImage->Modified();
  };

protected:
  MaximumPixelDisplacementCalculator() {}

private:
  MaximumPixelDisplacementCalculator( const Self & );
  void operator=( const Self & );
};

} // end namespace itk
} // end namespace BlockMatching

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlockMatchingMaximumPixelDisplacementCalculator.txx"
#endif

#endif
