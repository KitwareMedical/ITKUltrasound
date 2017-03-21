#ifndef __itkBlockMatchingParabolicInterpolationDisplacementCalculator_h
#define __itkBlockMatchingParabolicInterpolationDisplacementCalculator_h

#include "itkBlockMatchingMetricImageToDisplacementCalculator.h"

namespace itk
{
namespace BlockMatching
{

/** \class ParabolicInterpolationDisplacementCalculator
 *
 * \brief The displacement around the maximum pixel is interpolated by a
 * parabola in each direction and the peak of the parabola is used.
 *
 * Cespedes et. al. Methods for estimation of subsample time delays of digitized
 * echo signals.  Ultrasonic Imaging 17. 142-171.  1995.  
 *
 * \ingroup Ultrasound
 */
template < class TMetricImage, class TDisplacementImage, class TCoordRep=double > 
class ITK_EXPORT ParabolicInterpolationDisplacementCalculator:
  public MetricImageToDisplacementCalculator< TMetricImage, TDisplacementImage >
{
public:
  /** Standard class typedefs. */
  typedef ParabolicInterpolationDisplacementCalculator     Self; typedef
    MetricImageToDisplacementCalculator< TMetricImage, TDisplacementImage >
    Superclass;

  typedef SmartPointer< Self >           Pointer;
  typedef SmartPointer< const Self >     ConstPointer;

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      Superclass::ImageDimension);

  /** Method for creation through the object factory. */
  itkNewMacro( Self );
  
  /** Run-time type information (and related methods). */
  itkTypeMacro( ParabolicInterpolationDisplacementCalculator, MetricImageToDisplacementCalculator );

  typedef typename Superclass::MetricImageType        MetricImageType;
  typedef typename Superclass::MetricImagePointerType MetricImagePointerType;
  typedef typename MetricImageType::PixelType         PixelType;
  typedef typename MetricImageType::SpacingType       SpacingType;
  typedef typename MetricImageType::RegionType        RegionType;

  typedef typename Superclass::MetricImageImageType MetricImageImageType;
  typedef typename Superclass::MetricImageImagePointerType
    MetricImageImagePointerType;
  typedef ImageRegionIterator< MetricImageImageType >
    MetricImageImageIteratorType;

  typedef typename Superclass::CenterPointsImageType CenterPointsImageType;

  typedef typename Superclass::DisplacementImageType DisplacementImageType;

  typedef typename Superclass::PointType PointType;
  typedef typename Superclass::IndexType IndexType;

  virtual void Compute();

protected:
  ParabolicInterpolationDisplacementCalculator();

  typedef typename Superclass::ThreadFunctor  ThreadFunctor;
  typedef typename Superclass::ThreadStruct   ThreadStruct;

  /** Use a parabolic fit to find the subsample peak. */
  class ParabolicInterpolationThreadFunctor : public ThreadFunctor
    {
  public:
    virtual ITK_THREAD_RETURN_TYPE operator()( Superclass *superclass,
                                               RegionType & region, int threadId );
    };
  ParabolicInterpolationThreadFunctor m_ParabolicInterpolationThreadFunctor;

private:
  ParabolicInterpolationDisplacementCalculator( const Self & );
  void operator=( const Self & );
};

} // end namespace itk
} // end namespace BlockMatching

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlockMatchingParabolicInterpolationDisplacementCalculator.txx"
#endif

#endif
