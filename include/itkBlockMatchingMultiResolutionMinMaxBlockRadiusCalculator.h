#ifndef __itkBlockMatchingMultiResolutionMinMaxBlockRadiusCalculator_h
#define __itkBlockMatchingMultiResolutionMinMaxBlockRadiusCalculator_h

#include "itkBlockMatchingMultiResolutionBlockRadiusCalculator.h"

namespace itk
{
namespace BlockMatching
{
/** \class MultiResolutionMinMaxBlockRadiusCalculator
 *
 * \brief The radius varies linearly from a minimum at the lowest level to a
 * maximum at the highest level.
 *
 * \todo rename this from 'MinMax' to 'TopBottom'.
 *
 * \sa MultiResolutionBlockRadiusCalculator
 */
template< class TFixedImage >
class ITK_EXPORT MultiResolutionMinMaxBlockRadiusCalculator :
  public MultiResolutionBlockRadiusCalculator< TFixedImage >
{
public:
  /** Standard class typedefs. */
  typedef MultiResolutionMinMaxBlockRadiusCalculator          Self;
  typedef MultiResolutionBlockRadiusCalculator< TFixedImage > Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( MultiResolutionMinMaxBlockRadiusCalculator, MultiResolutionBlockRadiusCalculator );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** ImageDimension enumeration. */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       Superclass::ImageDimension );

  /** Type of the fixed image. */
  typedef typename Superclass::FixedImageType FixedImageType;
  typedef typename Superclass::RadiusType     RadiusType;

  /** Set the radius used for at the bottom level. */
  virtual void SetMinRadius( const RadiusType& radius )
  {
    m_MinRadius = radius;
  }
  itkGetConstReferenceMacro( MinRadius, RadiusType );

  /** Set the radius used at the top level. */
  virtual void SetMaxRadius( const RadiusType& radius )
  {
    m_MaxRadius = radius;
  }
  itkGetConstReferenceMacro( Radius, RadiusType );

  virtual const RadiusType& Compute( unsigned long level )
    {
    double slope;
    double distance = static_cast< double >( this->m_PyramidSchedule.rows() - 1 );
    for( unsigned int i = 0; i < this->m_PyramidSchedule.cols(); ++i )
      {
      slope = ( static_cast< double >( m_MinRadius[i] ) - static_cast< double >( m_MaxRadius[i] ) ) / distance;
      m_Radius[i] = static_cast< typename RadiusType::SizeValueType >( slope * level + m_MaxRadius[i] );
      }
    return m_Radius;
    }

protected:
  MultiResolutionMinMaxBlockRadiusCalculator()
   {
   m_MinRadius.Fill( 1 );
   m_MaxRadius.Fill( 1 );
   }

  RadiusType m_MinRadius;
  RadiusType m_MaxRadius;

  RadiusType m_Radius;

private:
  MultiResolutionMinMaxBlockRadiusCalculator( const Self & );
  void operator=( const Self& );

};

} // end namespace itk
} // end namespace BlockMatching


#endif
