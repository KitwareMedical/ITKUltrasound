#ifndef __itkBlockMatchingMultiResolutionFixedBlockRadiusCalculator_h
#define __itkBlockMatchingMultiResolutionFixedBlockRadiusCalculator_h

#include "itkBlockMatchingMultiResolutionBlockRadiusCalculator.h"

namespace itk
{
namespace BlockMatching
{
/** \class MultiResolutionFixedBlockRadiusCalculator
 *
 * \brief A fixed radius is used for every level.
 *
 * This class generates the fixed image matching kernel radius in a
 * BlockMatching::MultiResolutionImageRegistrationMethod.  A fixed block radius
 * is used, which means the size of the block in physical coordinates then
 * scales with the pyramid schedule.
 *
 * \sa MultiResolutionBlockRadiusCalculator
 * */
template< class TFixedImage >
class ITK_EXPORT MultiResolutionFixedBlockRadiusCalculator :
  public  MultiResolutionBlockRadiusCalculator< TFixedImage >
{
public:
  /** Standard class typedefs. */
  typedef MultiResolutionFixedBlockRadiusCalculator           Self;
  typedef MultiResolutionBlockRadiusCalculator< TFixedImage > Superclass;

  typedef SmartPointer< Self >            Pointer;
  typedef SmartPointer< const Self >      ConstPointer;

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( MultiResolutionFixedBlockRadiusCalculator, MultiResolutionBlockRadiusCalculator );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** ImageDimension enumeration. */
  itkStaticConstMacro( ImageDimension, unsigned int,
                      Superclass::ImageDimension );

  /** Type of the fixed image. */
  typedef typename Superclass::FixedImageType  FixedImageType;
  typedef typename Superclass::RadiusType      RadiusType;

  /** Set the radius used for all pyramid levels. */
  void SetRadius( const RadiusType& radius )
    { m_Radius = radius; }
  itkGetConstReferenceMacro( Radius, RadiusType );

  virtual const RadiusType& Compute( unsigned long level ) 
    {
    return m_Radius;
    }

protected:
  MultiResolutionFixedBlockRadiusCalculator()
   {
   m_Radius.Fill( 1 );
   }

  RadiusType m_Radius;

private:
  MultiResolutionFixedBlockRadiusCalculator( const Self& );
  void operator=( const Self& );
};

} // end namespace itk
} // end namespace BlockMatching

#endif 
