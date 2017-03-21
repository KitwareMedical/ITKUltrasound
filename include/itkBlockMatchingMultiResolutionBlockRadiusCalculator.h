#ifndef __itkBlockMatchingMultiResolutionBlockRadiusCalculator_h
#define __itkBlockMatchingMultiResolutionBlockRadiusCalculator_h

#include "itkArray2D.h"
#include "itkObject.h"

namespace itk
{
namespace BlockMatching
{

/** \class MultiResolutionBlockRadiusCalculator
 *
 * \brief Base class calculate the fixed image matching kernel radius for each level in a
 * BlockMatching::MultiResolutionImageRegistrationMethod.
 *
 * This must be able to produce the fixed image block radius for every level of
 * the MultiResolutionPyramidImage filter.
 *
 * \ingroup Ultrasound
 */
template< class TFixedImage >
class ITK_EXPORT MultiResolutionBlockRadiusCalculator :
  public ::itk::Object
{
public:
  /** Standard class typedefs. */
  typedef MultiResolutionBlockRadiusCalculator  Self;
  typedef Object                                Superclass;

  typedef SmartPointer< Self >           Pointer;
  typedef SmartPointer< const Self >     ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( MultiResolutionBlockRadiusCalculator, Object );

  /** ImageDimension enumeration. */
  itkStaticConstMacro( ImageDimension, unsigned int,
                      TFixedImage::ImageDimension );

  /** Type of the fixed image type. */
  typedef TFixedImage                            FixedImageType;
  typedef typename FixedImageType::ConstPointer  FixedImageConstPointer;
  typedef typename FixedImageType::SizeType      RadiusType;

  /** ScheduleType typedef support. */
  typedef Array2D<unsigned int>  ScheduleType;

  /** SetPyramidSchedule() gets called with the pyramid schedule after the pyramid has
   * been generated.  This information is the available for child classes if
   * they choose to use it.  */
  virtual void SetPyramidSchedule( const ScheduleType& schedule )
    {
    m_PyramidSchedule = schedule;
    this->Modified();
    }

  /** Get the multi-resolution schedule. */
  itkGetConstReferenceMacro( PyramidSchedule, ScheduleType );

  /** Set/Get the Fixed image. The fixed image is set during the
   * BlockMatching::MultiResolutionImageRegistrationMethod and is available for
   * child classes if the choose to use it.  */
  itkSetConstObjectMacro( FixedImage, FixedImageType );
  itkGetConstObjectMacro( FixedImage, FixedImageType ); 

  /** This abstract method must be implemented by child classes.  Given the
   * current level in the image pyramid, return the fixed block radius. */
  virtual const RadiusType& Compute( unsigned long current_level ) = 0;

protected:
  MultiResolutionBlockRadiusCalculator() {}

  ScheduleType m_PyramidSchedule;

  FixedImageConstPointer m_FixedImage;

private:
  MultiResolutionBlockRadiusCalculator( const Self & );
  void operator=( const Self & );
};

}
}

#endif
