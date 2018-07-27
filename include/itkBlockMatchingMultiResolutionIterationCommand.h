#ifndef __itkBlockMatchingMultiResolutionIterationCommand_h
#define __itkBlockMatchingMultiResolutionIterationCommand_h

#include "itkCommand.h"

#include "itkBlockMatchingMultiResolutionImageRegistrationMethod.h"

namespace itk
{
namespace BlockMatching
{

/** This is a base class for classes that want to observe/adjust a
 * BlockMatching::MultiResolutionImageRegistrationMethod at every iteration.
 * */
template< class TMultiResolutionMethod >
class MultiResolutionIterationCommand :
  public itk::Command
{
public:
  typedef MultiResolutionIterationCommand  Self;
  typedef Command                          Superclass;
  typedef SmartPointer< Self >             Pointer;

  virtual void Execute(itk::Object *caller, const itk::EventObject & event)
    {
    Execute( (const itk::Object *)caller, event);
    }

  virtual void Execute(const itk::Object * object, const itk::EventObject & event);

  typedef TMultiResolutionMethod MultiResolutionMethodType;
  typedef typename MultiResolutionMethodType::Pointer MultiResolutionMethodPointer;

  void SetMultiResolutionMethod( MultiResolutionMethodPointer& method )
    {
    m_MultiResolutionMethod = method;
    }

  typedef typename MultiResolutionMethodType::FixedImagePyramidType
    FixedImagePyramidType;
  typedef typename FixedImagePyramidType::Pointer FixedImagePyramidPointer;
  typedef typename MultiResolutionMethodType::MovingImagePyramidType
    MovingImagePyramidType;
  typedef typename MovingImagePyramidType::Pointer MovingImagePyramidPointer;

  typedef typename MultiResolutionMethodType::BlockRadiusCalculatorType
    BlockRadiusCalculatorType;
  typedef typename BlockRadiusCalculatorType::Pointer BlockRadiusCalculatorPointer;

  typedef typename MultiResolutionMethodType::SearchRegionImageSourceType
    SearchRegionImageSourceType;
  typedef typename SearchRegionImageSourceType::Pointer SearchRegionImageSourcePointer;

protected:
  MultiResolutionIterationCommand()
    {
    m_MultiResolutionMethod = NULL;
    };
  virtual ~MultiResolutionIterationCommand() {};

  MultiResolutionMethodPointer m_MultiResolutionMethod;

  FixedImagePyramidPointer  m_FixedImagePyramid;
  MovingImagePyramidPointer m_MovingImagePyramid;

  BlockRadiusCalculatorPointer m_BlockRadiusCalculator;

  SearchRegionImageSourcePointer m_SearchRegionImageSource;

private:
  MultiResolutionIterationCommand( const Self& );
  void operator=( const Self & );
};

} // end namespace BlockMatching
} // end namespace itk

#include "itkBlockMatchingMultiResolutionIterationCommand.txx"

#endif
