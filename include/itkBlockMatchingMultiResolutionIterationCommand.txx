#ifndef __itkBlockMatchingMultiResolutionIterationCommand_txx
#define __itkBlockMatchingMultiResolutionIterationCommand_txx

#include "itkBlockMatchingMultiResolutionIterationCommand.h"

namespace itk
{
namespace BlockMatching
{

template< class TMultiResolutionMethod >
void
MultiResolutionIterationCommand< TMultiResolutionMethod >
::Execute( const itk::Object * object , const itk::EventObject & event )
{
  if( !(IterationEvent().CheckEvent( &event )) )
    {
    return;
    }

  if( m_MultiResolutionMethod.GetPointer() == NULL )
    {
    itkExceptionMacro( << "The associated MultiResolutionMethod must be set." );
    }

  m_FixedImagePyramid     = m_MultiResolutionMethod->GetFixedImagePyramid();
  m_MovingImagePyramid    = m_MultiResolutionMethod->GetMovingImagePyramid();

  m_BlockRadiusCalculator = m_MultiResolutionMethod->GetBlockRadiusCalculator();

  m_SearchRegionImageSource = m_MultiResolutionMethod->GetSearchRegionImageSource();
}

} // end namespace BlockMatching
} // end namespace itk

#endif
