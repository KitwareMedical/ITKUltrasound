#ifndef __itkBlockMatchingMultiResolutionIterationDisplacementCalculatorCommand_txx
#define __itkBlockMatchingMultiResolutionIterationDisplacementCalculatorCommand_txx

#include "itkBlockMatchingMultiResolutionIterationDisplacementCalculatorCommand.h"

namespace itk
{
namespace BlockMatching
{

template< class TMultiResolutionMethod >
void
MultiResolutionIterationDisplacementCalculatorCommand< TMultiResolutionMethod >
::Execute( const itk::Object * object , const itk::EventObject & event )
{
  Superclass::Execute( object, event );

  const unsigned long level = this->m_MultiResolutionMethod->GetCurrentLevel();

  if( m_Level0ToNMinus1DisplacementCalculator.GetPointer() == NULL )
    {
    itkExceptionMacro( << "Level0ToNMinus1ToNMinus1DisplacementCalculator is not present." );
    }

  if( m_LevelNDisplacementCalculator.GetPointer() == NULL )
    {
    itkExceptionMacro( << "LevelNDisplacementCalculator is not present." );
    }

  if( m_Regularizer.GetPointer() == NULL )
    {
    if( level == this->m_MultiResolutionMethod->GetNumberOfLevels() - 1 )
      {
      this->m_MultiResolutionMethod->GetImageRegistrationMethod()->SetMetricImageToDisplacementCalculator(
        this->m_LevelNDisplacementCalculator );
      }
    else
      {
      this->m_MultiResolutionMethod->GetImageRegistrationMethod()->SetMetricImageToDisplacementCalculator(
        this->m_Level0ToNMinus1DisplacementCalculator );
      }
    }
  else
    {
    if( level == this->m_MultiResolutionMethod->GetNumberOfLevels() - 1 )
      {
      m_Regularizer->SetDisplacementCalculator( this->m_LevelNDisplacementCalculator );
      m_Regularizer->SetMaximumIterations( m_LevelNRegularizerIterations );
      }
    else
      {
      m_Regularizer->SetDisplacementCalculator( this->m_Level0ToNMinus1DisplacementCalculator );
      m_Regularizer->SetMaximumIterations( m_Level0ToNMinus1RegularizerIterations );
      }
    //this->m_MultiResolutionMethod->GetImageRegistrationMethod()->SetMetricImageToDisplacementCalculator(
      //this->m_Regularizer );
    }
}

}
}

#endif
