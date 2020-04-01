/*=========================================================================
 *
 *  Copyright NumFOCUS
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkBlockMatchingStrainWindowBlockAffineTransformCommand_hxx
#define itkBlockMatchingStrainWindowBlockAffineTransformCommand_hxx

#include "itkBlockMatchingStrainWindowBlockAffineTransformCommand.h"

namespace itk
{
namespace BlockMatching
{

template<typename TStrainWindowDisplacemenCalculator,
         typename TBlockAffineTransformMetricImageFilter,
         typename TStrainImageFilter>
StrainWindowBlockAffineTransformCommand<TStrainWindowDisplacemenCalculator,
                                        TBlockAffineTransformMetricImageFilter,
                                        TStrainImageFilter>
::StrainWindowBlockAffineTransformCommand():
  m_UseStrainWindowStrain( true )
{
  m_BlockAffineTransformMetricImageFilter = nullptr;

  m_StrainImageFilter = StrainImageFilterType::New();
  m_LeastSquaresFilter = LeastSquaresFilterType::New();
  typename LeastSquaresFilterType::RadiusType radius; // @todo make this settable.
  radius.Fill( 3 );
  m_LeastSquaresFilter->SetRadius( radius );
  m_StrainImageFilter->SetGradientFilter( m_LeastSquaresFilter );
}

template<typename TStrainWindowDisplacemenCalculator,
         typename TBlockAffineTransformMetricImageFilter,
         typename TStrainImageFilter>
void
StrainWindowBlockAffineTransformCommand<TStrainWindowDisplacemenCalculator,
                                        TBlockAffineTransformMetricImageFilter,
                                        TStrainImageFilter>
::Execute( const itk::Object * object, const itk::EventObject & event )
{
  if( !(itk::EndEvent().CheckEvent( &event )) )
    {
    return;
    }

  if( m_BlockAffineTransformMetricImageFilter.GetPointer() == nullptr )
    {
    itkExceptionMacro(<< "The BlockAffineTransformMetricImageFilter has not been set." );
    }

  StrainWindowDisplacementCalculatorType * strainWindower = const_cast< StrainWindowDisplacementCalculatorType* >(
    dynamic_cast< const StrainWindowDisplacementCalculatorType * >( object ) );
  if( !strainWindower )
    {
    itkExceptionMacro(<<"Could not downcast to a StrainWindowDisplacementCalculator.");
    }

  if( m_UseStrainWindowStrain )
    {
    m_BlockAffineTransformMetricImageFilter->SetStrainImage( strainWindower->GetStrainImageFilter()->GetOutput() );
    }
  else
    {
    m_StrainImageFilter->SetInput( strainWindower->GetStrainImageFilter()->GetInput() );
    m_StrainImageFilter->Update();
    m_BlockAffineTransformMetricImageFilter->SetStrainImage( m_StrainImageFilter->GetOutput() );
    }
}

} // end namespace BlockMatching
} // end namespace itk


#endif
