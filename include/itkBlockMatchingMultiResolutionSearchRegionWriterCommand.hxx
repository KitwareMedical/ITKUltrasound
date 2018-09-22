/*=========================================================================
 *
 *  Copyright Insight Software Consortium
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
#ifndef itkBlockMatchingMultiResolutionSearchRegionWriterCommand_hxx
#define itkBlockMatchingMultiResolutionSearchRegionWriterCommand_hxx

#include "itkBlockMatchingMultiResolutionSearchRegionWriterCommand.h"

#include <fstream>
#include <iostream>
#include <sstream>

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

namespace itk
{
namespace BlockMatching
{

template< typename TMultiResolutionMethod >
MultiResolutionSearchRegionWriterCommand< TMultiResolutionMethod >
::MultiResolutionSearchRegionWriterCommand()
{
  m_SearchRegionImageComponent  = SearchRegionImageComponentType::New();
  m_SearchRegionComponentWriter = SearchRegionComponentWriterType::New();
  m_SearchRegionComponentWriter->SetInput( m_SearchRegionImageComponent );
}


template< typename TMultiResolutionMethod >
MultiResolutionSearchRegionWriterCommand< TMultiResolutionMethod >
::~MultiResolutionSearchRegionWriterCommand()
{
}


template< typename TMultiResolutionMethod >
void
MultiResolutionSearchRegionWriterCommand< TMultiResolutionMethod >
::Execute( const itk::Object * object , const itk::EventObject & event )
{
  Superclass::Execute( object, event );

  const unsigned long level = this->m_MultiResolutionMethod->GetCurrentLevel();

  typename SearchRegionImageType::ConstPointer searchRegionImage = this->m_SearchRegionImageSource->GetOutput();
  m_SearchRegionImageComponent->CopyInformation( searchRegionImage );
  m_SearchRegionImageComponent->SetRegions( searchRegionImage->GetBufferedRegion() );
  m_SearchRegionImageComponent->Allocate();

  if( m_OutputFilePrefix.size() == 0 )
    {
    itkExceptionMacro( << "OutputFilePrefix is not set." );
    }

  std::ostringstream ostr;
  itk::ImageRegionConstIterator< SearchRegionImageType > regionIt( searchRegionImage,
    searchRegionImage->GetBufferedRegion() );
  itk::ImageRegionIterator< SearchRegionImageComponentType > compIt( m_SearchRegionImageComponent,
    m_SearchRegionImageComponent->GetBufferedRegion() );
  for( unsigned int dim = 0; dim < 2; ++dim )
    {
    for( regionIt.GoToBegin(), compIt.GoToBegin();
         !regionIt.IsAtEnd();
         ++regionIt, ++compIt )
      {
      compIt.Set( regionIt.Get().GetSize()[dim] );
      }

    ostr.str( "" );
    ostr << m_OutputFilePrefix << "_Level_" << level << "_SearchRegionImageComponent" << dim << ".mha";
    m_SearchRegionComponentWriter->SetFileName( ostr.str() );
    m_SearchRegionComponentWriter->Update();
    }
}

} // end namespace BlockMatching
} // end namespace itk

#endif
