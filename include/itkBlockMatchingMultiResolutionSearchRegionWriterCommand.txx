#ifndef __itkBlockMatchingMultiResolutionSearchRegionWriterCommand_txx
#define __itkBlockMatchingMultiResolutionSearchRegionWriterCommand_txx

#include <fstream>
#include <iostream>
#include <sstream>

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

namespace itk
{
namespace BlockMatching
{

template< class TMultiResolutionMethod >
MultiResolutionSearchRegionWriterCommand< TMultiResolutionMethod >
::MultiResolutionSearchRegionWriterCommand()
{
  m_SearchRegionImageComponent  = SearchRegionImageComponentType::New();
  m_SearchRegionComponentWriter = SearchRegionComponentWriterType::New();
  m_SearchRegionComponentWriter->SetInput( m_SearchRegionImageComponent );
}

template< class TMultiResolutionMethod >
MultiResolutionSearchRegionWriterCommand< TMultiResolutionMethod >
::~MultiResolutionSearchRegionWriterCommand()
{
}

template< class TMultiResolutionMethod >
void
MultiResolutionSearchRegionWriterCommand< TMultiResolutionMethod >
::Execute( const itk::Object * object , const itk::EventObject & event )
{
  Superclass::Execute( object, event );

  const unsigned long level = this->m_MultiResolutionMethod->GetCurrentLevel();

  this->m_SearchRegionImageSource->UpdateLargestPossibleRegion();
  typename SearchRegionImageType::ConstPointer searchRegionImage = this->m_SearchRegionImageSource->GetOutput();
  m_SearchRegionImageComponent->CopyInformation( searchRegionImage );
  m_SearchRegionImageComponent->SetRegions( searchRegionImage->GetBufferedRegion() );
  m_SearchRegionImageComponent->Allocate();

  std::cout << "Writing SearchRegionImage images..." << std::endl;


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


}
}

#endif
