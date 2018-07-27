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
#ifndef itkBlockMatchingMultiResolutionIterationObserver_hxx
#define itkBlockMatchingMultiResolutionIterationObserver_hxx

#include "itkBlockMatchingMultiResolutionIterationObserver.h"

#include <iostream>
#include <sstream>

namespace itk
{
namespace BlockMatching
{

template< typename TMultiResolutionMethod >
MultiResolutionIterationObserver< TMultiResolutionMethod >
::MultiResolutionIterationObserver()
{
  m_FixedImageWriter             = FixedImageWriterType::New();
  m_MovingImageWriter            = MovingImageWriterType::New();

  m_DisplacementWriter           = DisplacementWriterType::New();
  m_DisplacementComponentsFilter = DisplacementComponentsFilterType::New();
  m_DisplacementComponentsWriter = DisplacementComponentsWriterType::New();

  m_StrainFilter                 = StrainFilterType::New();
  m_StrainWriter                 = StrainWriterType::New();
  m_StrainWriter->SetInput( m_StrainFilter->GetOutput() );
  m_StrainComponentsFilter       = StrainComponentsFilterType::New();
  m_StrainComponentsFilter->SetInput( m_StrainFilter->GetOutput() );
  m_StrainComponentsWriter       = StrainComponentsWriterType::New();
}


template< typename TMultiResolutionMethod >
MultiResolutionIterationObserver< TMultiResolutionMethod >
::~MultiResolutionIterationObserver()
{
  if( m_CSVFile.is_open() )
    {
    m_CSVFile.flush();
    m_CSVFile.close();
    }
}


template< typename TMultiResolutionMethod >
void
MultiResolutionIterationObserver< TMultiResolutionMethod >
::Execute( const itk::Object * object , const itk::EventObject & event )
{
  Superclass::Execute( object, event );

  const unsigned long level = this->m_MultiResolutionMethod->GetCurrentLevel();
  std::cout << "Current Level: " << level + 1;
  std::cout << " / "             << this->m_MultiResolutionMethod->GetNumberOfLevels() << std::endl;

  if( !m_CSVFile.is_open() )
    {
    m_CSVFile.open( ( m_OutputFilePrefix + "_BlockRadius.csv" ).c_str() );
    if( !m_CSVFile.is_open() )
      throw std::runtime_error( "Could not open multi-level status file." );
    m_CSVFile << "Level, Block Radius" << std::endl;
    }
  m_CSVFile << level << ",  ";
  m_CSVFile << this->m_BlockRadiusCalculator->Compute( level ) << std::endl;

  // Skip the base level where the multilevel pyramid filter is not used.
  if( level < this->m_MultiResolutionMethod->GetNumberOfLevels() )
    {
    std::cout << "Writing fixed image..." << std::endl;
    std::ostringstream ostr;
    ostr << m_OutputFilePrefix << "_Level_" << level << "_FixedImage.mha";
    this->m_FixedImageWriter->SetInput( this->m_FixedImagePyramid->GetOutput( level ) );
    this->m_FixedImageWriter->SetFileName( ostr.str() );
    this->m_FixedImageWriter->Update();

    std::cout << "Writing moving image..." << std::endl;
    ostr.str( "" );
    ostr << m_OutputFilePrefix << "_Level_" << level << "_MovingImage.mha";
    this->m_MovingImageWriter->SetInput( this->m_MovingImagePyramid->GetOutput( level ) );
    this->m_MovingImageWriter->SetFileName( ostr.str() );
    this->m_MovingImageWriter->Update();

    if( level > 0 )
      {
      std::cout << "Writing displacement image..." << std::endl;
      ostr.str( "" );
      ostr << m_OutputFilePrefix << "_Level_" << level << "_PreviousDisplacements.mha";
      m_DisplacementWriter->SetInput( this->m_SearchRegionImageSource->GetPreviousDisplacements() );
      m_DisplacementWriter->SetFileName( ostr.str() );
      m_DisplacementWriter->Update();
      m_DisplacementComponentsFilter->SetInput( this->m_SearchRegionImageSource->GetPreviousDisplacements() );
      ostr.str( "" );
      ostr << m_OutputFilePrefix << "_Level_" << level << "_PreviousDisplacementComponent0.mha";
      m_DisplacementComponentsWriter->SetFileName( ostr.str() );
      m_DisplacementComponentsWriter->SetInput( m_DisplacementComponentsFilter->GetOutput( 0 ) );
      m_DisplacementComponentsWriter->Update();
      ostr.str( "" );
      ostr << m_OutputFilePrefix << "_Level_" << level << "_PreviousDisplacementComponent1.mha";
      m_DisplacementComponentsWriter->SetFileName( ostr.str() );
      m_DisplacementComponentsWriter->SetInput( m_DisplacementComponentsFilter->GetOutput( 1 ) );
      m_DisplacementComponentsWriter->Update();

      std::cout << "Writing strain image..." << std::endl;
      m_StrainFilter->SetInput( this->m_SearchRegionImageSource->GetPreviousDisplacements() );
      ostr.str( "" );
      ostr << m_OutputFilePrefix << "_Level_" << level << "_PreviousStrains.mha";
      m_StrainWriter->SetFileName( ostr.str() );
      m_StrainWriter->Update();
      for( unsigned int i = 0; i < 3; i++ )
        {
        m_StrainComponentsWriter->SetInput( m_StrainComponentsFilter->GetOutput( i ) );
        ostr.str( "" );
        ostr << m_OutputFilePrefix << "_Level_" << level << "_PreviousStrainComponent";
        ostr << i << ".mha";
        m_StrainComponentsWriter->SetFileName( ostr.str() );
        m_StrainComponentsWriter->Update();
        }
      }
    }
}

} // end namespace BlockMatching
} // end namespace itk

#endif
