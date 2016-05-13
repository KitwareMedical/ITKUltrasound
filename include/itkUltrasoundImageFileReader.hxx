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
#ifndef itkUltrasoundImageFileReader_hxx
#define itkUltrasoundImageFileReader_hxx

#include "itkUltrasoundImageFileReader.h"

#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"

#include "itkCurvilinearArraySpecialCoordinatesImage.h"

namespace
{

double
GetMetaDataAsDouble(const itk::MetaDataDictionary & dictionary, const std::string & key)
{
  std::string valueAsStr;
  itk::ExposeMetaData< std::string >( dictionary, key, valueAsStr );
  std::istringstream istrm( valueAsStr );
  double value;
  istrm >> value;
  return value;
}


template< typename TImage >
class Dispatch
{
public:
  typedef TImage ImageType;
  static void ExtractMetaData( TImage * image ) {};
};


template< typename TPixel, unsigned int VDimension >
class Dispatch< itk::CurvilinearArraySpecialCoordinatesImage< TPixel, VDimension > >
{
public:
  typedef itk::CurvilinearArraySpecialCoordinatesImage< TPixel, VDimension > ImageType;
  static void ExtractMetaData( ImageType * image )
  {
    const itk::MetaDataDictionary & dictionary = image->GetMetaDataDictionary();
    if( dictionary.HasKey( "LateralAngularSeparation" ) )
      {
      image->SetLateralAngularSeparation( GetMetaDataAsDouble( dictionary, "LateralAngularSeparation" ) );
      }
    if( dictionary.HasKey( "FirstSampleDistance" ) )
      {
      image->SetFirstSampleDistance( GetMetaDataAsDouble( dictionary, "FirstSampleDistance" ) );
      }
    if( dictionary.HasKey( "RadiusSampleSize" ) )
      {
      image->SetRadiusSampleSize( GetMetaDataAsDouble( dictionary, "RadiusSampleSize" ) );
      }
  }
};

} // end anonymous namespace

namespace itk
{

template < typename TOutputImage >
UltrasoundImageFileReader< TOutputImage >
::UltrasoundImageFileReader()
{
}


template < typename TOutputImage >
void
UltrasoundImageFileReader< TOutputImage >
::GenerateOutputInformation()
{
  Superclass::GenerateOutputInformation();
  Dispatch< TOutputImage >::ExtractMetaData( this->GetOutput() );
}


} // end namespace itk

#endif
