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
#include "itkArray.h"
#include "itkHDF5UltrasoundImageIO.h"
#include "itkHDF5UltrasoundImageIOFactory.h"
#include "itkMetaDataObject.h"
#include "itkTestingMacros.h"

int
itkHDF5UltrasoundImageIOTest( int argc, char * argv [] )
{
  if( argc < 2 )
    {
    std::cerr << "Usage: " << argv[0] << " inputImage" << std::endl;
    return EXIT_FAILURE;
    }
  const char * inputImageFileName = argv[1];

  itk::HDF5UltrasoundImageIO::Pointer imageIO = itk::HDF5UltrasoundImageIO::New();
  EXERCISE_BASIC_OBJECT_METHODS( imageIO, HDF5UltrasoundImageIO, StreamingImageIOBase );

  TEST_EXPECT_TRUE( !imageIO->CanReadFile( "AMINCFile.mnc" ) );
  TEST_EXPECT_TRUE( imageIO->CanReadFile( inputImageFileName ) );

  imageIO->SetFileName( inputImageFileName );

  imageIO->ReadImageInformation();

  const unsigned int Dimension = 3;

  itk::SizeValueType dimensions[Dimension];
  for( unsigned int ii = 0; ii < Dimension; ++ii )
    {
    dimensions[ii] = imageIO->GetDimensions( ii );
    }
  std::cout << "Dimensions: [ " << dimensions[0] << ", " << dimensions[1] << ", " << dimensions[2] << " ]" << std::endl;
  TEST_EXPECT_EQUAL( dimensions[0], 240 );
  TEST_EXPECT_EQUAL( dimensions[1], 328 );
  TEST_EXPECT_EQUAL( dimensions[2], 125 );

  std::cout << "ComponentType: " << imageIO->GetComponentTypeAsString( imageIO->GetComponentType() ) << std::endl;
  TEST_EXPECT_EQUAL( imageIO->GetComponentType(), 9 );

  const itk::MetaDataDictionary & metaDataDict = imageIO->GetMetaDataDictionary();
  std::string sliceType;
  itk::ExposeMetaData< std::string >( metaDataDict, "SliceType", sliceType );
  std::cout << "SliceType: " << sliceType << std::endl;
  TEST_EXPECT_EQUAL( sliceType, "Image" );

  return EXIT_SUCCESS;
}
