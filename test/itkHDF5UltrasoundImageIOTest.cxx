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
#include "itkMath.h"
#include "itkImageIORegion.h"

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
  ITK_EXERCISE_BASIC_OBJECT_METHODS( imageIO, HDF5UltrasoundImageIO, StreamingImageIOBase );

  ITK_TEST_EXPECT_TRUE( !imageIO->CanReadFile( "AMINCFile.mnc" ) );
  ITK_TEST_EXPECT_TRUE( imageIO->CanReadFile( inputImageFileName ) );

  imageIO->SetFileName( inputImageFileName );

  imageIO->ReadImageInformation();

  const unsigned int Dimension = 3;

  itk::SizeValueType dimensions[Dimension];
  for( unsigned int ii = 0; ii < Dimension; ++ii )
    {
    dimensions[ii] = imageIO->GetDimensions( ii );
    }
  std::cout << "Dimensions: [ " << dimensions[0] << ", " << dimensions[1] << ", " << dimensions[2] << " ]" << std::endl;
  ITK_TEST_EXPECT_EQUAL( dimensions[0], 240 );
  ITK_TEST_EXPECT_EQUAL( dimensions[1], 328 );
  ITK_TEST_EXPECT_EQUAL( dimensions[2], 125 );

  std::cout << "ComponentType: " << imageIO->GetComponentTypeAsString( imageIO->GetComponentType() ) << std::endl;
  ITK_TEST_EXPECT_EQUAL( imageIO->GetComponentType(), itk::IOComponentEnum::FLOAT );

  const itk::MetaDataDictionary & metaDataDict = imageIO->GetMetaDataDictionary();
  std::string sliceType;
  itk::ExposeMetaData< std::string >( metaDataDict, "SliceType", sliceType );
  std::cout << "SliceType: " << sliceType << std::endl;
  ITK_TEST_EXPECT_EQUAL( sliceType, "Image" );

  typedef itk::Array< double > SliceSpacingType;
  SliceSpacingType sliceSpacing( 2 );
  itk::ExposeMetaData< SliceSpacingType >( metaDataDict, "SliceSpacing", sliceSpacing );
  std::cout << "SliceSpacing: [ " << sliceSpacing[0] << ", " << sliceSpacing[1] << " ]" << std::endl;
  ITK_TEST_EXPECT_TRUE( itk::Math::FloatAlmostEqual( sliceSpacing[0], 0.1925 ) );
  ITK_TEST_EXPECT_TRUE( itk::Math::FloatAlmostEqual( sliceSpacing[1], 0.167811, 10, 1e-6 ) );

  typedef itk::Array< double > SliceOriginType;
  SliceOriginType sliceOrigin( 2 );
  itk::ExposeMetaData< SliceOriginType >( metaDataDict, "SliceOrigin", sliceOrigin );
  std::cout << "SliceOrigin: [ " << sliceOrigin[0] << ", " << sliceOrigin[1] << " ]" << std::endl;
  ITK_TEST_EXPECT_TRUE( itk::Math::FloatAlmostEqual( sliceOrigin[0], 0.0 ) );
  ITK_TEST_EXPECT_TRUE( itk::Math::FloatAlmostEqual( sliceOrigin[1], -27.2693, 10, 1e-3 ) );

  typedef itk::Array< double > ElevationalSliceAnglesType;
  ElevationalSliceAnglesType elevationalSliceAngles( imageIO->GetDimensions( 2 ) );
  itk::ExposeMetaData< ElevationalSliceAnglesType >( metaDataDict, "ElevationalSliceAngles", elevationalSliceAngles );
  std::cout << "ElevationalSliceAngles: [ " << elevationalSliceAngles[0] << ", " << elevationalSliceAngles[1] << ", " << elevationalSliceAngles[2] << " ..." << std::endl;
  ITK_TEST_EXPECT_TRUE( itk::Math::FloatAlmostEqual( elevationalSliceAngles[0], -1.0821, 10, 1e-3 ) );
  ITK_TEST_EXPECT_TRUE( itk::Math::FloatAlmostEqual( elevationalSliceAngles[1], -1.06465, 10, 1e-4 ) );


  float * buffer = new float[10 * 10 * 10];

  itk::ImageIORegion ioRegion( Dimension );
  for( unsigned int ii = 0; ii < Dimension; ++ii )
    {
    ioRegion.SetIndex( ii, 0 );
    ioRegion.SetSize( ii, 10 );
    }
  imageIO->SetIORegion( ioRegion );
  imageIO->Read( static_cast< void * >( buffer ) );
  ITK_TEST_EXPECT_EQUAL( buffer[0], 88.0 );
  ITK_TEST_EXPECT_EQUAL( buffer[1], 78.0 );
  ITK_TEST_EXPECT_EQUAL( buffer[2], 77.0 );

  delete[] buffer;

  return EXIT_SUCCESS;
}
