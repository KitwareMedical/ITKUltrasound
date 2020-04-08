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

#include "itkCurvilinearArraySpecialCoordinatesImage.h"

#include "itkUltrasoundImageFileReader.h"

#include "itkTestingMacros.h"

int itkCurvilinearArrayUltrasoundImageFileReaderTest( int argc, char * argv [] )
{
  if( argc < 2 )
    {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }
  const char * inputImageFileName = argv[1];

  const unsigned int Dimension = 3;
  using PixelType = unsigned char;

  using SpecialCoordinatesImageType = itk::CurvilinearArraySpecialCoordinatesImage< PixelType, Dimension >;

  using ReaderType = itk::UltrasoundImageFileReader< SpecialCoordinatesImageType >;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImageFileName );
  ITK_TRY_EXPECT_NO_EXCEPTION( reader->Update() );

  SpecialCoordinatesImageType::ConstPointer image = reader->GetOutput();
  ITK_TEST_EXPECT_EQUAL( image->GetLateralAngularSeparation(), 0.00862832 );
  ITK_TEST_EXPECT_EQUAL( image->GetRadiusSampleSize(), 0.0513434294 );
  ITK_TEST_EXPECT_EQUAL( image->GetFirstSampleDistance(), 26.4 );

  return EXIT_SUCCESS;
}
