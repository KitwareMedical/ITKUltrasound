/*=========================================================================
 *
 *  Copyright NumFOCUS
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         https://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#include "itkHDF5UltrasoundImageIO.h"
#include "itkTestingMacros.h"

int
itkHDF5UltrasoundImageIOCanReadITKImageTest(int argc, char * argv[])
{
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " inputImage" << std::endl;
    return EXIT_FAILURE;
  }
  const char * inputImageFileName = argv[1];

  itk::HDF5UltrasoundImageIO::Pointer imageIO = itk::HDF5UltrasoundImageIO::New();

  // Do not read an HDF5 image in the ITK format -- let the ITK HDF5ImageIO do
  // that instead.
  ITK_TEST_EXPECT_TRUE(!imageIO->CanReadFile(inputImageFileName));

  return EXIT_SUCCESS;
}
