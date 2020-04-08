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

#include "itkUltrasoundImageFileReader.h"
#include "itkHDF5UltrasoundImageIO.h"
#include "itkHDF5UltrasoundImageIOFactory.h"

#include "itkTestingMacros.h"
#include "itkSliceSeriesSpecialCoordinatesImage.h"
#include "itkEuler3DTransform.h"

int
itkHDF5BModeUltrasoundImageFileReaderTest(int argc, char * argv[])
{
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }
  const char * inputImageFileName = argv[1];

  itk::HDF5UltrasoundImageIOFactory::RegisterOneFactory();

  const unsigned int Dimension = 3;
  using PixelType = unsigned char;
  using ParametersValueType = double;

  using SliceImageType = itk::Image<PixelType, Dimension - 1>;
  using TransformType = itk::Euler3DTransform<ParametersValueType>;
  using SpecialCoordinatesImageType =
    itk::SliceSeriesSpecialCoordinatesImage<SliceImageType, TransformType, PixelType, Dimension>;

  using ReaderType = itk::UltrasoundImageFileReader<SpecialCoordinatesImageType>;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inputImageFileName);
  ITK_TRY_EXPECT_NO_EXCEPTION(reader->Update());

  SpecialCoordinatesImageType::ConstPointer image = reader->GetOutput();
  std::cout << image << std::endl;

  SpecialCoordinatesImageType::SliceImageType::ConstPointer sliceImage = image->GetSliceImage();

  const SpecialCoordinatesImageType::SliceImageType::SpacingType sliceSpacing = sliceImage->GetSpacing();
  ITK_TEST_EXPECT_TRUE(itk::Math::FloatAlmostEqual(sliceSpacing[0], 0.1925));
  ITK_TEST_EXPECT_TRUE(itk::Math::FloatAlmostEqual(sliceSpacing[1], 0.167811, 10, 1e-6));

  const SpecialCoordinatesImageType::SliceImageType::PointType sliceOrigin = sliceImage->GetOrigin();
  ITK_TEST_EXPECT_TRUE(itk::Math::FloatAlmostEqual(sliceOrigin[0], 0.0));
  ITK_TEST_EXPECT_TRUE(itk::Math::FloatAlmostEqual(sliceOrigin[1], -27.2693, 10, 1e-3));

  const TransformType * transform = image->GetSliceTransform(0);
  transform->Print(std::cout);
  ITK_TEST_EXPECT_TRUE(itk::Math::FloatAlmostEqual(transform->GetAngleY(), -1.0821, 10, 1e-3));

  return EXIT_SUCCESS;
}
