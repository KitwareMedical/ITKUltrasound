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

#include "itkCurvilinearArraySpecialCoordinatesImage.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkTestingMacros.h"


int
itkCurvilinearArraySpecialCoordinatesImageTest(int argc, char * argv[])
{
  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage outputImage";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }
  const char * inputImageFileName = argv[1];
  const char * outputImageFileName = argv[2];

  const unsigned int Dimension = 3;
  using PixelType = unsigned char;

  using ImageType = itk::Image<PixelType, Dimension>;
  using SpecialCoordinatesImageType = itk::CurvilinearArraySpecialCoordinatesImage<PixelType, Dimension>;

  using ReaderType = itk::ImageFileReader<SpecialCoordinatesImageType>;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inputImageFileName);
  ITK_TRY_EXPECT_NO_EXCEPTION(reader->Update());

  SpecialCoordinatesImageType::Pointer        curvilinearArrayImage = reader->GetOutput();
  const SpecialCoordinatesImageType::SizeType inputSize = curvilinearArrayImage->GetLargestPossibleRegion().GetSize();
  const double lateralAngularSeparation = (vnl_math::pi / 2.0 + 0.5) / (inputSize[1] - 1);
  std::cout << "lateral angular separation: " << lateralAngularSeparation << std::endl;
  curvilinearArrayImage->SetLateralAngularSeparation(lateralAngularSeparation);
  const double radiusStart = 26.4;
  const double radiusStop = 131.5;
  curvilinearArrayImage->SetFirstSampleDistance(radiusStart);
  curvilinearArrayImage->SetRadiusSampleSize((radiusStop - radiusStart) / (inputSize[0] - 1));

  curvilinearArrayImage->Print(std::cout);

  using ContinuousIndexType = itk::ContinuousIndex<double, Dimension>;
  ContinuousIndexType continuousIndex;
  continuousIndex[0] = 0.0;
  continuousIndex[1] = 0.0;
  continuousIndex[2] = 0.0;
  std::cout << "Continuous index: " << continuousIndex << std::endl;
  SpecialCoordinatesImageType::PointType point;
  curvilinearArrayImage->TransformContinuousIndexToPhysicalPoint(continuousIndex, point);
  std::cout << "Transformed point: " << point << std::endl;
  ContinuousIndexType continuousTransformedIndex;
  curvilinearArrayImage->TransformPhysicalPointToContinuousIndex(point, continuousTransformedIndex);
  std::cout << "Transformed continuous index: " << continuousTransformedIndex << std::endl;
  std::cout << std::endl;

  continuousIndex[0] = 2047.0;
  continuousIndex[1] = 240.0;
  continuousIndex[2] = 0.0;
  std::cout << "Continuous index: " << continuousIndex << std::endl;
  curvilinearArrayImage->TransformContinuousIndexToPhysicalPoint(continuousIndex, point);
  std::cout << "Transformed point: " << point << std::endl;
  curvilinearArrayImage->TransformPhysicalPointToContinuousIndex(point, continuousTransformedIndex);
  std::cout << "Transformed continuous index: " << continuousTransformedIndex << std::endl;
  std::cout << std::endl;

  using ResamplerType = itk::ResampleImageFilter<SpecialCoordinatesImageType, ImageType>;
  ResamplerType::Pointer resampler = ResamplerType::New();
  resampler->SetInput(reader->GetOutput());
  SpecialCoordinatesImageType::SizeType outputSize = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
  outputSize[0] = 800;
  outputSize[1] = 800;
  resampler->SetSize(outputSize);
  SpecialCoordinatesImageType::SpacingType outputSpacing;
  outputSpacing.Fill(0.15);
  resampler->SetOutputSpacing(outputSpacing);
  SpecialCoordinatesImageType::PointType outputOrigin;
  outputOrigin[0] = outputSize[0] * outputSpacing[0] / -2.0;
  outputOrigin[1] = radiusStart * std::cos(vnl_math::pi / 4.0);
  outputOrigin[2] = reader->GetOutput()->GetOrigin()[2];
  resampler->SetOutputOrigin(outputOrigin);

  using WindowedSincInterpolatorType = itk::WindowedSincInterpolateImageFunction<SpecialCoordinatesImageType, 3>;
  WindowedSincInterpolatorType::Pointer sincInterpolator = WindowedSincInterpolatorType::New();
  sincInterpolator->SetInputImage(curvilinearArrayImage);

  using WriterType = itk::ImageFileWriter<ImageType>;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(outputImageFileName);
  writer->SetInput(resampler->GetOutput());
  ITK_TRY_EXPECT_NO_EXCEPTION(writer->Update());

  // Check CopyInformation
  SpecialCoordinatesImageType::Pointer curvilinearArrayCopiedInformation = SpecialCoordinatesImageType::New();
  curvilinearArrayCopiedInformation->CopyInformation(curvilinearArrayImage);

  std::cout << "With copied information: " << curvilinearArrayCopiedInformation << std::endl;
  ITK_TEST_EXPECT_TRUE(itk::Math::FloatAlmostEqual(
    curvilinearArrayCopiedInformation->GetLateralAngularSeparation(), 0.00862832, 10, 1e-6));
  ITK_TEST_EXPECT_TRUE(
    itk::Math::FloatAlmostEqual(curvilinearArrayCopiedInformation->GetRadiusSampleSize(), 0.0513434, 10, 1e-6));
  ITK_TEST_EXPECT_TRUE(
    itk::Math::FloatAlmostEqual(curvilinearArrayCopiedInformation->GetFirstSampleDistance(), 26.4, 10, 1e-6));

  return EXIT_SUCCESS;
}
