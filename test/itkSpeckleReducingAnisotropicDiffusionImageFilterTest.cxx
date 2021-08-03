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

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMirrorPadImageFilter.h"
#include "itkTestingMacros.h"

#include "itkSpeckleReducingAnisotropicDiffusionImageFilter.h"

int
itkSpeckleReducingAnisotropicDiffusionImageFilterTest(int argc, char * argv[])
{
  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage outputImagePrefix";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }

  using PixelType = float;
  const unsigned int Dimension = 2;

  using ImageType = itk::Image<PixelType, Dimension>;

  ImageType::Pointer image = itk::ReadImage<ImageType>(argv[1]);

  // SRAD requires input image to be padded with 1 mirrored pixel
  // at each image border.
  using PadFilterType = itk::MirrorPadImageFilter<ImageType, ImageType>;
  PadFilterType::Pointer padFilter = PadFilterType::New();
  padFilter->SetInput(image);
  ImageType::SizeType padSize;
  padSize[0] = 1;
  padSize[1] = 1;
  padFilter->SetPadLowerBound(padSize);
  padFilter->SetPadUpperBound(padSize);
  padFilter->Update();

  using SpeckleFilterType = itk::SpeckleReducingAnisotropicDiffusionImageFilter<ImageType>;
  SpeckleFilterType::Pointer speckleFilter = SpeckleFilterType::New();
  speckleFilter->SetInput(padFilter->GetOutput());
  speckleFilter->SetNumberOfIterations(100);
  speckleFilter->SetTimeStep(0.002);

  ITK_EXERCISE_BASIC_OBJECT_METHODS(
    speckleFilter, SpeckleReducingAnisotropicDiffusionImageFilter, AnisotropicDiffusionImageFilter);

  ITK_TRY_EXPECT_NO_EXCEPTION(speckleFilter->Update());
  ImageType::Pointer result = speckleFilter->GetOutput();

  ITK_TEST_EXPECT_EQUAL(result->GetLargestPossibleRegion(), padFilter->GetOutput()->GetLargestPossibleRegion());

  ITK_TRY_EXPECT_NO_EXCEPTION(itk::WriteImage(result, argv[2]));

  return EXIT_SUCCESS;
}
