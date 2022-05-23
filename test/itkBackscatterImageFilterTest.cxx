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

#include <string>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTestingMacros.h"

#include "itkBackscatterImageFilter.h"

int
itkBackscatterImageFilterTest(int argc, char * argv[])
{
  if (argc < 5)
  {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " spectraImage averageImage slopeImage interceptImage <numWorkUnits>";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }

  using RealType = float;
  const unsigned int Dimension = 3;

  using SpectraImageType = itk::VectorImage<RealType, Dimension>;
  using OutputImageType = itk::Image<RealType, Dimension>;

  SpectraImageType::Pointer inputImage = itk::ReadImage<SpectraImageType>(std::string(argv[1]));

  const std::string averageImage = argv[2];
  const std::string slopeImage = argv[3];
  const std::string interceptImage = argv[4];

  unsigned int numWorkUnits = (argc > 5 ? std::stoi(argv[5]) : 1);

  // Initialize the filter
  using BackscatterFilterType = itk::BackscatterImageFilter<SpectraImageType, OutputImageType>;
  BackscatterFilterType::Pointer backscatterFilter = BackscatterFilterType::New();

  backscatterFilter->SetInput(inputImage);
  ITK_TEST_SET_GET_VALUE(inputImage, backscatterFilter->GetInput());

  // Missing sampling frequency produces an error
  ITK_TRY_EXPECT_EXCEPTION(backscatterFilter->Update());

  backscatterFilter->SetSamplingFrequencyMHz(60);
  ITK_TEST_SET_GET_VALUE(60.0, backscatterFilter->GetSamplingFrequencyMHz());
  ITK_TRY_EXPECT_EXCEPTION(backscatterFilter->UpdateOutputInformation());

  backscatterFilter->SetFrequencyBandStartMHz(5.0);
  ITK_TEST_SET_GET_VALUE(5.0, backscatterFilter->GetFrequencyBandStartMHz());
  backscatterFilter->SetFrequencyBandEndMHz(20.0);
  ITK_TEST_SET_GET_VALUE(20.0, backscatterFilter->GetFrequencyBandEndMHz());

  backscatterFilter->SetNumberOfWorkUnits(numWorkUnits);

  ITK_EXERCISE_BASIC_OBJECT_METHODS(backscatterFilter, BackscatterImageFilter, ImageToImageFilter);

  ITK_TRY_EXPECT_NO_EXCEPTION(backscatterFilter->Update());

  itk::WriteImage(backscatterFilter->GetOutput(0), averageImage, false);
  itk::WriteImage(backscatterFilter->GetOutput(1), slopeImage, false);
  itk::WriteImage(backscatterFilter->GetOutput(2), interceptImage, false);

  return EXIT_SUCCESS;
}
