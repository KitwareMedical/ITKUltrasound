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

#include <string>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTestingMacros.h"

#include "itkAttenuationImageFilter.h"
#include "itkMaskedImageToHistogramFilter.h"
#include "itkMinimumMaximumImageCalculator.h"

int
itkAttenuationImageFilterTest(int argc, char * argv[])
{
  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " spectraImage maskImage outputImage";
    std::cerr << " <outputMaskImage> <numWorkUnits> <fixedEstimationDepthMM>";
    std::cerr << " <considerNegativeAttenuations> <computationMode> <expectedAttenuation>";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }

  using RealType = float;
  using LabelType = unsigned char;
  const unsigned int Dimension = 3;

  using SpectraImageType = itk::VectorImage<RealType, Dimension>;
  using OutputImageType = itk::Image<RealType, Dimension>;
  using MaskImageType = itk::Image<LabelType, Dimension>;

  SpectraImageType::Pointer inputImage = itk::ReadImage<SpectraImageType>(std::string(argv[1]));
  MaskImageType::Pointer    maskImage;
  std::string               maskFilename = argv[2];
  if (maskFilename != "nul")
  {
    maskImage = itk::ReadImage<MaskImageType>(maskFilename);
  }

  const std::string outputImagePath = argv[3];
  const std::string outputMaskImagePath = (argc > 4 ? argv[4] : "");

  unsigned int numWorkUnits = (argc > 5 ? std::stoi(argv[5]) : 1);
  float        fixedEstimationDepthMM = (argc > 6 ? std::stof(argv[6]) : 0.0);
  bool         considerNegativeAttenuations = (argc > 7 ? std::stoul(argv[7]) : false);
  unsigned int computationMode = (argc > 8 ? std::stoul(argv[8]) : 0);

  // Initialize the filter
  using AttenuationFilterType = itk::AttenuationImageFilter<SpectraImageType, OutputImageType, MaskImageType>;
  AttenuationFilterType::Pointer attenuationFilter = AttenuationFilterType::New();

  // Verify spatial distances are not settable without input image for reference
  ITK_TRY_EXPECT_EXCEPTION(attenuationFilter->SetFixedEstimationDepthMM(3.0));
  ITK_TRY_EXPECT_EXCEPTION(attenuationFilter->SetPadUpperBoundsMM(3.0));
  ITK_TRY_EXPECT_EXCEPTION(attenuationFilter->SetPadLowerBoundsMM(3.0));

  attenuationFilter->SetInput(inputImage);
  ITK_TEST_SET_GET_VALUE(inputImage, attenuationFilter->GetInput());

  attenuationFilter->SetInputMaskImage(maskImage);
  ITK_TEST_SET_GET_VALUE(maskImage, attenuationFilter->GetInputMaskImage());

  // Bad scanline direction produces an error
  attenuationFilter->SetScanDirection(Dimension);
  ITK_TRY_EXPECT_EXCEPTION(attenuationFilter->Update());

  attenuationFilter->SetScanDirection(0);
  ITK_TEST_SET_GET_VALUE(0, attenuationFilter->GetScanDirection());

  // Missing sampling frequency produces an error
  ITK_TRY_EXPECT_EXCEPTION(attenuationFilter->Update());

  attenuationFilter->SetSamplingFrequencyMHz(60);
  ITK_TEST_SET_GET_VALUE(60.0, attenuationFilter->GetSamplingFrequencyMHz());
  ITK_TRY_EXPECT_NO_EXCEPTION(attenuationFilter->UpdateOutputInformation());

  attenuationFilter->SetFixedEstimationDepth(3);
  ITK_TEST_SET_GET_VALUE(3, attenuationFilter->GetFixedEstimationDepth());

  // Verify spatial distance is estimated to the nearest pixel center
  attenuationFilter->SetFixedEstimationDepthMM(fixedEstimationDepthMM);
  ITK_TEST_EXPECT_TRUE(fabs(attenuationFilter->GetFixedEstimationDepthMM() - fixedEstimationDepthMM) <
                       (inputImage->GetSpacing()[0] / 2));

  attenuationFilter->SetFrequencyBandStartMHz(5.0);
  ITK_TEST_SET_GET_VALUE(5.0, attenuationFilter->GetFrequencyBandStartMHz());
  attenuationFilter->SetFrequencyBandEndMHz(20.0);
  ITK_TEST_SET_GET_VALUE(20.0, attenuationFilter->GetFrequencyBandEndMHz());

  attenuationFilter->SetConsiderNegativeAttenuations(considerNegativeAttenuations);
  ITK_TEST_SET_GET_VALUE(considerNegativeAttenuations, attenuationFilter->GetConsiderNegativeAttenuations());

  attenuationFilter->SetPadLowerBounds(5);
  ITK_TEST_SET_GET_VALUE(5, attenuationFilter->GetPadLowerBounds());
  attenuationFilter->SetPadUpperBounds(2);
  ITK_TEST_SET_GET_VALUE(2, attenuationFilter->GetPadUpperBounds());

  attenuationFilter->SetPadLowerBoundsMM(0.1);
  ITK_TEST_EXPECT_TRUE(fabs(attenuationFilter->GetPadLowerBoundsMM() - 0.1) < (inputImage->GetSpacing()[0] / 2));
  attenuationFilter->SetPadUpperBoundsMM(0.0);
  ITK_TEST_EXPECT_TRUE(fabs(attenuationFilter->GetPadUpperBoundsMM() - 0.0) < (inputImage->GetSpacing()[0] / 2));

  attenuationFilter->SetNumberOfWorkUnits(numWorkUnits);
  attenuationFilter->SetComputationMode(computationMode);

  ITK_EXERCISE_BASIC_OBJECT_METHODS(attenuationFilter, AttenuationImageFilter, ImageToImageFilter);

  // Run
  ITK_TRY_EXPECT_NO_EXCEPTION(attenuationFilter->Update());

  // Verify output
  itk::WriteImage(attenuationFilter->GetOutput(), outputImagePath, false);

  // Verify mask output after erosion
  if (outputMaskImagePath != "")
  {
    itk::WriteImage(attenuationFilter->GetOutputMaskImage(), outputMaskImagePath, true);
  }

  LabelType maxMaskValue = 1;
  if (maskFilename != "nul")
  {
    // Discover mask's inside value
    using MaxType = itk::MinimumMaximumImageCalculator<MaskImageType>;
    auto maxCalculator = MaxType::New();
    maxCalculator->SetImage(maskImage);
    maxCalculator->ComputeMaximum();
    LabelType maxMaskValue = maxCalculator->GetMaximum();
    if (maxMaskValue == 0)
    {
      std::cerr << "The mask is empty!" << std::endl;
      return EXIT_FAILURE;
    }
  }

  // Now boil it down to a single attenuation value
  using HistogramFilterType = itk::Statistics::MaskedImageToHistogramFilter<OutputImageType, MaskImageType>;
  auto histogramFilter = HistogramFilterType::New();
  histogramFilter->SetInput(attenuationFilter->GetOutput());
  histogramFilter->SetMaskImage(attenuationFilter->GetOutputMaskImage());
  histogramFilter->SetMarginalScale(10);

  HistogramFilterType::HistogramSizeType histogramSize{ 1 };
  histogramSize[0] = 1e5;
  histogramFilter->SetHistogramSize(histogramSize);

  histogramFilter->SetMaskValue(maxMaskValue); // This is most likely the only label present
  histogramFilter->Update();
  auto histogram = histogramFilter->GetOutput();

  // We will use median as a robust estimate of the mean
  float median = histogram->Quantile(0, 0.50);
  std::cout << "Median attenuation: " << median << " dB/(MHz*cm)" << std::endl;

  if (!std::isfinite(median))
  {
    std::cerr << "The median attenuation if not a finite number! It is: " << median << std::endl;
    return EXIT_FAILURE;
  }

  if (argc > 9) // Expected value is provided on the command line
  {
    float expectedAttenuation = std::stof(argv[9]);
    if (!itk::Math::FloatAlmostEqual(expectedAttenuation, median, 4, 1e-4))
    {
      std::cerr << "Regression test failure: the expected attenuation is: " << expectedAttenuation << std::endl;
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}
