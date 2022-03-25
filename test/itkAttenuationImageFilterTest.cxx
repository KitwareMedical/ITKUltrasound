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
    std::cerr << " spectraImage maskImage outputImage <outputMaskImage> <numWorkUnits> <fixedEstimationDepthMM>";
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
  MaskImageType::Pointer    maskImage = itk::ReadImage<MaskImageType>(std::string(argv[2]));

  const std::string outputImagePath = argv[3];
  const std::string outputMaskImagePath = (argc > 4 ? argv[4] : "");

  unsigned int numWorkUnits = (argc > 5 ? std::atoi(argv[5]) : 1);
  float        fixedEstimationDepthMM = (argc > 6 ? std::atof(argv[6]) : 0.0);

  // Initialize the filter
  using AttenuationFilterType = itk::AttenuationImageFilter<SpectraImageType, OutputImageType, MaskImageType>;
  AttenuationFilterType::Pointer attenuationFilter = AttenuationFilterType::New();

  // Verify spatial distances are not settable without input image for reference
  ITK_TRY_EXPECT_EXCEPTION(attenuationFilter->SetFixedEstimationDepthMM(3.0));
  ITK_TRY_EXPECT_EXCEPTION(attenuationFilter->SetPadUpperBoundsMM(3.0));
  ITK_TRY_EXPECT_EXCEPTION(attenuationFilter->SetPadLowerBoundsMM(3.0));

  attenuationFilter->SetInput(inputImage);
  ITK_TEST_SET_GET_VALUE(inputImage, attenuationFilter->GetInput());

  // Verify filter does not run without mask
  ITK_TRY_EXPECT_EXCEPTION(attenuationFilter->Update());

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
  ITK_TRY_EXPECT_NO_EXCEPTION(attenuationFilter->Update());

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

  attenuationFilter->SetConsiderNegativeAttenuations(false);
  ITK_TEST_SET_GET_VALUE(false, attenuationFilter->GetConsiderNegativeAttenuations());

  attenuationFilter->SetPadLowerBounds(5);
  ITK_TEST_SET_GET_VALUE(5, attenuationFilter->GetPadLowerBounds());
  attenuationFilter->SetPadUpperBounds(2);
  ITK_TEST_SET_GET_VALUE(2, attenuationFilter->GetPadUpperBounds());

  attenuationFilter->SetPadLowerBoundsMM(0.1);
  ITK_TEST_EXPECT_TRUE(fabs(attenuationFilter->GetPadLowerBoundsMM() - 0.1) < (inputImage->GetSpacing()[0] / 2));
  attenuationFilter->SetPadUpperBoundsMM(0.0);
  ITK_TEST_EXPECT_TRUE(fabs(attenuationFilter->GetPadUpperBoundsMM() - 0.0) < (inputImage->GetSpacing()[0] / 2));

  attenuationFilter->SetNumberOfWorkUnits(numWorkUnits);

  ITK_EXERCISE_BASIC_OBJECT_METHODS(attenuationFilter, AttenuationImageFilter, ImageToImageFilter);

  // Run
  ITK_TRY_EXPECT_NO_EXCEPTION(attenuationFilter->Update());

  // Verify output
  itk::WriteImage(attenuationFilter->GetOutput(), outputImagePath);

  // Verify mask output after erosion
  if (outputMaskImagePath != "")
  {
    itk::WriteImage(attenuationFilter->GetOutputMaskImage(), outputMaskImagePath);
  }

  // discover mask's inside value
  using MaxType = itk::MinimumMaximumImageCalculator<MaskImageType>;
  auto maxCalculator = MaxType::New();
  maxCalculator->SetImage(attenuationFilter->GetOutputMaskImage());
  maxCalculator->ComputeMaximum();
  LabelType maxMaskValue = maxCalculator->GetMaximum();
  if (maxMaskValue == 0)
  {
    std::cerr << "The mask has been completely eroded away!" << std::endl;
    return EXIT_FAILURE;
  }

  // now boil it down to a single attenuation value
  using HistogramFilterType = itk::Statistics::MaskedImageToHistogramFilter<OutputImageType, MaskImageType>;
  auto histogramFilter = HistogramFilterType::New();
  histogramFilter->SetInput(attenuationFilter->GetOutput());
  histogramFilter->SetMaskImage(attenuationFilter->GetOutputMaskImage());
  histogramFilter->SetMarginalScale(10);

  HistogramFilterType::HistogramSizeType histogramSize{ 1 };
  histogramSize[0] = 1e5;
  histogramFilter->SetHistogramSize(histogramSize);

  histogramFilter->SetMaskValue(maxMaskValue); // this is most likely the only label present
  histogramFilter->Update();
  auto histogram = histogramFilter->GetOutput();

  // we will use median as a robust estimate of the mean
  auto median = histogram->Quantile(0, 0.50);
  std::cout << "Median attenuation: " << median << " dB/(MHz*cm)" << std::endl;

  if (!std::isfinite(median))
  {
    std::cerr << "The median attenuation if not a finite number! It is: " << median << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
