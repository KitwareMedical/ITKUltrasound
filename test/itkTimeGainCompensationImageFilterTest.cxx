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

#include "itkTimeGainCompensationImageFilter.h"
#include "itkCurvilinearArraySpecialCoordinatesImage.h"
#include "itkBModeImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkTestingMacros.h"

int
itkTimeGainCompensationImageFilterTest(int argc, char * argv[])
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

  const unsigned int Dimension = 2;
  using IntegerPixelType = signed short;
  using IntegerImageType = itk::CurvilinearArraySpecialCoordinatesImage<IntegerPixelType, Dimension>;
  using RealPixelType = float;
  using RealImageType = itk::CurvilinearArraySpecialCoordinatesImage<RealPixelType, Dimension>;
  using ScanConvertedImageType = itk::Image<RealPixelType, Dimension>;

  using ReaderType = itk::ImageFileReader<IntegerImageType>;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inputImageFileName);

  using TGCFilterType = itk::TimeGainCompensationImageFilter<IntegerImageType>;
  TGCFilterType::Pointer tgcFilter = TGCFilterType::New();
  tgcFilter->SetInput(reader->GetOutput());

  using GainType = TGCFilterType::GainType;

  GainType gain = tgcFilter->GetGain();
  ITK_TEST_SET_GET_VALUE(1.0, gain(0, 1));

  // Invalid number of columns
  gain.SetSize(4, 3);
  tgcFilter->SetGain(gain);
  ITK_TRY_EXPECT_EXCEPTION(tgcFilter->Update());

  // Invalid number of rows
  gain.SetSize(1, 2);
  tgcFilter->SetGain(gain);
  ITK_TRY_EXPECT_EXCEPTION(tgcFilter->Update());

  // Depths are not ascending
  gain.SetSize(3, 2);
  gain(0, 0) = 0.0;
  gain(0, 1) = 1.0;
  gain(1, 0) = 2000.0;
  gain(1, 1) = 3.0;
  gain(2, 0) = 1000.0;
  gain(2, 1) = 5.0;
  tgcFilter->SetGain(gain);
  ITK_TRY_EXPECT_EXCEPTION(tgcFilter->Update());

  gain(1, 0) = 1000.0;
  gain(2, 0) = 2000.0;
  tgcFilter->SetGain(gain);

  using CasterType = itk::CastImageFilter<IntegerImageType, RealImageType>;
  CasterType::Pointer caster = CasterType::New();
  caster->SetInput(tgcFilter->GetOutput());

  using BModeFilterType = itk::BModeImageFilter<RealImageType, RealImageType>;
  BModeFilterType::Pointer bmodeFilter = BModeFilterType::New();
  bmodeFilter->SetInput(caster->GetOutput());

  ITK_TRY_EXPECT_NO_EXCEPTION(bmodeFilter->Update());

  tgcFilter->Print(std::cout);

  RealImageType::Pointer curvilinearArrayImage = bmodeFilter->GetOutput();
  curvilinearArrayImage->DisconnectPipeline();
  const RealImageType::SizeType inputSize = curvilinearArrayImage->GetLargestPossibleRegion().GetSize();
  const double                  lateralAngularSeparation = (vnl_math::pi / 2.0 + 0.5) / (inputSize[1] - 1);
  curvilinearArrayImage->SetLateralAngularSeparation(lateralAngularSeparation);
  const double radiusStart = 7.0;
  const double radiusStop = 112.5;
  curvilinearArrayImage->SetFirstSampleDistance(radiusStart);
  curvilinearArrayImage->SetRadiusSampleSize((radiusStop - radiusStart) / (inputSize[0] - 1));

  using ResamplerType = itk::ResampleImageFilter<RealImageType, ScanConvertedImageType>;
  ResamplerType::Pointer resampler = ResamplerType::New();
  resampler->SetInput(curvilinearArrayImage);
  RealImageType::SizeType outputSize;
  outputSize[0] = 400;
  outputSize[1] = 400;
  resampler->SetSize(outputSize);
  RealImageType::SpacingType outputSpacing;
  outputSpacing.Fill(0.30);
  resampler->SetOutputSpacing(outputSpacing);
  RealImageType::PointType outputOrigin;
  outputOrigin[0] = outputSize[0] * outputSpacing[0] / -2.0;
  outputOrigin[1] = radiusStart * std::cos(vnl_math::pi / 4.0);
  resampler->SetOutputOrigin(outputOrigin);

  using OutputImageType = itk::Image<unsigned char, Dimension>;
  using RescalerType = itk::RescaleIntensityImageFilter<ScanConvertedImageType, OutputImageType>;
  RescalerType::Pointer rescaler = RescalerType::New();
  rescaler->SetInput(resampler->GetOutput());

  using WriterType = itk::ImageFileWriter<OutputImageType>;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(outputImageFileName);
  writer->SetInput(rescaler->GetOutput());
  ITK_TRY_EXPECT_NO_EXCEPTION(writer->Update());

  return EXIT_SUCCESS;
}
