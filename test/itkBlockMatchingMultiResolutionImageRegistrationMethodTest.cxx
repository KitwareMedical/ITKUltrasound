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
#include "itkVector.h"

#include "itkBlockMatchingImageRegistrationMethod.h"
#include "itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter.h"
#include "itkBlockMatchingMultiResolutionFixedBlockRadiusCalculator.h"
#include "itkBlockMatchingMultiResolutionFixedSearchRegionImageSource.h"
#include "itkBlockMatchingMultiResolutionImageRegistrationMethod.h"

int
itkBlockMatchingMultiResolutionImageRegistrationMethodTest(int argc, char * argv[])
{
  if (argc < 4)
  {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixedImage movingImage displacementImage";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }

  const unsigned int Dimension = 2;
  using InputPixelType = signed short;
  using InputImageType = itk::Image<InputPixelType, Dimension>;
  using RadiusType = InputImageType::SizeType;

  using MetricPixelType = double;
  using MetricImageType = itk::Image<MetricPixelType, Dimension>;

  using VectorType = itk::Vector<MetricPixelType, Dimension>;
  using DisplacementImageType = itk::Image<VectorType, Dimension>;

  using CoordRepType = double;

  using ReaderType = itk::ImageFileReader<InputImageType>;
  ReaderType::Pointer fixedReader = ReaderType::New();
  fixedReader->SetFileName(argv[1]);
  ReaderType::Pointer movingReader = ReaderType::New();
  movingReader->SetFileName(argv[2]);

  using BlockRadiusCalculatorType = itk::BlockMatching::MultiResolutionFixedBlockRadiusCalculator<InputImageType>;
  BlockRadiusCalculatorType::Pointer blockRadiusCalculator = BlockRadiusCalculatorType::New();
  RadiusType                         blockRadius;
  blockRadius[0] = 12;
  blockRadius[1] = 4;
  blockRadiusCalculator->SetRadius(blockRadius);


  using SearchRegionImageSourceType = itk::BlockMatching::
    MultiResolutionFixedSearchRegionImageSource<InputImageType, InputImageType, DisplacementImageType>;
  SearchRegionImageSourceType::Pointer             searchRegionSource = SearchRegionImageSourceType::New();
  SearchRegionImageSourceType::PyramidScheduleType pyramidSchedule(3, Dimension);
  pyramidSchedule(0, 0) = 3;
  pyramidSchedule(0, 1) = 2;
  pyramidSchedule(1, 0) = 2;
  pyramidSchedule(1, 1) = 1;
  pyramidSchedule(2, 0) = 1;
  pyramidSchedule(2, 1) = 1;
  searchRegionSource->SetPyramidSchedule(pyramidSchedule);
  RadiusType searchRadius;
  searchRadius[0] = 50;
  searchRadius[1] = 6;
  searchRegionSource->SetSearchRegionRadiusSchedule(searchRadius);
  searchRegionSource->SetOverlapSchedule(1.0);

  using MetricImageFilterType = itk::BlockMatching::
    NormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter<InputImageType, InputImageType, MetricImageType>;
  MetricImageFilterType::Pointer metricImageFilter = MetricImageFilterType::New();

  using LevelRegistrationMethodType = itk::BlockMatching::
    ImageRegistrationMethod<InputImageType, InputImageType, MetricImageType, DisplacementImageType, CoordRepType>;
  LevelRegistrationMethodType::Pointer levelRegistrationMethod = LevelRegistrationMethodType::New();
  levelRegistrationMethod->SetMetricImageFilter(metricImageFilter);

  using RegistrationMethodType = itk::BlockMatching::MultiResolutionImageRegistrationMethod<InputImageType,
                                                                                            InputImageType,
                                                                                            MetricImageType,
                                                                                            DisplacementImageType,
                                                                                            CoordRepType>;
  RegistrationMethodType::Pointer multiResRegistrationMethod = RegistrationMethodType::New();
  multiResRegistrationMethod->SetFixedImage(fixedReader->GetOutput());
  multiResRegistrationMethod->SetMovingImage(movingReader->GetOutput());
  multiResRegistrationMethod->SetBlockRadiusCalculator(blockRadiusCalculator);
  multiResRegistrationMethod->SetSearchRegionImageSource(searchRegionSource);
  multiResRegistrationMethod->SetSchedules(pyramidSchedule, pyramidSchedule);
  multiResRegistrationMethod->SetImageRegistrationMethod(levelRegistrationMethod);

  using WriterType = itk::ImageFileWriter<DisplacementImageType>;
  WriterType::Pointer displacementWriter = WriterType::New();
  displacementWriter->SetFileName(argv[3]);
  displacementWriter->SetInput(multiResRegistrationMethod->GetOutput());
  try
  {
    displacementWriter->Update();
  }
  catch (itk::ExceptionObject & ex)
  {
    std::cerr << "Exception caught!" << std::endl;
    std::cerr << ex << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
