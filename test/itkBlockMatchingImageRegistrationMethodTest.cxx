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
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVector.h"

#include "itkBlockMatchingImageRegistrationMethod.h"
#include "itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter.h"
#include "itkBlockMatchingSearchRegionImageInitializer.h"

int
itkBlockMatchingImageRegistrationMethodTest(int argc, char * argv[])
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

  using SearchRegionInitializerType = itk::BlockMatching::SearchRegionImageInitializer<InputImageType, InputImageType>;
  SearchRegionInitializerType::Pointer searchRegions = SearchRegionInitializerType::New();
  searchRegions->SetFixedImage(fixedReader->GetOutput());
  searchRegions->SetMovingImage(movingReader->GetOutput());
  RadiusType blockRadius;
  blockRadius[0] = 20;
  blockRadius[1] = 4;
  RadiusType searchRadius;
  searchRadius[0] = 130;
  searchRadius[1] = 5;
  searchRegions->SetFixedBlockRadius(blockRadius);
  searchRegions->SetSearchRegionRadius(searchRadius);

  using RegistrationMethodType = itk::BlockMatching::
    ImageRegistrationMethod<InputImageType, InputImageType, MetricImageType, DisplacementImageType, CoordRepType>;
  RegistrationMethodType::Pointer registrationMethod = RegistrationMethodType::New();
  registrationMethod->SetFixedImage(fixedReader->GetOutput());
  registrationMethod->SetMovingImage(movingReader->GetOutput());
  registrationMethod->SetInput(searchRegions->GetOutput());
  registrationMethod->SetRadius(blockRadius);

  using MetricImageFilterType = itk::BlockMatching::
    NormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter<InputImageType, InputImageType, MetricImageType>;
  MetricImageFilterType::Pointer metricImageFilter = MetricImageFilterType::New();

  registrationMethod->SetMetricImageFilter(metricImageFilter);

  using WriterType = itk::ImageFileWriter<DisplacementImageType>;
  WriterType::Pointer displacementWriter = WriterType::New();
  displacementWriter->SetFileName(argv[3]);
  displacementWriter->SetInput(registrationMethod->GetOutput());
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
