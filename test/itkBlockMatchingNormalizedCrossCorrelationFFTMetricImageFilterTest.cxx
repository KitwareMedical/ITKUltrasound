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

#include "itkBlockMatchingNormalizedCrossCorrelationFFTMetricImageFilter.h"

int
itkBlockMatchingNormalizedCrossCorrelationFFTMetricImageFilterTest(int argc, char * argv[])
{
  if (argc < 4)
  {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputFixedImage inputMovingImage metricImage ";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }

  const unsigned int Dimension = 2;
  using InputPixelType = signed short;
  using InputImageType = itk::Image<InputPixelType, Dimension>;

  using MetricPixelType = double;
  using MetricImageType = itk::Image<MetricPixelType, Dimension>;

  using ReaderType = itk::ImageFileReader<InputImageType>;
  using FilterType =
    itk::BlockMatching::NormalizedCrossCorrelationFFTMetricImageFilter<InputImageType, InputImageType, MetricImageType>;
  using WriteType = itk::ImageFileWriter<MetricImageType>;

  ReaderType::Pointer readerFixed = ReaderType::New();
  ReaderType::Pointer readerMoving = ReaderType::New();
  FilterType::Pointer filter = FilterType::New();
  WriteType::Pointer  writer = WriteType::New();

  readerFixed->SetFileName(argv[1]);
  readerMoving->SetFileName(argv[2]);

  filter->SetFixedImage(readerFixed->GetOutput());
  filter->SetMovingImage(readerMoving->GetOutput());
  using RegionType = MetricImageType::RegionType;
  RegionType           fixedRegion;
  RegionType::SizeType fixedSize;
  fixedSize[0] = 31;
  fixedSize[1] = 5;
  fixedRegion.SetSize(fixedSize);
  RegionType::IndexType fixedIndex;
  fixedIndex[0] = 999;
  fixedIndex[1] = 99;
  fixedRegion.SetIndex(fixedIndex);
  filter->SetFixedImageRegion(fixedRegion);
  RegionType           movingRegion;
  RegionType::SizeType movingSize;
  movingSize[0] = 100;
  movingSize[1] = 30;
  movingRegion.SetSize(movingSize);
  RegionType::IndexType movingIndex;
  movingIndex[0] = fixedIndex[0] - movingSize[0] / 2;
  movingIndex[1] = fixedIndex[1] - movingSize[1] / 2;
  movingRegion.SetIndex(movingIndex);
  filter->SetMovingImageRegion(movingRegion);

  writer->SetInput(filter->GetOutput());
  writer->SetFileName(argv[3]);

  try
  {
    writer->Update();
  }
  catch (itk::ExceptionObject & ex)
  {
    std::cerr << "Exception caught!" << std::endl;
    std::cerr << ex << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
