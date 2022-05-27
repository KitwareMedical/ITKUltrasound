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
// Test the performance of the BModeImageFilter
//

#include "itkImageFileReader.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkTimeProbe.h"
#include "itkImageFileWriter.h"

#include "itkBModeImageFilter.h"

int
itkBModeImageFilterTestTiming(int argc, char * argv[])
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

  using PixelType = float;
  const unsigned int Dimension = 3;
  using ImageType = itk::Image<PixelType, Dimension>;

  using ReaderType = itk::ImageFileReader<ImageType>;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inputImageFileName);
  try
  {
    reader->Update();
  }
  catch (itk::ExceptionObject & e)
  {
    std::cerr << "Error: " << e << std::endl;
    return EXIT_FAILURE;
  }
  ImageType::Pointer input = reader->GetOutput();
  input->DisconnectPipeline();

  using BModeFilterType = itk::BModeImageFilter<ImageType, ImageType>;
  BModeFilterType::Pointer bMode = BModeFilterType::New();
  bMode->SetInput(input);

  using WindowingType = itk::IntensityWindowingImageFilter<ImageType, ImageType>;
  WindowingType::Pointer window = WindowingType::New();
  window->SetInput(bMode->GetOutput());

  itk::TimeProbe clock;

  const int runs = 1000;
  for (int i = 0; i < runs; i++)
  {
    bMode->Modified();
    clock.Start();
    window->Update();
    clock.Stop();
  }

  ImageType::SizeType size = bMode->GetOutput()->GetLargestPossibleRegion().GetSize();

  double frame_rate = static_cast<double>(size[2]) / clock.GetMean();

  std::cout << "Frame rate achieved over " << clock.GetNumberOfStarts() << " runs was " << frame_rate << " fp"
            << clock.GetUnit() << "." << std::endl;

  using WriterType = itk::ImageFileWriter<ImageType>;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(outputImageFileName);
  writer->SetInput(window->GetOutput());
  try
  {
    writer->Update();
  }
  catch (itk::ExceptionObject & error)
  {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
