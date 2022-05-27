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

#include "itkImageFileWriter.h"
#include "itkTestingMacros.h"

#include "itkReplaceNonFiniteImageFilter.h"

#include <limits>

int
itkReplaceNonFiniteImageFilterTest(int argc, char * argv[])
{
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " outputImage";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }
  const char * outputImageFileName = argv[1];

  const unsigned int Dimension = 2;
  using PixelType = float;

  using ImageType = itk::Image<PixelType, Dimension>;

  ImageType::RegionType region;
  ImageType::SizeType   size;
  size.Fill(10);
  region.SetSize(size);

  ImageType::Pointer image = ImageType::New();
  image->SetRegions(region);
  image->Allocate();
  image->FillBuffer(7.5f);

  ImageType::IndexType index;
  index.Fill(3);
  image->SetPixel(index, std::numeric_limits<PixelType>::infinity());
  index.Fill(4);
  image->SetPixel(index, -std::numeric_limits<PixelType>::infinity());
  index.Fill(5);
  image->SetPixel(index, std::numeric_limits<PixelType>::quiet_NaN());

  using FilterType = itk::ReplaceNonFiniteImageFilter<ImageType>;
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput(image);

  using WriterType = itk::ImageFileWriter<ImageType>;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(outputImageFileName);
  writer->SetInput(filter->GetOutput());
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
