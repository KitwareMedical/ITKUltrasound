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

#include "itkSpectra1DSupportWindowImageFilter.h"
#include "itkSpectra1DSupportWindowToMaskImageFilter.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkPermuteAxesImageFilter.h"

int
itkSpectra1DSupportWindowToMaskImageFilterTest(int argc, char * argv[])
{
  const unsigned int Dimension = 2;
  using PixelType = short;
  using ImageType = itk::Image<PixelType, Dimension>;
  using IndexType = ImageType::IndexType;

  if (argc < 4)
  {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage outputImage index0 index1";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }
  const char * inputImageFileName = argv[1];
  const char * outputImageFileName = argv[2];
  IndexType    windowIndex;
  windowIndex[0] = atoi(argv[3]);
  windowIndex[1] = atoi(argv[4]);

  using ReaderType = itk::ImageFileReader<ImageType>;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(inputImageFileName);

  // Want RF to be along direction 0
  using PermuteAxesFilterType = itk::PermuteAxesImageFilter<ImageType>;
  PermuteAxesFilterType::Pointer               permuteAxesFilter = PermuteAxesFilterType::New();
  PermuteAxesFilterType::PermuteOrderArrayType permuteOrder;
  permuteOrder[0] = 1;
  permuteOrder[1] = 0;
  permuteAxesFilter->SetOrder(permuteOrder);
  permuteAxesFilter->SetInput(reader->GetOutput());
  try
  {
    permuteAxesFilter->UpdateOutputInformation();
  }
  catch (itk::ExceptionObject & error)
  {
    std::cout << "Error: " << error << std::endl;
    return EXIT_FAILURE;
  }
  ImageType::ConstPointer rfImage = permuteAxesFilter->GetOutput();

  ImageType::Pointer sideLines = ImageType::New();
  sideLines->CopyInformation(rfImage);
  sideLines->SetRegions(rfImage->GetLargestPossibleRegion());
  sideLines->Allocate();
  sideLines->FillBuffer(5);

  using SpectraSupportWindowFilterType = itk::Spectra1DSupportWindowImageFilter<ImageType>;
  SpectraSupportWindowFilterType::Pointer spectraSupportWindowFilter = SpectraSupportWindowFilterType::New();
  spectraSupportWindowFilter->SetInput(sideLines);

  using SupportWindowImageType = SpectraSupportWindowFilterType::OutputImageType;
  using MaskImageType = itk::Image<unsigned char, Dimension>;
  using SpectraSupportWindowMaskFilterType =
    itk::Spectra1DSupportWindowToMaskImageFilter<SupportWindowImageType, MaskImageType>;
  SpectraSupportWindowMaskFilterType::Pointer spectraSupportWindowMaskFilter =
    SpectraSupportWindowMaskFilterType::New();
  spectraSupportWindowMaskFilter->SetInput(spectraSupportWindowFilter->GetOutput());
  spectraSupportWindowMaskFilter->SetMaskIndex(windowIndex);

  using WriterType = itk::ImageFileWriter<MaskImageType>;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(spectraSupportWindowMaskFilter->GetOutput());
  writer->SetFileName(outputImageFileName);

  try
  {
    writer->Update();
  }
  catch (itk::ExceptionObject & error)
  {
    std::cout << "Error: " << error << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
