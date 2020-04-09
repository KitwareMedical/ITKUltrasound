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

#include <complex>
#include <iostream>
#include <string>

#include "itkComplexToRealImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkComposeImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkOpenCLComplexToComplex1DFFTImageFilter.h"

int
itkOpenCLComplexToComplex1DFFTImageFilterTest(int argc, char * argv[])
{
  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImagePrefix outputImage";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }

  using PixelType = float;
  const unsigned int Dimension = 2;

  using ImageType = itk::Image<PixelType, Dimension>;
  using ComplexImageType = itk::Image<std::complex<PixelType>, Dimension>;

  using ReaderType = itk::ImageFileReader<ImageType>;
  using FFTType = itk::OpenCLComplexToComplex1DFFTImageFilter<ComplexImageType, ComplexImageType>;
  using JoinFilterType = itk::ComposeImageFilter<ImageType, ComplexImageType>;
  using ToRealFilterType = itk::ComplexToRealImageFilter<ComplexImageType, ImageType>;
  using ExtractType = itk::ExtractImageFilter<ImageType, ImageType>;
  using WriterType = itk::ImageFileWriter<ImageType>;

  ReaderType::Pointer       readerReal = ReaderType::New();
  ReaderType::Pointer       readerImag = ReaderType::New();
  FFTType::Pointer          fft = FFTType::New();
  JoinFilterType::Pointer   joinFilter = JoinFilterType::New();
  ToRealFilterType::Pointer toReal = ToRealFilterType::New();
  ExtractType::Pointer      extractor = ExtractType::New();
  WriterType::Pointer       writer = WriterType::New();

  readerReal->SetFileName(std::string(argv[1]) + "Real128.mhd");
  readerImag->SetFileName(std::string(argv[1]) + "Imaginary128.mhd");
  joinFilter->SetInput1(readerReal->GetOutput());
  joinFilter->SetInput2(readerImag->GetOutput());
  fft->SetTransformDirection(FFTType::INVERSE);
  fft->SetInput(joinFilter->GetOutput());
  fft->SetDirection(0);
  toReal->SetInput(fft->GetOutput());
  extractor->SetInput(toReal->GetOutput());
  ImageType::RegionType            outputRegion;
  ImageType::RegionType::IndexType outputIndex;
  outputIndex[0] = 0;
  outputIndex[1] = 0;
  outputRegion.SetIndex(outputIndex);
  ImageType::RegionType::SizeType outputSize;
  outputSize[0] = 100;
  outputSize[1] = 100;
  outputRegion.SetSize(outputSize);
  extractor->SetExtractionRegion(outputRegion);
  writer->SetInput(extractor->GetOutput());
  writer->SetFileName(argv[2]);

  try
  {
    writer->Update();
  }
  catch (itk::ExceptionObject & excep)
  {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }

  fft.Print(std::cout);

  return EXIT_SUCCESS;
}
