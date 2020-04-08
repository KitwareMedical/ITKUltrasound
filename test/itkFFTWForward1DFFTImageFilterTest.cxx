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
#include <string>

#include "itkComplexToImaginaryImageFilter.h"
#include "itkComplexToRealImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkFFTWForward1DFFTImageFilter.h"

int
itkFFTWForward1DFFTImageFilterTest(int argc, char * argv[])
{
  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage outputImagePrefix";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }

  using PixelType = double;
  const unsigned int Dimension = 2;

  using ImageType = itk::Image<PixelType, Dimension>;
  using ComplexImageType = itk::Image<std::complex<PixelType>, Dimension>;

  using ReaderType = itk::ImageFileReader<ImageType>;
  using FFTType = itk::FFTWForward1DFFTImageFilter<ImageType, ComplexImageType>;
  using RealFilterType = itk::ComplexToRealImageFilter<ComplexImageType, ImageType>;
  using ImaginaryFilterType = itk::ComplexToImaginaryImageFilter<ComplexImageType, ImageType>;
  using WriterType = itk::ImageFileWriter<ImageType>;

  ReaderType::Pointer          reader = ReaderType::New();
  FFTType::Pointer             fft = FFTType::New();
  RealFilterType::Pointer      realFilter = RealFilterType::New();
  ImaginaryFilterType::Pointer imaginaryFilter = ImaginaryFilterType::New();
  WriterType::Pointer          writer = WriterType::New();

  reader->SetFileName(argv[1]);
  fft->SetInput(reader->GetOutput());
  realFilter->SetInput(fft->GetOutput());
  imaginaryFilter->SetInput(fft->GetOutput());

  try
  {
    writer->SetInput(realFilter->GetOutput());
    writer->SetFileName(std::string(argv[2]) + "Real.mha");
    writer->Update();

    writer->SetInput(imaginaryFilter->GetOutput());
    writer->SetFileName(std::string(argv[2]) + "Imaginary.mha");
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
