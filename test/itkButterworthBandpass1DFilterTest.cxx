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

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTestingMacros.h"

#include "itkForward1DFFTImageFilter.h"
#include "itkInverse1DFFTImageFilter.h"
#include "itkFrequencyDomain1DImageFilter.h"
#include "itkButterworthBandpass1DFilterFunction.h"

int itkButterworthBandpass1DFilterTest( int argc, char* argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage outputImage";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }
  const char * inputImage = argv[1];
  const char * outputImage = argv[2];

  using PixelType = float;
  const unsigned int Dimension = 2;
  const unsigned int direction = 0;

  using ImageType = itk::Image< PixelType, Dimension >;
  using ComplexImageType = itk::Image< std::complex< PixelType >, Dimension >;

  using ReaderType = itk::ImageFileReader< ImageType >;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImage );

  using FFTForwardType = itk::Forward1DFFTImageFilter< ImageType, ComplexImageType >;
  FFTForwardType::Pointer fftForward = FFTForwardType::New();
  fftForward->SetInput( reader->GetOutput() );
  fftForward->SetDirection( direction );

  using FilterFunctionType = itk::ButterworthBandpass1DFilterFunction;
  FilterFunctionType::Pointer filterFunction = FilterFunctionType::New();
  filterFunction->SetLowerFrequency( 0.12 );
  filterFunction->SetOrder( 7 );
  filterFunction->Print( std::cout );

  using FrequencyFilterType = itk::FrequencyDomain1DImageFilter< ComplexImageType, ComplexImageType >;
  FrequencyFilterType::Pointer frequencyFilter = FrequencyFilterType::New();
  frequencyFilter->SetInput( fftForward->GetOutput() );
  frequencyFilter->SetDirection( direction );
  frequencyFilter->SetFilterFunction( filterFunction.GetPointer() );
  frequencyFilter->Print( std::cout );

  using FFTInverseType = itk::Inverse1DFFTImageFilter< ComplexImageType, ImageType >;
  FFTInverseType::Pointer fftInverse = FFTInverseType::New();
  fftInverse->SetInput( frequencyFilter->GetOutput() );
  fftInverse->SetDirection( direction );

  using WriterType = itk::ImageFileWriter< ImageType >;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( fftInverse->GetOutput() );
  writer->SetFileName( outputImage );

  ITK_TRY_EXPECT_NO_EXCEPTION( writer->Update() );

  return EXIT_SUCCESS;
}
