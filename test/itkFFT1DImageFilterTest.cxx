/*=========================================================================
 *
 *  Copyright Insight Software Consortium
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

#include "itkForward1DFFTImageFilter.h"
#include "itkInverse1DFFTImageFilter.h"

int itkFFT1DImageFilterTest( int argc, char* argv[] )
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

  typedef double PixelType;
  const unsigned int Dimension = 2;
  const unsigned int direction = 1;

  typedef itk::Image< PixelType, Dimension >                                         ImageType;
  typedef itk::Image< std::complex< PixelType >, Dimension >                         ComplexImageType;

  typedef itk::ImageFileReader< ImageType >                                          ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImage );

  typedef itk::Forward1DFFTImageFilter< ImageType, ComplexImageType > FFTForwardType;
  FFTForwardType::Pointer fftForward = FFTForwardType::New();
  fftForward->SetInput( reader->GetOutput() );
  fftForward->SetDirection( direction );

  typedef itk::Inverse1DFFTImageFilter< ComplexImageType, ImageType > FFTInverseType;
  FFTInverseType::Pointer fftInverse = FFTInverseType::New();
  fftInverse->SetInput( fftForward->GetOutput() );
  fftInverse->SetDirection( direction );

  typedef itk::ImageFileWriter< ImageType >                                          WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( fftInverse->GetOutput() );
  writer->SetFileName( outputImage );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }

  fftForward.Print( std::cout );
  fftInverse.Print( std::cout );

  return EXIT_SUCCESS;
}
