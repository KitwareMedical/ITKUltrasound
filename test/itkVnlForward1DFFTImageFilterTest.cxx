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

#include "itkComplexToImaginaryImageFilter.h"
#include "itkComplexToRealImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkVnlForward1DFFTImageFilter.h"

int itkVnlForward1DFFTImageFilterTest( int argc, char* argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage outputImagePrefix";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }

  typedef double PixelType;
  const unsigned int Dimension = 2;

  typedef itk::Image< PixelType, Dimension >                                            ImageType;
  typedef itk::Image< std::complex< PixelType >, Dimension >                            ComplexImageType;

  typedef itk::ImageFileReader< ImageType >                                             ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );

  typedef itk::VnlForward1DFFTImageFilter< ImageType, ComplexImageType > FFTType;
  FFTType::Pointer    fft    = FFTType::New();
  fft->SetInput( reader->GetOutput() );

  typedef itk::ComplexToRealImageFilter< ComplexImageType, ImageType >                  RealFilterType;
  RealFilterType::Pointer realFilter = RealFilterType::New();
  realFilter->SetInput( fft->GetOutput() );

  typedef itk::ComplexToImaginaryImageFilter< ComplexImageType, ImageType >             ImaginaryFilterType;
  ImaginaryFilterType::Pointer imaginaryFilter = ImaginaryFilterType::New();
  imaginaryFilter->SetInput( fft->GetOutput() );

  typedef itk::ImageFileWriter< ImageType >                                             WriterType;
  WriterType::Pointer writer = WriterType::New();

  try
    {
    writer->SetInput( realFilter->GetOutput() );
    writer->SetFileName( std::string( argv[2] ) + "Real.mha" );
    writer->Update();

    writer->SetInput( imaginaryFilter->GetOutput() );
    writer->SetFileName( std::string( argv[2] ) + "Imaginary.mha" );
    writer->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }

  fft.Print( std::cout );

  return EXIT_SUCCESS;
}
