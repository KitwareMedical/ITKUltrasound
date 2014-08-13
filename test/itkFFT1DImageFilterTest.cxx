/** 
 * @file itkFFT1DImageFilterTest.cxx
 * @brief Test for FFT1DImageFilter
 * @author Matthew McCormick (thewtex) <matt@mmmccormick.com>
 * @date 2009-12-06
 */

#include <complex>
#include <iostream>
#include <string>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkFFT1DRealToComplexConjugateImageFilter.h"
#include "itkFFT1DComplexConjugateToRealImageFilter.h"

int itkFFT1DImageFilterTest( int argc, char* argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage outputImage";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }

  typedef double PixelType;
  const unsigned int Dimension = 2;
  const unsigned int direction = 1;

  typedef itk::Image< PixelType, Dimension >                                  ImageType;
  typedef itk::Image< std::complex< PixelType >, Dimension >                  ComplexImageType;

  typedef itk::ImageFileReader< ImageType >                                   ReaderType;
  typedef itk::FFT1DRealToComplexConjugateImageFilter< PixelType, Dimension > FFTForwardType;
  typedef itk::FFT1DComplexConjugateToRealImageFilter< PixelType, Dimension > FFTInverseType;
  typedef itk::ImageFileWriter< ImageType >                                   WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  FFTForwardType::Pointer fftForward = FFTForwardType::New();
  FFTInverseType::Pointer fftInverse = FFTInverseType::New();
  WriterType::Pointer writer = WriterType::New();

  reader->SetFileName( argv[1] );
  fftForward->SetInput( reader->GetOutput() );
  fftForward->SetDirection( direction );
  fftInverse->SetInput( fftForward->GetOutput() );
  fftInverse->SetDirection( direction );
  writer->SetInput( fftInverse->GetOutput() );
  writer->SetFileName( argv[2] );

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
