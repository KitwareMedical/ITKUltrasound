/** 
 * @file itkFFT1DRealToComplexConjugateImageFilterTest.cxx
 * @brief Test for FFT1DRealToComplexConjugateImageFilter
 * @author Matthew McCormick (thewtex) <matt@mmmccormick.com>
 * @date 2010-02-09
 */

#include "itkTestMain.h"

void RegisterTests()
{
  REGISTER_TEST( itkOpenCL1DRealToComplexConjugateImageFilterTest );
}

#include <complex>
#include <iostream>
#include <string>

#include "itkComplexToImaginaryImageFilter.h"
#include "itkComplexToRealImageFilter.h"
#include "itkConstantPadImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkOpenCL1DRealToComplexConjugateImageFilter.h"

int itkOpenCL1DRealToComplexConjugateImageFilterTest( int argc, char* argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage outputImagePrefix";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }

  typedef float PixelType;
  const unsigned int Dimension = 2;

  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::Image< std::complex< PixelType >, Dimension > ComplexImageType;

  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ConstantPadImageFilter< ImageType, ImageType > PadType;
  typedef itk::OpenCL1DRealToComplexConjugateImageFilter< PixelType, Dimension > FFTType;
  typedef itk::ComplexToRealImageFilter< ComplexImageType, ImageType > RealFilterType;
  typedef itk::ComplexToImaginaryImageFilter< ComplexImageType, ImageType > ImaginaryFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  PadType::Pointer    pad    = PadType::New();
  FFTType::Pointer    fft    = FFTType::New();
  RealFilterType::Pointer realFilter = RealFilterType::New();
  ImaginaryFilterType::Pointer imaginaryFilter = ImaginaryFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  reader->SetFileName( argv[1] );
  pad->SetInput( reader->GetOutput() );
  ImageType::SizeType padSize;
  padSize[0] = 28;
  padSize[1] = 0;
  pad->SetPadUpperBound( padSize );
  pad->SetConstant( 0. );
  fft->SetInput( pad->GetOutput() );
  fft->SetDirection( 0 );
  realFilter->SetInput( fft->GetOutput() );
  imaginaryFilter->SetInput( fft->GetOutput() );

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
