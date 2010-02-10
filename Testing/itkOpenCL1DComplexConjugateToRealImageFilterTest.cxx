/** 
 * @file itkFFT1DComplexConjugateToRealImageFilterTest.cxx
 * @brief Test for FFT1DComplexConjugateToRealImageFilter
 * @author Matthew McCormick (thewtex) <matt@mmmccormick.com>
 * @date 2010-02-10
 */

#include "itkTestMain.h"

void RegisterTests()
{
  REGISTER_TEST( itkOpenCL1DComplexConjugateToRealImageFilterTest );
}

#include <complex>
#include <iostream>
#include <string>

#include "itkExtractImageFilter.h"
#include "itkRealAndImaginaryToComplexImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkOpenCL1DComplexConjugateToRealImageFilter.h"

int itkOpenCL1DComplexConjugateToRealImageFilterTest( int argc, char* argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImagePrefix outputImage";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }

  typedef float PixelType;
  const unsigned int Dimension = 2;

  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::Image< std::complex< PixelType >, Dimension > ComplexImageType;

  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::OpenCL1DComplexConjugateToRealImageFilter< PixelType, Dimension > FFTType;
  typedef itk::RealAndImaginaryToComplexImageFilter< PixelType, PixelType, PixelType, Dimension > JoinFilterType;
  typedef itk::ExtractImageFilter< ImageType, ImageType > ExtractType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  ReaderType::Pointer readerReal = ReaderType::New();
  ReaderType::Pointer readerImag = ReaderType::New();
  FFTType::Pointer    fft    = FFTType::New();
  JoinFilterType::Pointer joinFilter = JoinFilterType::New();
  ExtractType::Pointer extractor = ExtractType::New();
  WriterType::Pointer writer = WriterType::New();

  readerReal->SetFileName( std::string( argv[1] ) + "Real128.mhd" );
  readerImag->SetFileName( std::string( argv[1] ) + "Imaginary128.mhd" );
  joinFilter->SetInput1( readerReal->GetOutput() );
  joinFilter->SetInput2( readerImag->GetOutput() );
  fft->SetInput( joinFilter->GetOutput() );
  fft->SetDirection( 0 );
  extractor->SetInput( fft->GetOutput() );
  ImageType::RegionType outputRegion;
  ImageType::RegionType::IndexType outputIndex;
  outputIndex[0] = 0;
  outputIndex[1] = 0;
  outputRegion.SetIndex( outputIndex );
  ImageType::RegionType::SizeType outputSize;
  outputSize[0] = 100;
  outputSize[1] = 100;
  outputRegion.SetSize( outputSize );
  extractor->SetExtractionRegion( outputRegion );
  writer->SetInput( extractor->GetOutput() );
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

  fft.Print( std::cout );

  return EXIT_SUCCESS;
}
