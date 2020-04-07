/** 
 * @file itkOpenCL1DComplexToComplexImageFilterTest.cxx
 * @author Matthew McCormick (thewtex) <matt@mmmccormick.com>
 * @date 2010-02-15
 */

#include <complex>
#include <iostream>
#include <string>

#include "itkComplexToRealImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkRealAndImaginaryToComplexImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkOpenCLComplexToComplex1DFFTImageFilter.h"

int itkOpenCLComplexToComplex1DFFTImageFilterTest( int argc, char* argv[] )
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
  typedef itk::OpenCL1DComplexToComplexImageFilter< PixelType, Dimension > FFTType;
  typedef itk::RealAndImaginaryToComplexImageFilter< PixelType, PixelType, PixelType, Dimension > JoinFilterType;
  typedef itk::ComplexToRealImageFilter< ComplexImageType, ImageType > ToRealFilterType;
  typedef itk::ExtractImageFilter< ImageType, ImageType > ExtractType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  ReaderType::Pointer readerReal = ReaderType::New();
  ReaderType::Pointer readerImag = ReaderType::New();
  FFTType::Pointer    fft    = FFTType::New();
  JoinFilterType::Pointer joinFilter = JoinFilterType::New();
  ToRealFilterType::Pointer toReal = ToRealFilterType::New();
  ExtractType::Pointer extractor = ExtractType::New();
  WriterType::Pointer writer = WriterType::New();

  readerReal->SetFileName( std::string( argv[1] ) + "Real128.mhd" );
  readerImag->SetFileName( std::string( argv[1] ) + "Imaginary128.mhd" );
  joinFilter->SetInput1( readerReal->GetOutput() );
  joinFilter->SetInput2( readerImag->GetOutput() );
  fft->SetTransformDirection( FFTType::INVERSE );
  fft->SetInput( joinFilter->GetOutput() );
  fft->SetDirection( 0 );
  toReal->SetInput( fft->GetOutput() );
  extractor->SetInput( toReal->GetOutput() );
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
