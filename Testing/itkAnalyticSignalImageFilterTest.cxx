/** 
 * @file itkAnalyticSignalImageFilterTest.cxx
 * @author Matthew McCormick (thewtex) <matt@mmmccormick.com>
 * @date 2010-02-11
 */

#include "itkTestMain.h"

void RegisterTests()
{
  REGISTER_TEST( itkAnalyticSignalImageFilterTest );
}

#include <complex>
#include <iostream>
#include <string>

#include "itkComplexToImaginaryImageFilter.h"
#include "itkComplexToRealImageFilter.h"
#include "itkConstantPadImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkAnalyticSignalImageFilter.h"

int itkAnalyticSignalImageFilterTest( int argc, char* argv[] )
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
  typedef itk::AnalyticSignalImageFilter< PixelType, Dimension > AnalyticType;
  typedef itk::ComplexToRealImageFilter< ComplexImageType, ImageType > RealFilterType;
  typedef itk::ComplexToImaginaryImageFilter< ComplexImageType, ImageType > ImaginaryFilterType;
  typedef itk::ImageFileWriter< ImageType > WriterType;

  ReaderType::Pointer reader = ReaderType::New();
  PadType::Pointer    pad    = PadType::New();
  AnalyticType::Pointer    analytic    = AnalyticType::New();
  RealFilterType::Pointer realFilter = RealFilterType::New();
  ImaginaryFilterType::Pointer imaginaryFilter = ImaginaryFilterType::New();
  WriterType::Pointer writer = WriterType::New();

  reader->SetFileName( argv[1] );
  pad->SetInput( reader->GetOutput() );
  ImageType::SizeType padSize;
  padSize[0] = 0;
  padSize[1] = 28;
  pad->SetPadUpperBound( padSize );
  pad->SetConstant( 0. );
  analytic->SetInput( pad->GetOutput() );
  analytic->SetDirection( 1 );
  realFilter->SetInput( analytic->GetOutput() );
  imaginaryFilter->SetInput( analytic->GetOutput() );

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

  analytic.Print( std::cout );

  return EXIT_SUCCESS;
}
