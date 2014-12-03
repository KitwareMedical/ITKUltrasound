// Test the performance of the BModeImageFilter
//

#include "itkImageFileReader.h"
#include "itkIntensityWindowingImageFilter.h"
#include "itkTimeProbe.h"
#include "itkImageFileWriter.h"

#include "itkBModeImageFilter.h"

int itkBModeImageFilterTestTiming( int argc, char* argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage outputImage";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }
  const char * inputImageFileName = argv[1];
  const char * outputImageFileName = argv[2];

  typedef float PixelType;
  const unsigned int Dimension = 3;
  typedef itk::Image< PixelType, Dimension > ImageType;

  typedef itk::ImageFileReader               < ImageType            > ReaderType;
  typedef itk::BModeImageFilter              < ImageType, ImageType > BModeFilterType;
  typedef itk::IntensityWindowingImageFilter < ImageType, ImageType > WindowingType;

  ReaderType::Pointer reader     = ReaderType::New();
  BModeFilterType::Pointer bMode = BModeFilterType::New();
  WindowingType::Pointer window  = WindowingType::New();

  reader->SetFileName( inputImageFileName );

  try
    {
    reader->Update();
    }
  catch ( itk::ExceptionObject& e )
    {
    std::cerr << "Error: " << e << std::endl;
    return EXIT_FAILURE;
    }
  ImageType::Pointer input = reader->GetOutput();
  input->DisconnectPipeline();
  bMode->SetInput( input );
  window->SetInput( bMode->GetOutput() );

  itk::TimeProbe clock;

  const int runs = 1000;
  for(int i = 0; i < runs; i++)
    {
    bMode->Modified();
    clock.Start();
    window->Update();
    clock.Stop();
    }

  ImageType::SizeType size = bMode->GetOutput()->GetLargestPossibleRegion().GetSize();

  double frame_rate = static_cast< double >( size[2] ) / clock.GetMean();

  std::cout << "Frame rate achieved over " << clock.GetNumberOfStarts() << " runs was " << frame_rate << " fp" << clock.GetUnit() << "." << std::endl;

  typedef itk::ImageFileWriter< ImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputImageFileName );
  writer->SetInput( window->GetOutput() );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
