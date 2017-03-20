#include "itkTestMain.h"

void RegisterTests()
{
  REGISTER_TEST( itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilterTest );
}

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter.h"

int itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilterTest( int argc, char* argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputFixedImage inputMovingImage metricImage ";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 2;
  typedef signed short InputPixelType;
  typedef itk::Image< InputPixelType, Dimension > InputImageType;

  typedef double MetricPixelType;
  typedef itk::Image< MetricPixelType, Dimension > MetricImageType;

  typedef itk::ImageFileReader< InputImageType > ReaderType;
  typedef itk::BlockMatching::NormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter< InputImageType, InputImageType, MetricImageType >  FilterType;
  typedef itk::ImageFileWriter< MetricImageType >  WriteType;

  ReaderType::Pointer readerFixed  = ReaderType::New();
  ReaderType::Pointer readerMoving = ReaderType::New();
  FilterType::Pointer filter = FilterType::New();
  WriteType::Pointer  writer = WriteType::New();

  readerFixed->SetFileName( argv[1] );
  readerMoving->SetFileName( argv[2] );

  filter->SetFixedImage( readerFixed->GetOutput() );
  filter->SetMovingImage( readerMoving->GetOutput() );
  typedef MetricImageType::RegionType RegionType;
  RegionType fixedRegion;
  RegionType::SizeType fixedSize;
  fixedSize[0] = 31;
  fixedSize[1] = 5;
  fixedRegion.SetSize( fixedSize );
  RegionType::IndexType fixedIndex;
  fixedIndex[0] = 999;
  fixedIndex[1] = 99;
  fixedRegion.SetIndex( fixedIndex );
  filter->SetFixedImageRegion( fixedRegion );
  RegionType movingRegion;
  RegionType::SizeType movingSize;
  movingSize[0] = 100;
  movingSize[1] = 30;
  movingRegion.SetSize( movingSize );
  RegionType::IndexType movingIndex;
  movingIndex[0] = fixedIndex[0] - movingSize[0] / 2;
  movingIndex[1] = fixedIndex[1] - movingSize[1] / 2;
  movingRegion.SetIndex( movingIndex );
  filter->SetMovingImageRegion( movingRegion );

  writer->SetInput( filter->GetOutput() );
  writer->SetFileName( argv[3] );

  try
    {
    writer->Update();
    }
  catch (itk::ExceptionObject& ex)
    {
      std::cerr << "Exception caught!" << std::endl;
      std::cerr << ex << std::endl;
      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
