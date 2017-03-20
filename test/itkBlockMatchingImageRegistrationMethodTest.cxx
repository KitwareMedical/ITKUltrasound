#include "itkTestMain.h"

void RegisterTests()
{
  REGISTER_TEST( itkBlockMatchingImageRegistrationMethodTest );
}

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVector.h"

#include "itkBlockMatchingImageRegistrationMethod.h"
#include "itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter.h"
#include "itkBlockMatchingSearchRegionImageInitializer.h"

int itkBlockMatchingImageRegistrationMethodTest( int argc, char* argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixedImage movingImage displacementImage";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }

  const unsigned int Dimension = 2;
  typedef signed short InputPixelType;
  typedef itk::Image< InputPixelType, Dimension > InputImageType;
  typedef InputImageType::SizeType RadiusType;

  typedef double MetricPixelType;
  typedef itk::Image< MetricPixelType, Dimension > MetricImageType;

  typedef itk::Vector< MetricPixelType, Dimension > VectorType;
  typedef itk::Image< VectorType, Dimension > DisplacementImageType;

  typedef double CoordRepType;

  typedef itk::ImageFileReader< InputImageType > ReaderType;
  ReaderType::Pointer fixedReader = ReaderType::New();
  fixedReader->SetFileName( argv[1] );
  ReaderType::Pointer movingReader = ReaderType::New();
  movingReader->SetFileName( argv[2] );

  typedef itk::BlockMatching::SearchRegionImageInitializer< InputImageType, InputImageType > SearchRegionInitializerType;
  SearchRegionInitializerType::Pointer searchRegions = SearchRegionInitializerType::New();
  searchRegions->SetFixedImage( fixedReader->GetOutput() );
  searchRegions->SetMovingImage( movingReader->GetOutput() );
  RadiusType blockRadius;
  blockRadius[0] = 20;
  blockRadius[1] = 4;
  RadiusType searchRadius;
  searchRadius[0] = 130;
  searchRadius[1] = 5;
  searchRegions->SetFixedBlockRadius( blockRadius );
  searchRegions->SetSearchRegionRadius( searchRadius );

  typedef itk::BlockMatching::ImageRegistrationMethod< InputImageType, InputImageType, MetricImageType, DisplacementImageType, CoordRepType > RegistrationMethodType;
  RegistrationMethodType::Pointer registrationMethod = RegistrationMethodType::New();
  registrationMethod->SetFixedImage( fixedReader->GetOutput() );
  registrationMethod->SetMovingImage( movingReader->GetOutput() );
  registrationMethod->SetInput( searchRegions->GetOutput() );
  registrationMethod->SetRadius( blockRadius );

  typedef itk::BlockMatching::NormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter< InputImageType, InputImageType, MetricImageType > MetricImageFilterType;
  MetricImageFilterType::Pointer metricImageFilter = MetricImageFilterType::New();

  registrationMethod->SetMetricImageFilter( metricImageFilter );

  typedef itk::ImageFileWriter< DisplacementImageType > WriterType;
  WriterType::Pointer displacementWriter = WriterType::New();
  displacementWriter->SetFileName( argv[3] );
  displacementWriter->SetInput( registrationMethod->GetOutput() );
  try
    {
    displacementWriter->Update();
    }
  catch (itk::ExceptionObject& ex)
    {
      std::cerr << "Exception caught!" << std::endl;
      std::cerr << ex << std::endl;
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
