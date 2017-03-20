#include "itkTestMain.h"

void RegisterTests()
{
  REGISTER_TEST( itkBlockMatchingBayesianRegularizationDisplacementCalculatorTest );
}

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVector.h"

#include "itkBlockMatchingImageRegistrationMethod.h"
#include "itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter.h"
#include "itkBlockMatchingSearchRegionImageInitializer.h"
#include "itkBlockMatchingBayesianRegularizationDisplacementCalculator.h"
#include "itkSplitComponentsImageFilter.h"

int itkBlockMatchingBayesianRegularizationDisplacementCalculatorTest( int argc, char* argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " fixedImage movingImage displacementImagePrefix";
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

  // Input file readers.
  typedef itk::ImageFileReader< InputImageType > ReaderType;
  ReaderType::Pointer fixedReader = ReaderType::New();
  fixedReader->SetFileName( argv[1] );
  ReaderType::Pointer movingReader = ReaderType::New();
  movingReader->SetFileName( argv[2] );

  // Make the search region image.
  typedef itk::BlockMatching::SearchRegionImageInitializer< InputImageType, InputImageType > SearchRegionInitializerType;
  SearchRegionInitializerType::Pointer searchRegions = SearchRegionInitializerType::New();
  searchRegions->SetFixedImage( fixedReader->GetOutput() );
  searchRegions->SetMovingImage( movingReader->GetOutput() );
  RadiusType blockRadius;
  blockRadius[0] = 20;
  blockRadius[1] = 4;
  RadiusType searchRadius;
  searchRadius[0] = 130;
  searchRadius[1] = 6;
  searchRegions->SetFixedBlockRadius( blockRadius );
  searchRegions->SetSearchRegionRadius( searchRadius );
  // For speed...
  searchRegions->SetOverlap( 3.0 );

  // The image registration method.
  typedef itk::BlockMatching::ImageRegistrationMethod< InputImageType, InputImageType, MetricImageType, DisplacementImageType, CoordRepType > RegistrationMethodType;
  RegistrationMethodType::Pointer registrationMethod = RegistrationMethodType::New();
  registrationMethod->SetFixedImage( fixedReader->GetOutput() );
  registrationMethod->SetMovingImage( movingReader->GetOutput() );
  registrationMethod->SetInput( searchRegions->GetOutput() );
  registrationMethod->SetRadius( blockRadius );

  // Our similarity metric.
  typedef itk::BlockMatching::NormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter< InputImageType, InputImageType, MetricImageType > MetricImageFilterType;
  MetricImageFilterType::Pointer metricImageFilter = MetricImageFilterType::New();

  registrationMethod->SetMetricImageFilter( metricImageFilter );

  // Perform regularization.
  typedef itk::BlockMatching::BayesianRegularizationDisplacementCalculator<
    MetricImageType, DisplacementImageType > DisplacmentRegularizerType;
  DisplacmentRegularizerType::Pointer regularizer = DisplacmentRegularizerType::New();
  regularizer->SetMetricLowerBound( -1.0 );
  MetricImageType::SpacingType strainSigma;
  strainSigma[0] = 0.08;
  strainSigma[1] = 0.04;
  regularizer->SetStrainSigma( strainSigma );
  regularizer->SetMaximumIterations( 3 );
  registrationMethod->SetMetricImageToDisplacementCalculator( regularizer );

  // Break the displacement vector image into components.
  typedef itk::SplitComponentsImageFilter< DisplacementImageType,
          MetricImageType > TensorComponentsFilterType;
  TensorComponentsFilterType::Pointer componentsFilter = TensorComponentsFilterType::New();
  componentsFilter->SetInput( registrationMethod->GetOutput() );

  typedef itk::ImageFileWriter< MetricImageType > WriterType;
  WriterType::Pointer displacementWriter = WriterType::New();
  try
    {
    displacementWriter->SetInput( componentsFilter->GetOutput( 0 ) );
    displacementWriter->SetFileName( std::string( argv[3] ) + "Component0.mha" );
    displacementWriter->Update();
    displacementWriter->SetInput( componentsFilter->GetOutput( 1 ) );
    displacementWriter->SetFileName( std::string( argv[3] ) + "Component1.mha" );
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
