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

#include "itkTimeGainCompensationImageFilter.h"
#include "itkCurvilinearArraySpecialCoordinatesImage.h"
#include "itkBModeImageFilter.h"
#include "itkCastImageFilter.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"

int itkTimeGainCompensationImageFilterTest( int argc, char * argv[] )
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

  const unsigned int Dimension = 2;
  typedef signed short                                                                IntegerPixelType;
  typedef itk::CurvilinearArraySpecialCoordinatesImage< IntegerPixelType, Dimension > IntegerImageType;
  typedef float                                                                       RealPixelType;
  typedef itk::CurvilinearArraySpecialCoordinatesImage< RealPixelType, Dimension >    RealImageType;
  typedef itk::Image< RealPixelType, Dimension >                                      ScanConvertedImageType;

  typedef itk::ImageFileReader< IntegerImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImageFileName );

  typedef itk::TimeGainCompensationImageFilter< IntegerImageType > TGCFilterType;
  TGCFilterType::Pointer tgcFilter = TGCFilterType::New();
  tgcFilter->SetInput( reader->GetOutput() );

  typedef itk::CastImageFilter< IntegerImageType, RealImageType > CasterType;
  CasterType::Pointer caster = CasterType::New();
  caster->SetInput( tgcFilter->GetOutput() );

  typedef itk::BModeImageFilter< RealImageType, RealImageType > BModeFilterType;
  BModeFilterType::Pointer bmodeFilter = BModeFilterType::New();
  bmodeFilter->SetInput( caster->GetOutput() );

  try
    {
    bmodeFilter->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
    }

  tgcFilter->Print( std::cout );

  RealImageType::Pointer curvilinearArrayImage = bmodeFilter->GetOutput();
  curvilinearArrayImage->DisconnectPipeline();
  curvilinearArrayImage->Print( std::cout );
  std::cout << "Size: " << curvilinearArrayImage->GetPixelContainer()->Size() << std::endl;
  const RealImageType::SizeType inputSize = curvilinearArrayImage->GetLargestPossibleRegion().GetSize();
  const double lateralAngularSeparation = (vnl_math::pi / 2.0 + 0.5) /
    (inputSize[1] - 1);
  curvilinearArrayImage->SetLateralAngularSeparation( lateralAngularSeparation );
  const double radiusStart = 7.0;
  const double radiusStop = 112.5;
  curvilinearArrayImage->SetFirstSampleDistance( radiusStart );
  curvilinearArrayImage->SetRadiusSampleSize( (radiusStop - radiusStart) / (inputSize[0] -1) );

  typedef itk::ResampleImageFilter< RealImageType, ScanConvertedImageType > ResamplerType;
  ResamplerType::Pointer resampler = ResamplerType::New();
  resampler->SetInput( curvilinearArrayImage );
  RealImageType::SizeType outputSize;
  outputSize[0] = 800;
  outputSize[1] = 800;
  resampler->SetSize( outputSize );
  RealImageType::SpacingType outputSpacing;
  outputSpacing.Fill( 0.15 );
  resampler->SetOutputSpacing( outputSpacing );
  RealImageType::PointType outputOrigin;
  outputOrigin[0] = outputSize[0] * outputSpacing[0] / -2.0;
  outputOrigin[1] = radiusStart * std::cos( vnl_math::pi / 4.0 );
  resampler->SetOutputOrigin( outputOrigin );

  typedef itk::ImageFileWriter< ScanConvertedImageType > WriterType;
  //typedef itk::ImageFileWriter< RealImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputImageFileName );
  writer->SetInput( resampler->GetOutput() );
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
