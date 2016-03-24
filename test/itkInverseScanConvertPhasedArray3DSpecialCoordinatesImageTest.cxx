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

#include "itkMath.h"
#include "itkPhasedArray3DSpecialCoordinatesImage.h"
#include "itkResampleImageFilter.h"
#include "itkGaborImageSource.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkMetaDataObject.h"
#include "itkImageFileWriter.h"

int itkInverseScanConvertPhasedArray3DSpecialCoordinatesImageTest( int argc, char* argv[] )
{
  if( argc < 2 )
    {
    std::cerr << "Usage: " << argv[0] << "outputImageFile" << std::endl;
    return EXIT_FAILURE;
    }
  const char * outputImageFile = argv[1];

  const unsigned int Dimension = 3;
  typedef float                                                  RealPixelType;
  typedef unsigned char                                          PixelType;
  typedef itk::Image< RealPixelType, Dimension >                 RealImageType;
  typedef itk::Image< PixelType, Dimension >                     InputImageType;
  typedef itk::PhasedArray3DSpecialCoordinatesImage< PixelType > OutputImageType;
  typedef double                                                 CoordRepType;


  typedef itk::GaborImageSource< RealImageType > SourceType;
  SourceType::Pointer source = SourceType::New();

  SourceType::ArrayType sigma;
  sigma[0] = 5.0;
  sigma[1] = 5.0;
  sigma[2] = 8.0;
  source->SetSigma( sigma );

  SourceType::ArrayType mean;
  mean.Fill( 0.0 );
  source->SetMean( mean );

  RealImageType::SizeType size;
  size.Fill( 128 );
  source->SetSize( size );

  RealImageType::SpacingType spacing;
  spacing[0] = 0.2;
  spacing[1] = 0.2;
  spacing[2] = 0.4;
  source->SetSpacing( spacing );

  RealImageType::PointType origin;
  for( unsigned int ii = 0; ii < Dimension - 1; ++ii )
    {
    origin[ii] = -1 * spacing[ii] * size[ii] / 2;
    }
  origin[2] = 0.0;
  source->SetOrigin( origin );

  source->SetFrequency( 0.2 );

  source->SetCalculateImaginaryPart( true );

  typedef itk::RescaleIntensityImageFilter< RealImageType, InputImageType > RescalerType;
  RescalerType::Pointer rescaler = RescalerType::New();
  rescaler->SetInput( source->GetOutput() );

  try
    {
    rescaler->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error while generating Gabor input: " << error << std::endl;
    return EXIT_FAILURE;
    }
  InputImageType::Pointer gaborInput = rescaler->GetOutput();
  gaborInput->DisconnectPipeline();

  typedef itk::ResampleImageFilter< InputImageType, OutputImageType > ResamplerType;
  ResamplerType::Pointer resampler = ResamplerType::New();
  resampler->SetInput( gaborInput );

  OutputImageType::Pointer phasedArraySampledData = resampler->GetOutput();

  const double azimuthAngularSeparation = 5.0 * vnl_math::pi/180.0;
  const double elevationAngularSeparation = 1.0 * vnl_math::pi/180.0;
  const double radiusSampleSize = 0.2;
  const double firstSampleDistance = 8.0;
  phasedArraySampledData->SetAzimuthAngularSeparation( azimuthAngularSeparation );
  phasedArraySampledData->SetElevationAngularSeparation( elevationAngularSeparation );
  phasedArraySampledData->SetRadiusSampleSize( radiusSampleSize );
  phasedArraySampledData->SetFirstSampleDistance( firstSampleDistance );

  OutputImageType::SizeType outputSize;
  outputSize[0] = 16;
  outputSize[1] = 32;
  outputSize[2] = 64;
  resampler->SetSize( outputSize );

  try
    {
    resampler->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error while resampling data: " << error << std::endl;
    return EXIT_FAILURE;
    }
  phasedArraySampledData->DisconnectPipeline();

  // Since this is non-standard image metadata, make sure it is serialized
  // into the MetaDataDictionary
  itk::MetaDataDictionary & dictionary = phasedArraySampledData->GetMetaDataDictionary();
  itk::EncapsulateMetaData< double >( dictionary, "AzimuthAngularSeparation", azimuthAngularSeparation );
  itk::EncapsulateMetaData< double >( dictionary, "ElevationAngularSeparation", elevationAngularSeparation );
  itk::EncapsulateMetaData< double >( dictionary, "RadiusSampleSize", radiusSampleSize );
  itk::EncapsulateMetaData< double >( dictionary, "FirstSampleDistance", firstSampleDistance );
  phasedArraySampledData->Print( std::cout );
  phasedArraySampledData->GetMetaDataDictionary().Print( std::cout );

  typedef itk::ImageFileWriter< OutputImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputImageFile );
  writer->SetInput( phasedArraySampledData );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error while writing output: " << error << std::endl;
    }

  return EXIT_SUCCESS;
}
