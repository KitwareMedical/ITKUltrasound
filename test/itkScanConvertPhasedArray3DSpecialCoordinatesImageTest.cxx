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
#include "itkRescaleIntensityImageFilter.h"
#include "itkMetaDataObject.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

/** MetaImage reads in all non-standard MetaDataDictionary entries with the
 * std::string type. */
double GetMetaDataAsDouble(const itk::MetaDataDictionary & dictionary, const std::string & key)
{
  std::string valueAsStr;
  if( !itk::ExposeMetaData< std::string >( dictionary, key, valueAsStr ) )
    {
    std::cerr << "Error while reading: " << key << std::endl;
    throw std::runtime_error("Could not find MetaDataDictionary entry");
    }
  std::istringstream istrm( valueAsStr );
  double value;
  istrm >> value;
  return value;
}

int itkScanConvertPhasedArray3DSpecialCoordinatesImageTest( int argc, char* argv[] )
{
  if( argc < 2 )
    {
    std::cerr << "Usage: " << argv[0] << "inputImageFile outputImageFile" << std::endl;
    return EXIT_FAILURE;
    }
  const char * inputImageFile = argv[1];
  const char * outputImageFile = argv[2];

  const unsigned int Dimension = 3;
  typedef unsigned char                                          PixelType;
  typedef itk::PhasedArray3DSpecialCoordinatesImage< PixelType > InputImageType;
  typedef itk::Image< PixelType, Dimension >                     OutputImageType;
  typedef double                                                 CoordRepType;


  typedef itk::ImageFileReader< InputImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImageFile );

  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error while reading input: " << error << std::endl;
    return EXIT_FAILURE;
    }
  InputImageType::Pointer inputImage = reader->GetOutput();
  inputImage->DisconnectPipeline();

  // Load the metadata-from the MetaDataDictionary
  // into the MetaDataDictionary
  const itk::MetaDataDictionary & dictionary = inputImage->GetMetaDataDictionary();
  try
    {
    inputImage->SetAzimuthAngularSeparation( GetMetaDataAsDouble( dictionary, "AzimuthAngularSeparation" ) );
    inputImage->SetElevationAngularSeparation( GetMetaDataAsDouble( dictionary, "ElevationAngularSeparation" ) );
    inputImage->SetRadiusSampleSize( GetMetaDataAsDouble( dictionary, "RadiusSampleSize" ) );
    inputImage->SetFirstSampleDistance( GetMetaDataAsDouble( dictionary, "FirstSampleDistance" ) );
    }
  catch( std::exception & error )
    {
    std::cerr << "Error: " << error.what() << std::endl;
    return EXIT_FAILURE;
    }
  inputImage->Print( std::cout );


  typedef itk::ResampleImageFilter< InputImageType, OutputImageType > ResamplerType;
  ResamplerType::Pointer resampler = ResamplerType::New();
  resampler->SetInput( inputImage );

  OutputImageType::SizeType size;
  size.Fill( 128 );
  resampler->SetSize( size );

  OutputImageType::SpacingType spacing;
  spacing[0] = 0.2;
  spacing[1] = 0.2;
  spacing[2] = 0.2;
  resampler->SetOutputSpacing( spacing );

  OutputImageType::PointType origin;
  for( unsigned int ii = 0; ii < Dimension - 1; ++ii )
    {
    origin[ii] = -1 * spacing[ii] * size[ii]  / 2;
    }
  origin[2] = 0.0;
  resampler->SetOutputOrigin( origin );

  typedef itk::ImageFileWriter< OutputImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputImageFile );
  writer->SetInput( resampler->GetOutput() );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error while resampling data: " << error << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
