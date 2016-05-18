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

#include "itkCurvilinearArraySpecialCoordinatesImage.h"
#include "itkSliceSeriesSpecialCoordinatesImage.h"

#include "itkUltrasoundImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkResampleImageFilter.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkEuler3DTransform.h"
#include "vnl/vnl_math.h"


int itkSliceSeriesSpecialCoordinatesImageTest( int argc, char * argv[] )
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

  const unsigned int Dimension = 3;
  const unsigned int SliceDimension = Dimension - 1;

  typedef unsigned char PixelType;
  typedef double        ParametersValueType;

  typedef itk::Image< PixelType, Dimension >                                        ImageType;
  typedef itk::CurvilinearArraySpecialCoordinatesImage< PixelType, SliceDimension > SliceImageType;
  typedef itk::Euler3DTransform< ParametersValueType >                              TransformType;
  typedef itk::SliceSeriesSpecialCoordinatesImage< SliceImageType, TransformType >  SliceSeriesImageType;


  typedef itk::UltrasoundImageFileReader < SliceImageType > SliceReaderType;
  SliceReaderType::Pointer sliceReader = SliceReaderType::New();
  sliceReader->SetFileName( inputImageFileName );
  try
    {
    sliceReader->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
    }
  SliceImageType::Pointer sliceImage = sliceReader->GetOutput();
  sliceImage->Print( std::cout );

  typedef itk::ImageFileReader< SliceSeriesImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImageFileName );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
    }
  SliceSeriesImageType::Pointer image = reader->GetOutput();
  image->SetSliceImage( sliceImage );

  double elevationAngle = 0.0;
  const itk::SizeValueType numberOfSlices = image->GetLargestPossibleRegion().GetSize()[Dimension - 1];
  for( itk::SizeValueType sliceIndex = 0; sliceIndex < numberOfSlices; ++sliceIndex )
    {
    TransformType::Pointer transform = TransformType::New();
    transform->SetRotation( 0.0, 0.0, sliceIndex * 5.0 * vnl_math::pi / 180.0 );
    image->SetSliceTransform( sliceIndex, transform );
    }
  TransformType::ConstPointer obtainedTransform = image->GetSliceTransform( 0 );

  image->Print( std::cout );

  typedef itk::ResampleImageFilter< SliceSeriesImageType, ImageType > ResamplerType;
  ResamplerType::Pointer resampler = ResamplerType::New();
  resampler->SetInput( reader->GetOutput() );
  //SpecialCoordinatesImageType::SizeType outputSize = reader->GetOutput()->GetLargestPossibleRegion().GetSize();
  //outputSize[0] = 800;
  //outputSize[1] = 800;
  //resampler->SetSize( outputSize );
  //SpecialCoordinatesImageType::SpacingType outputSpacing;
  //outputSpacing.Fill( 0.15 );
  //resampler->SetOutputSpacing( outputSpacing );
  //SpecialCoordinatesImageType::PointType outputOrigin;
  //outputOrigin[0] = outputSize[0] * outputSpacing[0] / -2.0;
  //outputOrigin[1] = radiusStart * std::cos( vnl_math::pi / 4.0 );
  //outputOrigin[2] = reader->GetOutput()->GetOrigin()[2];
  //resampler->SetOutputOrigin( outputOrigin );

  //typedef itk::WindowedSincInterpolateImageFunction< SpecialCoordinatesImageType, 3 > WindowedSincInterpolatorType;
  //WindowedSincInterpolatorType::Pointer sincInterpolator = WindowedSincInterpolatorType::New();
  //sincInterpolator->SetInputImage( curvilinearArrayImage );

  //typedef itk::ImageFileWriter< ImageType > WriterType;
  //WriterType::Pointer writer = WriterType::New();
  //writer->SetFileName( outputImageFileName );
  //writer->SetInput( resampler->GetOutput() );
  //try
    //{
    //writer->Update();
    //}
  //catch( itk::ExceptionObject & error )
    //{
    //std::cerr << "Error: " << error << std::endl;
    //return EXIT_FAILURE;
    //}
  return EXIT_SUCCESS;
}
