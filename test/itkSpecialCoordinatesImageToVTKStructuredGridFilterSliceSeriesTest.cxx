/*=========================================================================
 *
 *  Copyright NumFOCUS
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

#include "itkUltrasoundImageFileReader.h"
#include "itkHDF5UltrasoundImageIO.h"
#include "itkHDF5UltrasoundImageIOFactory.h"
#include "itkSpecialCoordinatesImageToVTKStructuredGridFilter.h"

#include "itkTestingMacros.h"
#include "itkSliceSeriesSpecialCoordinatesImage.h"
#include "itkEuler3DTransform.h"

#include "vtkStructuredGridWriter.h"
#include "vtkNew.h"

int itkSpecialCoordinatesImageToVTKStructuredGridFilterSliceSeriesTest( int argc, char * argv [] )
{
  if( argc < 3 )
    {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage outputStructuredGridFileName";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }
  const char * inputImageFileName = argv[1];
  const char * outputStructuredGridFileName = argv[2];

  itk::HDF5UltrasoundImageIOFactory::RegisterOneFactory();

  const unsigned int Dimension = 3;
  using PixelType = unsigned char;
  using ParametersValueType = double;

  using SliceImageType = itk::Image< PixelType, Dimension - 1 >;
  using TransformType = itk::Euler3DTransform< ParametersValueType >;
  using SpecialCoordinatesImageType = itk::SliceSeriesSpecialCoordinatesImage< SliceImageType, TransformType, PixelType, Dimension >;

  using ReaderType = itk::UltrasoundImageFileReader< SpecialCoordinatesImageType >;
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

  using ConversionFilterType = itk::SpecialCoordinatesImageToVTKStructuredGridFilter< SpecialCoordinatesImageType >;
  ConversionFilterType::Pointer conversionFilter = ConversionFilterType::New();

  conversionFilter->SetInput( reader->GetOutput() );

  try
    {
    conversionFilter->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cerr << "Error: " << error << std::endl;
    return EXIT_FAILURE;
    }

  vtkSmartPointer< vtkStructuredGrid > structuredGrid = conversionFilter->GetOutput();
  structuredGrid->Print( std::cout );

  vtkNew< vtkStructuredGridWriter > writer;
  writer->SetInputData( structuredGrid );
  writer->SetFileName( outputStructuredGridFileName );
  writer->SetFileTypeToBinary();
  writer->Write();

  return EXIT_SUCCESS;
}
