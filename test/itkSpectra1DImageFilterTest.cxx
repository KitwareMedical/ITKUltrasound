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

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkPermuteAxesImageFilter.h"
#include "itkVectorImage.h"

#include "itkSpectra1DSupportWindowImageFilter.h"
#include "itkSpectra1DImageFilter.h"

int itkSpectra1DImageFilterTest( int argc, char* argv[] )
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
  typedef short PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;

  typedef itk::ImageFileReader< ImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImageFileName );

  // Want RF to be along direction 0
  typedef itk::PermuteAxesImageFilter< ImageType > PermuteAxesFilterType;
  PermuteAxesFilterType::Pointer permuteAxesFilter = PermuteAxesFilterType::New();
  PermuteAxesFilterType::PermuteOrderArrayType permuteOrder;
  permuteOrder[0] = 1;
  permuteOrder[1] = 0;
  permuteAxesFilter->SetOrder( permuteOrder );
  permuteAxesFilter->SetInput( reader->GetOutput() );
  try
    {
    permuteAxesFilter->UpdateOutputInformation();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cout << "Error: " << error << std::endl;
    return EXIT_FAILURE;
    }
  ImageType::ConstPointer rfImage = permuteAxesFilter->GetOutput();

  ImageType::Pointer sideLines = ImageType::New();
  sideLines->CopyInformation( rfImage );
  sideLines->SetRegions( rfImage->GetLargestPossibleRegion() );
  sideLines->Allocate();
  sideLines->FillBuffer( 5 );

  typedef itk::Spectra1DSupportWindowImageFilter< ImageType > SpectraSupportWindowFilterType;
  SpectraSupportWindowFilterType::Pointer spectraSupportWindowFilter = SpectraSupportWindowFilterType::New();
  spectraSupportWindowFilter->SetInput( sideLines );

  typedef SpectraSupportWindowFilterType::OutputImageType SupportWindowImageType;
  typedef float SpectraComponentType;
  typedef itk::VectorImage< SpectraComponentType, Dimension > SpectraImageType;
  typedef itk::Spectra1DImageFilter< ImageType, SupportWindowImageType, SpectraImageType > SpectraFilterType;
  SpectraFilterType::Pointer spectraFilter = SpectraFilterType::New();
  spectraFilter->SetInputImage( rfImage );
  spectraFilter->SetSupportWindowImage( spectraSupportWindowFilter->GetOutput() );

  typedef itk::ImageFileWriter< SpectraImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputImageFileName );
  writer->SetInput( spectraFilter->GetOutput() );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & error )
    {
    std::cout << "Error: " << error << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
