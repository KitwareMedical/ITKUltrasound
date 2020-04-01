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

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVectorImage.h"
#include "itkTestingMacros.h"

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
  const char * referenceSpectraImageFileName = argv[2];
  const char * outputImageFileName = argv[3];

  const unsigned int Dimension = 2;
  typedef short                              PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;

  typedef itk::ImageFileReader< ImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImageFileName );
  ITK_TRY_EXPECT_NO_EXCEPTION( reader->UpdateLargestPossibleRegion() );
  ImageType::ConstPointer rfImage = reader->GetOutput();

  ImageType::Pointer sideLines = ImageType::New();
  sideLines->CopyInformation( rfImage );
  sideLines->SetRegions( rfImage->GetLargestPossibleRegion() );
  sideLines->Allocate();
  sideLines->FillBuffer( 5 );

  typedef itk::Spectra1DSupportWindowImageFilter< ImageType > SpectraSupportWindowFilterType;
  SpectraSupportWindowFilterType::Pointer spectraSupportWindowFilter = SpectraSupportWindowFilterType::New();
  spectraSupportWindowFilter->SetInput( sideLines );
  spectraSupportWindowFilter->SetFFT1DSize( 128 );
  spectraSupportWindowFilter->SetStep( 16 );
  ITK_TRY_EXPECT_NO_EXCEPTION( spectraSupportWindowFilter->UpdateLargestPossibleRegion() );

  typedef SpectraSupportWindowFilterType::OutputImageType     SupportWindowImageType;
  SupportWindowImageType * supportWindowImage = spectraSupportWindowFilter->GetOutput();

  typedef float                                               SpectraComponentType;
  typedef itk::VectorImage< SpectraComponentType, Dimension > SpectraImageType;

  typedef itk::ImageFileReader< SpectraImageType > ReferenceSpectraReaderType;
  ReferenceSpectraReaderType::Pointer referenceSpectraReader = ReferenceSpectraReaderType::New();
  referenceSpectraReader->SetFileName( referenceSpectraImageFileName );
  ITK_TRY_EXPECT_NO_EXCEPTION( referenceSpectraReader->UpdateLargestPossibleRegion() );
  typedef SpectraImageType::PixelType SpectraPixelType;
  SpectraImageType::IndexType referenceIndex;
  referenceIndex.Fill( 0 );
  const SpectraPixelType referenceSpectra = referenceSpectraReader->GetOutput()->GetPixel( referenceIndex );
  SpectraImageType::Pointer referenceSpectraImage = SpectraImageType::New();
  referenceSpectraImage->CopyInformation( supportWindowImage );
  referenceSpectraImage->SetRegions( supportWindowImage->GetLargestPossibleRegion() );
  referenceSpectraImage->SetNumberOfComponentsPerPixel( referenceSpectraReader->GetOutput()->GetNumberOfComponentsPerPixel() );
  referenceSpectraImage->Allocate();
  referenceSpectraImage->FillBuffer( referenceSpectra );

  typedef itk::Spectra1DImageFilter< ImageType, SupportWindowImageType, SpectraImageType > SpectraFilterType;
  SpectraFilterType::Pointer spectraFilter = SpectraFilterType::New();
  spectraFilter->SetInput( rfImage );
  spectraFilter->SetSupportWindowImage( spectraSupportWindowFilter->GetOutput() );
  spectraFilter->SetReferenceSpectraImage( referenceSpectraImage );
  ITK_TRY_EXPECT_NO_EXCEPTION( spectraFilter->UpdateLargestPossibleRegion() );

  typedef itk::ImageFileWriter< SpectraImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputImageFileName );
  ITK_TRY_EXPECT_NO_EXCEPTION( writer->SetInput( spectraFilter->GetOutput() ) );

  ITK_TRY_EXPECT_NO_EXCEPTION( writer->UpdateLargestPossibleRegion() );

  return EXIT_SUCCESS;
}
