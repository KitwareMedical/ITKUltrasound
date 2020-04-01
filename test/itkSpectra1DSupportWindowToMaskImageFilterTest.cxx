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

#include "itkSpectra1DSupportWindowImageFilter.h"
#include "itkSpectra1DSupportWindowToMaskImageFilter.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkPermuteAxesImageFilter.h"

int itkSpectra1DSupportWindowToMaskImageFilterTest( int argc, char* argv[] )
{
  const unsigned int Dimension = 2;
  typedef short                              PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef ImageType::IndexType               IndexType;

  if( argc < 4 )
    {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage outputImage index0 index1";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }
  const char * inputImageFileName = argv[1];
  const char * outputImageFileName = argv[2];
  IndexType windowIndex;
  windowIndex[0] = atoi( argv[3] );
  windowIndex[1] = atoi( argv[4] );

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

  typedef SpectraSupportWindowFilterType::OutputImageType                                       SupportWindowImageType;
  typedef itk::Image< unsigned char, Dimension >                                                MaskImageType;
  typedef itk::Spectra1DSupportWindowToMaskImageFilter< SupportWindowImageType, MaskImageType > SpectraSupportWindowMaskFilterType;
  SpectraSupportWindowMaskFilterType::Pointer spectraSupportWindowMaskFilter = SpectraSupportWindowMaskFilterType::New();
  spectraSupportWindowMaskFilter->SetInput( spectraSupportWindowFilter->GetOutput() );
  spectraSupportWindowMaskFilter->SetMaskIndex( windowIndex );

  typedef itk::ImageFileWriter< MaskImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( spectraSupportWindowMaskFilter->GetOutput() );
  writer->SetFileName( outputImageFileName );

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
