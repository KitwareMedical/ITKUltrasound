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

#include "itkSpectra1DSupportWindowImageFilter.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkPermuteAxesImageFilter.h"
#include "itkTestingMacros.h"

int itkSpectra1DSupportWindowImageFilterTest( int argc, char* argv[] )
{
  if( argc < 2 )
    {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage";
    std::cerr << std::endl;
    return EXIT_FAILURE;
    }
  const char * inputImageFileName = argv[1];

  const unsigned int Dimension = 2;
  typedef short                              PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;

  typedef itk::ImageFileReader< ImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImageFileName );
  ITK_TRY_EXPECT_NO_EXCEPTION( reader->UpdateLargestPossibleRegion() );

  // Want RF to be along direction 0
  typedef itk::PermuteAxesImageFilter< ImageType > PermuteAxesFilterType;
  PermuteAxesFilterType::Pointer permuteAxesFilter = PermuteAxesFilterType::New();
  PermuteAxesFilterType::PermuteOrderArrayType permuteOrder;
  permuteOrder[0] = 1;
  permuteOrder[1] = 0;
  permuteAxesFilter->SetOrder( permuteOrder );
  permuteAxesFilter->SetInput( reader->GetOutput() );
  ITK_TRY_EXPECT_NO_EXCEPTION( permuteAxesFilter->UpdateOutputInformation() );
  ImageType::ConstPointer rfImage = permuteAxesFilter->GetOutput();

  ImageType::Pointer sideLines = ImageType::New();
  sideLines->CopyInformation( rfImage );
  sideLines->SetRegions( rfImage->GetLargestPossibleRegion() );
  sideLines->Allocate();
  sideLines->FillBuffer( 5 );

  typedef itk::Spectra1DSupportWindowImageFilter< ImageType > SpectraSupportWindowFilterType;
  SpectraSupportWindowFilterType::Pointer spectraSupportWindowFilter = SpectraSupportWindowFilterType::New();
  spectraSupportWindowFilter->SetInput( sideLines );
  ITK_TRY_EXPECT_NO_EXCEPTION( spectraSupportWindowFilter->Update() );

  spectraSupportWindowFilter->Print( std::cout );
  spectraSupportWindowFilter->GetOutput()->GetMetaDataDictionary().Print(std::cout);

  spectraSupportWindowFilter->SetStep( 10 );
  ITK_TRY_EXPECT_NO_EXCEPTION( spectraSupportWindowFilter->UpdateLargestPossibleRegion() );
  std::cout << "\n\nAfter setting the Step to 10: " << std::endl;
  spectraSupportWindowFilter->Print( std::cout );

  return EXIT_SUCCESS;
}
