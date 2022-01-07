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

#include <algorithm>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVectorImage.h"
#include "itkTestingMacros.h"
#include "itkSimpleFilterWatcher.h"
#include "itkSpectra1DAveragingImageFilter.h"
#include "itkMath.h"

// this does not work for VectorImage, only for Image<Vector> due to:
// error C2440: 'static_cast': cannot convert from 'float *' to 'itk::VariableLengthVector<float> *'
template <typename TInImage1, typename TInImage2>
void
compareBuffers(TInImage1 * baseline, TInImage2 * test)
{
  itk::IdentifierType pixelCount = baseline->GetLargestPossibleRegion().GetNumberOfPixels();
  itk::IdentifierType testPixelCount = test->GetLargestPossibleRegion().GetNumberOfPixels();
  if (pixelCount != testPixelCount)
  {
    itkGenericExceptionMacro(<< "Images have different number of pixels, " << pixelCount << " vs " << testPixelCount);
  }

  itk::IdentifierType differentPixels = 0;
  itk::IdentifierType linePixels = baseline->GetLargestPossibleRegion().GetSize()[0];
  itk::IdentifierType lineCount = pixelCount / linePixels;

  using InputPixelType1 = typename TInImage1::PixelType;
  using InputPixelType2 = typename TInImage2::PixelType;
  auto * bufferBaseline = static_cast<InputPixelType1 *>(baseline->GetBufferPointer());
  auto * bufferTest = static_cast<InputPixelType2 *>(test->GetBufferPointer());
  for (itk::IdentifierType p = 0; p < pixelCount; ++p)
  {
    for (unsigned d = 0; d < TInImage1::PixelType::Dimension; ++d)
    {
      if (!itk::Math::FloatAlmostEqual(bufferBaseline[p][d], bufferTest[p][d], lineCount, 1e-4))
      {
        ++differentPixels;
        std::cout << "Images differ at " << p << "[" << d << "]: " << bufferBaseline[p][d] << " vs " << bufferTest[p][d]
                  << "\n";
      }
    }
  }

  itk::IdentifierType maxDifferentPixels = 0;
  if (differentPixels > maxDifferentPixels)
  {
    itkGenericExceptionMacro(<< differentPixels << " pixels are different, maximum allowed is " << maxDifferentPixels);
  }
}

template <typename TInImage, typename TOutImage>
int
doTest(int argc, char * argv[], itk::SmartPointer<TOutImage> & output)
{
  using AverageFilterType = itk::Spectra1DAveragingImageFilter<TInImage, TOutImage>;
  typename AverageFilterType::Pointer averageFilter = AverageFilterType::New();

  itk::SimpleFilterWatcher watcher(averageFilter);

  for (int i = 2; i < argc; ++i)
  {
    using ReaderType = itk::ImageFileReader<TInImage>;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(argv[i]);
    averageFilter->SetInput(i - 2, reader->GetOutput());

    // it should not be necessary to do this here, the filter is supposed to call this when it needs this input
    ITK_TRY_EXPECT_NO_EXCEPTION(reader->Update());
  }

  ITK_TRY_EXPECT_NO_EXCEPTION(averageFilter->Update());
  output = averageFilter->GetOutput();
  return EXIT_SUCCESS;
}

template <typename TInImage, typename TOutImage>
itk::SmartPointer<TOutImage>
createReference(int argc, char * argv[])
{
  using AverageFilterType = itk::Spectra1DAveragingImageFilter<TInImage, TOutImage>;
  auto averageFilter = AverageFilterType::New();

  for (int i = 2; i < argc; ++i)
  {
    using ReaderType = itk::ImageFileReader<TInImage>;
    typename ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(argv[i]);
    averageFilter->SetInput(i - 2, reader->GetOutput());

    // it should not be necessary to do this here, the filter is supposed to call this when it needs this input
    reader->Update();
  }

  averageFilter->Update();
  return averageFilter->GetOutput();
}


int
itkSpectra1DAveragingImageFilterTest(int argc, char * argv[])
{
  if (argc < 3)
  {
    std::cerr << "Usage: " << itkNameOfTestExecutableMacro(argv);
    std::cerr << " outputImage inputImage [inputImage ...]";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }
  const char * inFileName = argv[2];

  const unsigned Channels = 31;

  int returnCount = 0;


  using IV3f = itk::Image<itk::Vector<float, Channels>, 3>;
  using IV2f = itk::Image<itk::Vector<float, Channels>, 2>;
  using IV1f = itk::Image<itk::Vector<float, Channels>, 1>;
  using VI3f = itk::VectorImage<float, 3>;
  using VI2f = itk::VectorImage<float, 2>;
  using VI1f = itk::VectorImage<float, 1>;

  IV2f::Pointer testIV2;
  IV1f::Pointer testIV1;
  VI2f::Pointer testVI2;
  VI1f::Pointer testVI1;

  // create reference result - all other images will be compared to this one
  IV2f::Pointer reference = createReference<IV2f, IV2f>(argc, argv);

  // now test different dimensionality
  returnCount += doTest<VI3f, VI2f>(argc, argv, testVI2);
  returnCount += doTest<VI3f, VI1f>(argc, argv, testVI1);

  // test all 4 combinations (VI x IV) with same dimensionality
  returnCount += doTest<VI2f, VI2f>(argc, argv, testVI2);
  // compareBuffers(reference.GetPointer(), testVI2.GetPointer()); // does not compile
  returnCount += doTest<VI2f, IV2f>(argc, argv, testIV2);
  compareBuffers(reference.GetPointer(), testIV2.GetPointer());
  returnCount += doTest<IV2f, VI2f>(argc, argv, testVI2);
  returnCount += doTest<IV2f, IV2f>(argc, argv, testIV2);
  compareBuffers(reference.GetPointer(), testIV2.GetPointer());

  // all (VI x IV) combinations with dimensionality 2 -> 1
  returnCount += doTest<VI2f, VI1f>(argc, argv, testVI1);
  returnCount += doTest<VI2f, IV1f>(argc, argv, testIV1);
  compareBuffers(reference.GetPointer(), testIV1.GetPointer());
  returnCount += doTest<IV2f, VI1f>(argc, argv, testVI1);
  returnCount += doTest<IV2f, IV1f>(argc, argv, testIV1);
  compareBuffers(reference.GetPointer(), testIV1.GetPointer());

  // test float/double
  returnCount += doTest<itk::VectorImage<double, 2>, IV1f>(argc, argv, testIV1);

  // the next line requires ITK 5.3 or later
  returnCount += doTest<itk::Image<itk::Vector<double, Channels>, 2>, IV1f>(argc, argv, testIV1);
  compareBuffers(reference.GetPointer(), testIV1.GetPointer());

  ITK_TRY_EXPECT_NO_EXCEPTION(itk::WriteImage(reference, argv[1], true));

  if (returnCount != 12 * EXIT_SUCCESS)
  {
    std::cerr << "Test failed";
    return EXIT_FAILURE;
  }

  std::cout << "Test finished";
  return EXIT_SUCCESS;
}
