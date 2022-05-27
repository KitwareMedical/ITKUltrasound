/*=========================================================================
 *
 *  Copyright NumFOCUS
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         https://www.apache.org/licenses/LICENSE-2.0.txt
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
#include "itkSpectra1DNormalizeImageFilter.h"
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
      if (!itk::Math::AlmostEquals(bufferBaseline[p][d], bufferTest[p][d]))
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

template <typename TInImage, typename TReferenceImage>
int
doTest(const char *                  inFileName,
       const char *                  refFileName,
       const char *                  outFileName,
       itk::SmartPointer<TInImage> & output)
{
  using AverageFilterType = itk::Spectra1DNormalizeImageFilter<TInImage, TReferenceImage>;
  auto averageFilter = AverageFilterType::New();

  itk::SimpleFilterWatcher watcher(averageFilter);

  auto inImage = itk::ReadImage<TInImage>(inFileName);
  auto refImage = itk::ReadImage<TReferenceImage>(refFileName);

  averageFilter->SetInput(0, inImage);
  averageFilter->SetReferenceImage(refImage);

  ITK_TRY_EXPECT_NO_EXCEPTION(averageFilter->Update());
  output = averageFilter->GetOutput();
  return EXIT_SUCCESS;
}

template <typename TInImage, typename TReferenceImage>
itk::SmartPointer<TInImage>
createReference(const char * inFileName, const char * refFileName, const char * outFileName)
{
  using AverageFilterType = itk::Spectra1DNormalizeImageFilter<TInImage, TReferenceImage>;
  auto averageFilter = AverageFilterType::New();

  auto       inImage = itk::ReadImage<TInImage>(inFileName);
  const auto refImage = itk::ReadImage<TReferenceImage>(refFileName);

  averageFilter->SetInput(inImage);
  averageFilter->SetInput("ReferenceImage", refImage.GetPointer());

  averageFilter->Update();
  return averageFilter->GetOutput();
}


int
itkSpectra1DNormalizeImageFilterTest(int argc, char * argv[])
{
  if (argc < 4)
  {
    std::cerr << "Usage: " << itkNameOfTestExecutableMacro(argv);
    std::cerr << " inputImage referenceImage outputImage";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }
  const char * inFileName = argv[1];
  const char * refFileName = argv[2];
  const char * outFileName = argv[3];

  const unsigned Channels = 31;

  int returnCount = 0;

  using VI2d = itk::VectorImage<double, 2>;
  using IV2d = itk::Image<itk::Vector<double, Channels>, 2>;

  using IV2f = itk::Image<itk::Vector<float, Channels>, 2>;
  using IV1f = itk::Image<itk::Vector<float, Channels>, 1>;
  using VI3f = itk::VectorImage<float, 3>;
  using VI2f = itk::VectorImage<float, 2>;
  using VI1f = itk::VectorImage<float, 1>;

  typename IV2f::Pointer testIV2;
  typename VI3f::Pointer testVI3;
  typename VI2f::Pointer testVI2;
  typename VI2d::Pointer testVI2d;
  typename IV2d::Pointer testIV2d;

  // create reference result - all other images will be compared to this one
  typename IV2f::Pointer reference = createReference<IV2f, IV1f>(inFileName, refFileName, outFileName);

  // now test different dimensionality
  returnCount += doTest<VI3f, VI2f>(inFileName, refFileName, outFileName, testVI3);
  returnCount += doTest<VI3f, VI1f>(inFileName, refFileName, outFileName, testVI3);

  // test all 4 combinations (VI x IV) with same dimensionality
  returnCount += doTest<VI2f, VI2f>(inFileName, refFileName, outFileName, testVI2);
  // compareBuffers(reference.GetPointer(), testVI2.GetPointer()); // does not compile
  returnCount += doTest<VI2f, IV2f>(inFileName, refFileName, outFileName, testVI2);
  // compareBuffers(reference.GetPointer(), testVI2.GetPointer()); // does not compile
  returnCount += doTest<IV2f, VI2f>(inFileName, refFileName, outFileName, testIV2);
  compareBuffers(reference.GetPointer(), testIV2.GetPointer());
  returnCount += doTest<IV2f, IV2f>(inFileName, refFileName, outFileName, testIV2);
  compareBuffers(reference.GetPointer(), testIV2.GetPointer());

  // all (VI x IV) combinations with dimensionality 2 -> 1
  returnCount += doTest<VI2f, VI1f>(inFileName, refFileName, outFileName, testVI2);
  ITK_TRY_EXPECT_NO_EXCEPTION(itk::WriteImage(testVI2, outFileName, false));
  returnCount += doTest<VI2f, IV1f>(inFileName, refFileName, outFileName, testVI2);
  // compareBuffers(reference.GetPointer(), testVI2.GetPointer()); // does not compile
  returnCount += doTest<IV2f, VI1f>(inFileName, refFileName, outFileName, testIV2);
  compareBuffers(reference.GetPointer(), testIV2.GetPointer());
  returnCount += doTest<IV2f, IV1f>(inFileName, refFileName, outFileName, testIV2);
  compareBuffers(reference.GetPointer(), testIV2.GetPointer());

  // test float/double
  returnCount += doTest<VI2d, IV1f>(inFileName, refFileName, outFileName, testVI2d);
  // compareBuffers(reference.GetPointer(), testVI2d.GetPointer()); // does not compile

  // the next line requires ITK 5.3 or later
  returnCount += doTest<IV2d, IV1f>(inFileName, refFileName, outFileName, testIV2d);
  compareBuffers(reference.GetPointer(), testIV2d.GetPointer());

  ITK_TRY_EXPECT_NO_EXCEPTION(itk::WriteImage(reference, outFileName, false));

  if (returnCount != 12 * EXIT_SUCCESS)
  {
    std::cerr << "Test failed";
    return EXIT_FAILURE;
  }

  std::cout << "Test finished";
  return EXIT_SUCCESS;
}
