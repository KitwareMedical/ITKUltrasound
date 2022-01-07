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
#ifndef itkSpectra1DAveragingImageFilter_hxx
#define itkSpectra1DAveragingImageFilter_hxx

#include "itkImageScanlineIterator.h"
#include "itkDivideImageFilter.h"
#include "itkTotalProgressReporter.h"

#include <iterator>

namespace itk
{
template <typename TInputImage, typename TOutputImage>
void
Spectra1DAveragingImageFilter<TInputImage, TOutputImage>::GenerateOutputInformation()
{
  auto *            input = const_cast<InputImageType *>(this->GetInput());
  OutputImageType * output = this->GetOutput();

  input->UpdateOutputInformation();

  constexpr unsigned commonDimension = std::min(InputImageType::ImageDimension, OutputImageType::ImageDimension);

  typename OutputImageType::SpacingType outSpacing{ 1.0 }; // 1.0 along all dimensions
  for (unsigned d = 0; d < commonDimension; ++d)
  {
    outSpacing[d] = input->GetSpacing()[d]; // keep as much spacing information as we can
  }

  typename OutputImageType::PointType outOrigin{ 0.0 }; // 0.0 along all dimensions
  for (unsigned d = 0; d < commonDimension; ++d)
  {
    outOrigin[d] = input->GetOrigin()[d]; // keep as much origin information as we can
  }

  // Copying part of the direction matrix is harder and usually not needed. If needed,
  // we could follow the logic from Modules/Core/Common/include/itkExtractImageFilter.hxx

  using OutputRegionType = typename OutputImageType::RegionType;
  using OutputSizeType = typename OutputRegionType::SizeType;

  OutputRegionType outRegion{ OutputSizeType::Filled(1) };            // 1-sized along all dimensions
  outRegion.SetSize(0, input->GetLargestPossibleRegion().GetSize(0)); // but keep the depth dimension

  output->SetSpacing(outSpacing);
  output->SetOrigin(outOrigin);
  output->SetRegions(outRegion);

  this->PrepareOutput(input, output);
}


template <typename TScalarA, typename TScalarB, unsigned VDimension>
Vector<TScalarA, VDimension>
operator+(const Vector<TScalarA, VDimension> & a, const VariableLengthVector<TScalarB> & b)
{
  Vector<TScalarA, VDimension> result{ a };
  for (unsigned i = 0; i < VDimension; ++i)
  {
    result[i] += b[i];
  }
  return result;
}

template <typename TScalarA, typename TScalarB, unsigned VDimension>
VariableLengthVector<TScalarA>
operator+(const VariableLengthVector<TScalarA> & a, const Vector<TScalarB, VDimension> & b)
{
  VariableLengthVector<TScalarA> result{ a };
  for (unsigned i = 0; i < VDimension; ++i)
  {
    result[i] += b[i];
  }
  return result;
}

template <typename TScalarA, typename TScalarB, unsigned VDimension>
Vector<TScalarA, VDimension>&
operator+=(Vector<TScalarA, VDimension> & a, const VariableLengthVector<TScalarB> & b)
{
  for (unsigned i = 0; i < VDimension; ++i)
  {
    a[i] += b[i];
  }
  return a;
}

template <typename TScalarA, typename TScalarB, unsigned VDimension>
VariableLengthVector<TScalarA>&
operator+=(VariableLengthVector<TScalarA> & a, const Vector<TScalarB, VDimension> & b)
{
  for (unsigned i = 0; i < VDimension; ++i)
  {
    a[i] += b[i];
  }
  return a;
}


template <typename TInputImage, typename TOutputImage>
void
Spectra1DAveragingImageFilter<TInputImage, TOutputImage>::GenerateData()
{
  this->UpdateProgress(0.0f);
  auto *            input = const_cast<InputImageType *>(this->GetInput());
  OutputImageType * output = this->GetOutput();

  using InputRegion = typename InputImageType::RegionType;
  InputRegion    inRegion = input->GetLargestPossibleRegion();
  IdentifierType depthSize = inRegion.GetSize(0);

  output->Allocate(true); // allocate and zero-initialize the image

  IdentifierType inputNumber = 0;
  IdentifierType lineCount = 0;

  // find out how many lines we need to process so we can provide progress
  for (InputDataObjectConstIterator it(this); !it.IsAtEnd(); ++it)
  {
    auto * modifiableInput = const_cast<DataObject *>(it.GetInput());
    input = dynamic_cast<TInputImage *>(modifiableInput);
    if (input) // this input is set
    {
      input->UpdateOutputInformation();
      inRegion = input->GetLargestPossibleRegion();
      if (inRegion.GetSize(0) != depthSize)
      {
        itkExceptionMacro(<< "Input " << inputNumber << " has size " << inRegion.GetSize() << ".\n"
                          << "Size along depth (0) dimension is " << inRegion.GetSize(0) << " but " << depthSize
                          << "was exepcted.");
      }
      if (input->GetSpacing()[0] != output->GetSpacing()[0])
      {
        itkExceptionMacro(<< "Input " << inputNumber << " has spacing " << input->GetSpacing() << ".\n"
                          << "Spacing along depth (0) dimension is " << input->GetSpacing()[0] << " but "
                          << output->GetSpacing()[0] << "was exepcted.");
      }

      lineCount += inRegion.GetNumberOfPixels() / depthSize;
      ++inputNumber;
    }
  }

  TotalProgressReporter progress(this, lineCount);

  // now go through all the inputs and add them to the output
  for (InputDataObjectConstIterator it(this); !it.IsAtEnd(); ++it)
  {
    auto * modifiableInput = const_cast<DataObject *>(it.GetInput());
    input = dynamic_cast<TInputImage *>(modifiableInput);
    if (input) // this input is set
    {
      inRegion = input->GetLargestPossibleRegion();
      input->Update(); // this will trigger reading from file if not already in memory

      typename OutputImageType::IndexType ind1{ 0 }; // initialize all indices zero

      // parallelizing while ensuring correct concurrent writes into the output is tricky
      // as this is probably going to be memory-access limited, just do it single-threaded
      itk::ImageScanlineIterator<InputImageType> iIt(input, inRegion);
      while (!iIt.IsAtEnd())
      {
        while (!iIt.IsAtEndOfLine())
        {
          ind1[0] = iIt.GetIndex()[0]; // set the index along the depth dimension

          OutputPixelType p = output->GetPixel(ind1);
          p += iIt.Get();
          output->SetPixel(ind1, p);

          ++iIt;
        }

        progress.CompletedPixel();
      }
    }
  }

  output->Modified();

  // divide by the number of contributing lines
  using DividerType =
    itk::DivideImageFilter<OutputImageType, itk::Image<IdentifierType, OutputImageDimension>, OutputImageType>;
  typename DividerType::Pointer divider = DividerType::New();
  divider->SetInput1(output);
  divider->SetConstant2(lineCount);
  divider->SetInPlace(true);
  divider->GraftOutput(this->GetOutput());
  divider->Update();
  this->GraftOutput(divider->GetOutput());

  this->UpdateProgress(1.0f);
}


} // end namespace itk

#endif // itkSpectra1DAveragingImageFilter_hxx
