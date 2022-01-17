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
#ifndef itkSpectra1DNormalizeImageFilter_hxx
#define itkSpectra1DNormalizeImageFilter_hxx


#include "itkImageScanlineIterator.h"
#include "itkTotalProgressReporter.h"

namespace itk
{
template <typename TInputImage, typename TReferenceImage>
void
Spectra1DNormalizeImageFilter<TInputImage, TReferenceImage>::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // ImageToImageFilter expects all inputs to be of TInputImage type
  // so we don't use this->GetInput(), instead we use ProcessObject directly
  auto * in1 = reinterpret_cast<const ReferenceImageType *>(this->ProcessObject::GetInput(1));
  auto * reference = const_cast<ReferenceImageType *>(in1);
  // instead of fiddling with cropping the requested region to reference line,
  // we simply request the entire reference line
  reference->SetRequestedRegionToLargestPossibleRegion();
}

#define DeclareVectorDivideOperator(TypeA, TypeB)                      \
  template <typename TScalarA, typename TScalarB, unsigned VDimension> \
  TypeA ElementWiseDivide(const TypeA & a, const TypeB & b)            \
  {                                                                    \
    TypeA result{ a };                                                 \
    for (unsigned i = 0; i < VDimension; ++i)                          \
    {                                                                  \
      if (b[i] != 0) /* else implicitly assume b[i] == 1.0 */          \
      {                                                                \
        result[i] /= b[i];                                             \
      }                                                                \
    }                                                                  \
    return result;                                                     \
  }                                                                    \
  ITK_NOOP_STATEMENT

// TypeA nor TypeB can have a comma in their specification
// because commas are used to separate macro arguments
#define ITK_COMMA ,

DeclareVectorDivideOperator(Vector<TScalarA ITK_COMMA VDimension>, VariableLengthVector<TScalarB>);
DeclareVectorDivideOperator(VariableLengthVector<TScalarA>, Vector<TScalarB ITK_COMMA VDimension>);
DeclareVectorDivideOperator(Vector<TScalarA ITK_COMMA VDimension>, Vector<TScalarB ITK_COMMA VDimension>);

// VI->VI case needs to be separate, because it does not have compile-time Dimension
template <typename TScalarA, typename TScalarB>
VariableLengthVector<TScalarA>
ElementWiseDivide(const VariableLengthVector<TScalarA> & a, const VariableLengthVector<TScalarB> & b)
{
  VariableLengthVector<TScalarA> result{ a };

  itkAssertInDebugAndIgnoreInReleaseMacro(a.GetNumberOfElements() == b.GetNumberOfElements());
  const unsigned aDimension = a.GetNumberOfElements();
  for (unsigned i = 0; i < aDimension; ++i)
  {
    if (b[i] != 0) /* else implicitly assume b[i] == 1.0 */
    {
      result[i] /= b[i];
    }
  }
  return result;
}


template <typename TInputImage, typename TReferenceImage>
void
Spectra1DNormalizeImageFilter<TInputImage, TReferenceImage>::DynamicThreadedGenerateData(
  const OutputImageRegionType & outputRegionForThread)
{
  const InputImageType * input = this->GetInput();
  const auto *           reference = reinterpret_cast<const ReferenceImageType *>(this->ProcessObject::GetInput(1));
  OutputImageType *      output = this->GetOutput();

  TotalProgressReporter progress(this, output->GetRequestedRegion().GetNumberOfPixels());

  MultiThreaderBase * multiThreader = this->GetMultiThreader();
  multiThreader->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());

  typename ReferenceImageType::IndexType ind1{ 0 }; // initialize all indices to zero

  ImageScanlineConstIterator<InputImageType> iIt(input, outputRegionForThread);
  ImageScanlineIterator<OutputImageType>     oIt(output, outputRegionForThread);

  while (!iIt.IsAtEnd())
  {
    while (!iIt.IsAtEndOfLine())
    {
      ind1[0] = iIt.GetIndex()[0]; // set the index along the depth dimension

      oIt.Set(ElementWiseDivide(iIt.Get(), reference->GetPixel(ind1)));

      ++iIt;
      ++oIt;
    }

    iIt.NextLine();
    oIt.NextLine();
    progress.Completed(outputRegionForThread.GetSize()[0]);
  }
}

} // end namespace itk

#endif // itkSpectra1DNormalizeImageFilter_hxx
