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
#ifndef itkBoxSigmaSqrtNMinusOneImageFilter_hxx
#define itkBoxSigmaSqrtNMinusOneImageFilter_hxx

#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkOffset.h"
#include "itkProgressAccumulator.h"
#include "itkNumericTraits.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkShapedNeighborhoodIterator.h"
#include "itkBoxUtilities.h"

namespace itk
{


template <class TInputImage, class TOutputImage>
BoxSigmaSqrtNMinusOneImageFilter<TInputImage, TOutputImage>::BoxSigmaSqrtNMinusOneImageFilter()
{}


template <class TInputImage, class TOutputImage>
void
BoxSigmaSqrtNMinusOneImageFilter<TInputImage, TOutputImage>
#if ITK_VERSION_MAJOR < 5
  ::ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId)
#else
  ::DynamicThreadedGenerateData(const OutputImageRegionType & outputRegionForThread)
#endif
{

  // Accumulate type is too small
  using AccValueType = typename itk::NumericTraits<PixelType>::RealType;
  using AccPixType = typename itk::Vector<AccValueType, 2>;
  using AccumImageType = typename itk::Image<AccPixType, TInputImage::ImageDimension>;

  typename TInputImage::SizeType internalRadius;
  for (unsigned int i = 0; i < TInputImage::ImageDimension; ++i)
  {
    internalRadius[i] = this->GetRadius()[i] + 1;
  }


  const InputImageType * inputImage = this->GetInput();
  OutputImageType *      outputImage = this->GetOutput();
  RegionType             accumRegion = outputRegionForThread;
  accumRegion.PadByRadius(internalRadius);
  accumRegion.Crop(inputImage->GetRequestedRegion());

#if ITK_VERSION_MAJOR < 5 || defined(ITKV4_COMPATIBILITY)
  // Dummy reporter for compatibility
  ProgressReporter progress(this, 1, 2 * accumRegion.GetNumberOfPixels());
#endif

  typename AccumImageType::Pointer accImage = AccumImageType::New();
  accImage->SetRegions(accumRegion);
  accImage->Allocate();

  BoxSquareAccumulateFunction<TInputImage, AccumImageType>(inputImage,
                                                           accImage.GetPointer(),
                                                           accumRegion,
                                                           accumRegion
#if ITK_VERSION_MAJOR < 5 || defined(ITKV4_COMPATIBILITY)
                                                           ,
                                                           progress);
#else
  );
#endif
  BoxSigmaSqrtNMinusOneCalculatorFunction<AccumImageType, TOutputImage>(
    accImage.GetPointer(), outputImage, accumRegion, outputRegionForThread, this->GetRadius());
}


} // end namespace itk

#endif
