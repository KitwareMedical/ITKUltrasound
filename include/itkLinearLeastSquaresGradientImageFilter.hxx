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
#ifndef itkLinearLeastSquaresGradientImageFilter_hxx
#define itkLinearLeastSquaresGradientImageFilter_hxx


#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"

#include "itkMath.h"
#include "vnl/algo/vnl_svd.h"

namespace itk
{

template <typename TInputImage, typename TOperatorValueType, typename TOutputValueType>
LinearLeastSquaresGradientImageFilter<TInputImage, TOperatorValueType, TOutputValueType>::
  LinearLeastSquaresGradientImageFilter()
  : m_UseImageSpacing(true)
  , m_UseImageDirection(true)
{
  m_Radius.Fill(1);
}


template <typename TInputImage, typename TOperatorValueType, typename TOutputValueType>
void
LinearLeastSquaresGradientImageFilter<TInputImage, TOperatorValueType, TOutputValueType>::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // get pointers to the input and output
  InputImagePointer  inputPtr = const_cast<InputImageType *>(this->GetInput());
  OutputImagePointer outputPtr = this->GetOutput();

  if (!inputPtr || !outputPtr)
  {
    return;
  }

  // get a copy of the input requested region (should equal the output
  // requested region)
  typename TInputImage::RegionType inputRequestedRegion;
  inputRequestedRegion = inputPtr->GetRequestedRegion();

  // pad the input requested region by the operator radius
  inputRequestedRegion.PadByRadius(m_Radius);

  // crop the input requested region at the input's largest possible region
  if (inputRequestedRegion.Crop(inputPtr->GetLargestPossibleRegion()))
  {
    inputPtr->SetRequestedRegion(inputRequestedRegion);
    return;
  }
  else
  {
    // Couldn't crop the region (requested region is outside the largest
    // possible region).  Throw an exception.

    // store what we tried to request (prior to trying to crop)
    inputPtr->SetRequestedRegion(inputRequestedRegion);

    // build an exception
    InvalidRequestedRegionError e(__FILE__, __LINE__);
    e.SetLocation(ITK_LOCATION);
    e.SetDescription("Requested region is (at least partially) outside the largest possible region.");
    e.SetDataObject(inputPtr);
    throw e;
  }
}


template <typename TInputImage, typename TOperatorValueType, typename TOutputValueType>
TOutputValueType
LinearLeastSquaresGradientImageFilter<TInputImage, TOperatorValueType, TOutputValueType>::GetDerivative(
  const std::slice &                             s,
  const ConstNeighborhoodIterator<TInputImage> & nit)
{
  unsigned int maxConsistentCount = 1;
  unsigned int currentConsistentCount = maxConsistentCount;
  int          previousSgn = itk::Math::sgn0(nit.GetPixel(s.start()));
  int          nextSgn;
  unsigned int maxConsistentIdx = 0;
  unsigned int currentConsistentIdx = maxConsistentIdx;
  unsigned int sliceIdx;
  unsigned int nitIdx;
  for (sliceIdx = 1, nitIdx = s.start() + s.stride(); sliceIdx < s.size(); ++sliceIdx, nitIdx += s.stride())
  {
    nextSgn = itk::Math::sgn0(nit.GetPixel(nitIdx));
    if (nextSgn == previousSgn)
    {
      ++currentConsistentCount;
      if (currentConsistentCount > maxConsistentCount)
      {
        if (currentConsistentCount != maxConsistentCount + 1)
        {
          maxConsistentIdx = currentConsistentIdx;
        }
        maxConsistentCount = currentConsistentCount;
      }
    }
    else
    {
      currentConsistentCount = 1;
      currentConsistentIdx = sliceIdx;
    }
    previousSgn = nextSgn;
  }

  // If the max consistent count is less than the radius, then use the entire
  // slice.
  if (maxConsistentCount < (s.size() - 1) / 2)
  {
    maxConsistentIdx = 0;
    maxConsistentCount = s.size();
  }

  vnl_matrix<OperatorValueType> A(maxConsistentCount, 2);
  A.set_column(1, NumericTraits<OperatorValueType>::One);
  vnl_vector<OperatorValueType> y(maxConsistentCount);
  for (unsigned int i = 0, nitIdx = s.start() + maxConsistentIdx * s.stride(); i < maxConsistentCount;
       ++i, nitIdx += s.stride())
  {
    A[i][0] = i;
    y[i] = nit.GetPixel(nitIdx);
  }

  return static_cast<TOutputValueType>(vnl_svd<OperatorValueType>(A).solve(y)[0]);
}

template <typename TInputImage, typename TOperatorValueType, typename TOutputValueType>
void
LinearLeastSquaresGradientImageFilter<TInputImage, TOperatorValueType, TOutputValueType>::DynamicThreadedGenerateData(
  const OutputImageRegionType & outputRegionForThread)
{
  unsigned int    i;
  OutputPixelType gradient;

  ZeroFluxNeumannBoundaryCondition<InputImageType> nbc;

  ConstNeighborhoodIterator<InputImageType> nit;
  ImageRegionIterator<OutputImageType>      it;

  // Get the input and output
  OutputImageType *      outputImage = this->GetOutput();
  const InputImageType * inputImage = this->GetInput();

  // Find the data-set boundary "faces"
  typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType faceList;
  NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>                        bC;
  faceList = bC(inputImage, outputRegionForThread, m_Radius);

  typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<InputImageType>::FaceListType::iterator fit;
  fit = faceList.begin();

  // Initialize the x_slice array
  nit = ConstNeighborhoodIterator<InputImageType>(m_Radius, inputImage, *fit);

  typename InputImageType::SpacingType spacingScale;
  if (m_UseImageSpacing)
  {
    spacingScale = inputImage->GetSpacing();
  }
  else
  {
    spacingScale.Fill(1);
  }

  // Process non-boundary face and then each of the boundary faces.
  // These are N-d regions which border the edge of the buffer.
  for (fit = faceList.begin(); fit != faceList.end(); ++fit)
  {
    nit = ConstNeighborhoodIterator<InputImageType>(m_Radius, inputImage, *fit);
    it = ImageRegionIterator<OutputImageType>(outputImage, *fit);
    nit.OverrideBoundaryCondition(&nbc);
    nit.GoToBegin();

    while (!nit.IsAtEnd())
    {
      for (i = 0; i < ImageDimension; ++i)
      {
        gradient[i] = this->GetDerivative(nit.GetSlice(i), nit) / spacingScale[i];
      }

      if (this->m_UseImageDirection)
      {
        inputImage->TransformLocalVectorToPhysicalVector(gradient, it.Value());
      }
      else
      {
        it.Value() = gradient;
      }
      ++nit;
      ++it;
    }
  }
}

} // end namespace itk

#endif
