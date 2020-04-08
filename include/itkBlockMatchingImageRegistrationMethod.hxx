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
#ifndef itkBlockMatchingImageRegistrationMethod_hxx
#define itkBlockMatchingImageRegistrationMethod_hxx

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkProgressReporter.h"

#include "itkBlockMatchingMaximumPixelDisplacementCalculator.h"
#include "itkBlockMatchingImageRegistrationMethod.h"


namespace itk
{
namespace BlockMatching
{

template <typename TFixedImage,
          typename TMovingImage,
          typename TMetricImage,
          typename TDisplacementImage,
          typename TCoordRep>
ImageRegistrationMethod<TFixedImage, TMovingImage, TMetricImage, TDisplacementImage, TCoordRep>::
  ImageRegistrationMethod()
  : m_UseStreaming(false)
{
  m_FixedImage = nullptr;
  m_MovingImage = nullptr;
  m_MetricImageFilter = nullptr;
  m_MetricImageToDisplacementCalculator = MaximumPixelDisplacementCalculator<TMetricImage, TDisplacementImage>::New();

  m_Radius.Fill(0);
}


template <typename TFixedImage,
          typename TMovingImage,
          typename TMetricImage,
          typename TDisplacementImage,
          typename TCoordRep>
void
ImageRegistrationMethod<TFixedImage, TMovingImage, TMetricImage, TDisplacementImage, TCoordRep>::SetFixedImage(
  FixedImageType * fixedImage)
{
  if (this->m_FixedImage.GetPointer() != fixedImage)
  {
    this->m_FixedImage = fixedImage;
    this->Modified();
  }
}


template <typename TFixedImage,
          typename TMovingImage,
          typename TMetricImage,
          typename TDisplacementImage,
          typename TCoordRep>
void
ImageRegistrationMethod<TFixedImage, TMovingImage, TMetricImage, TDisplacementImage, TCoordRep>::SetMovingImage(
  MovingImageType * movingImage)
{
  if (this->m_MovingImage.GetPointer() != movingImage)
  {
    this->m_MovingImage = movingImage;
    this->Modified();
  }
}


template <typename TFixedImage,
          typename TMovingImage,
          typename TMetricImage,
          typename TDisplacementImage,
          typename TCoordRep>
void
ImageRegistrationMethod<TFixedImage, TMovingImage, TMetricImage, TDisplacementImage, TCoordRep>::
  GenerateOutputInformation()
{
  Superclass::GenerateOutputInformation();

  typename ImageType::Pointer output = this->GetOutput();

  // We do this here instead of GenerateData() so the
  // MetricImageToDisplacementCalculator has a chance to modify the input
  // requested region based on the displacement requested region.
  this->m_MetricImageToDisplacementCalculator->SetDisplacementImage(output.GetPointer());
}


template <typename TFixedImage,
          typename TMovingImage,
          typename TMetricImage,
          typename TDisplacementImage,
          typename TCoordRep>
void
ImageRegistrationMethod<TFixedImage, TMovingImage, TMetricImage, TDisplacementImage, TCoordRep>::
  GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();

  typename SearchRegionImageType::Pointer input = const_cast<SearchRegionImageType *>(this->GetInput());
  if (!input)
  {
    itkExceptionMacro(<< "Input SearchRegionImage is not present.");
  }

  MovingRegionType region = input->GetRequestedRegion();
  this->m_MetricImageToDisplacementCalculator->ModifyGenerateInputRequestedRegion(region);
  input->SetRequestedRegion(region);
}


template <typename TFixedImage,
          typename TMovingImage,
          typename TMetricImage,
          typename TDisplacementImage,
          typename TCoordRep>
void
ImageRegistrationMethod<TFixedImage, TMovingImage, TMetricImage, TDisplacementImage, TCoordRep>::
  EnlargeOutputRequestedRegion(DataObject * data)
{
  this->m_MetricImageToDisplacementCalculator->ModifyEnlargeOutputRequestedRegion(data);
}


template <typename TFixedImage,
          typename TMovingImage,
          typename TMetricImage,
          typename TDisplacementImage,
          typename TCoordRep>
void
ImageRegistrationMethod<TFixedImage, TMovingImage, TMetricImage, TDisplacementImage, TCoordRep>::GenerateData()
{
  const SearchRegionImageType * input = this->GetInput();

  this->Initialize();

  ImageType * output = this->GetOutput();

  RegionType requestedRegion = output->GetRequestedRegion();
  using IteratorType = ImageRegionIteratorWithIndex<ImageType>;
  IteratorType it(output, requestedRegion);
  using SearchRegionImageIteratorType = ImageRegionConstIterator<SearchRegionImageType>;
  SearchRegionImageIteratorType searchIt(input, requestedRegion);

  // The fixed image region is the kernel block size, and its size is constant.
  FixedRegionType                     fixedRegion;
  typename FixedRegionType::SizeType  fixedSize;
  typename FixedRegionType::IndexType fixedIndex;
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    fixedSize[i] = m_Radius[i] * 2 + 1;
  }
  fixedRegion.SetSize(fixedSize);
  FixedRegionType fixedLargestPossibleRegion = m_FixedImage->GetLargestPossibleRegion();

  CoordRepType coord;

  if (!m_UseStreaming)
  {
    m_FixedImage->Update();
    m_MovingImage->Update();
    m_FixedImage->DisconnectPipeline();
    m_MovingImage->DisconnectPipeline();
  }

  // Note that this may not be accurate if
  // m_MetricImageToDisplacementCalculator->Compute() takes a long time.  In
  // that case one may want to monitor the progress of
  // m_MetricImageToDisplacementCalculator separately.
  ProgressReporter progress(this, 0, requestedRegion.GetNumberOfPixels());

  for (it.GoToBegin(), searchIt.GoToBegin(); !it.IsAtEnd(); ++it, ++searchIt)
  {
    output->TransformIndexToPhysicalPoint(it.GetIndex(), coord);
    m_FixedImage->TransformPhysicalPointToIndex(coord, fixedIndex);
    for (unsigned int i = 0; i < ImageDimension; ++i)
    {
      fixedIndex[i] -= m_Radius[i];
    }
    fixedRegion.SetIndex(fixedIndex);
    m_MetricImageFilter->SetFixedImageRegion(fixedRegion);
    m_MetricImageFilter->SetMovingImageRegion(searchIt.Get());
    m_MetricImageFilter->Update();
    m_MetricImageToDisplacementCalculator->SetMetricImagePixel(coord, it.GetIndex(), m_MetricImageFilter->GetOutput());
    progress.CompletedPixel();
  }

  m_MetricImageToDisplacementCalculator->Compute();
}


template <typename TFixedImage,
          typename TMovingImage,
          typename TMetricImage,
          typename TDisplacementImage,
          typename TCoordRep>
void
ImageRegistrationMethod<TFixedImage, TMovingImage, TMetricImage, TDisplacementImage, TCoordRep>::Initialize()
{
  if (!m_FixedImage)
  {
    itkExceptionMacro(<< "FixedImage is not present.");
  }
  m_FixedImage->UpdateOutputInformation();

  if (!m_MovingImage)
  {
    itkExceptionMacro(<< "MovingImage is not present.");
  }
  m_MovingImage->UpdateOutputInformation();

  if (!m_MetricImageFilter)
  {
    itkExceptionMacro(<< "BlockMatching::MetricImageFilter is not present.");
  }

  RadiusType nullRadius;
  nullRadius.Fill(0);
  if (nullRadius == m_Radius)
  {
    itkExceptionMacro(<< "The block radius has not been set.");
  }

  m_MetricImageFilter->SetFixedImage(m_FixedImage);
  m_MetricImageFilter->SetMovingImage(m_MovingImage);

  this->AllocateOutputs();
  typename ImageType::Pointer output = this->GetOutput();
  if (!output)
    return;

  typename SearchRegionImageType::ConstPointer input = this->GetInput();
  if (!input)
  {
    itkExceptionMacro(<< "Input SearchRegionImage is not present.");
  }
}

} // end namespace BlockMatching
} // end namespace itk

#endif
