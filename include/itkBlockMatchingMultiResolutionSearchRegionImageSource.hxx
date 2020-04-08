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
#ifndef itkBlockMatchingMultiResolutionSearchRegionImageSource_hxx
#define itkBlockMatchingMultiResolutionSearchRegionImageSource_hxx

#include "itkBlockMatchingMultiResolutionSearchRegionImageSource.h"

namespace itk
{
namespace BlockMatching
{

template <typename TFixedImage, typename TMovingImage, typename TDisplacementImage>
MultiResolutionSearchRegionImageSource<TFixedImage, TMovingImage, TDisplacementImage>::
  MultiResolutionSearchRegionImageSource()
  : m_CurrentLevel(0)
{
  m_PreviousDisplacements = DisplacementImageType::New();
  m_DisplacementDuplicator = DisplacementDuplicatorType::New();
  m_DisplacementResampler = DisplacementResamplerType::New();
}


template <typename TFixedImage, typename TMovingImage, typename TDisplacementImage>
void
MultiResolutionSearchRegionImageSource<TFixedImage, TMovingImage, TDisplacementImage>::SetPreviousDisplacements(
  const DisplacementImageType * displacements)
{
  // We make a copy because the ImageRegistrationMethods will mess with the
  // output information on GenerateOutputInformation before the
  // SearchRegionImageSource has a chance to use it in GenerateData.
  m_DisplacementDuplicator->SetInputImage(displacements);
  m_DisplacementDuplicator->Update();
  m_PreviousDisplacements = m_DisplacementDuplicator->GetOutput();
}


template <typename TFixedImage, typename TMovingImage, typename TDisplacementImage>
void
MultiResolutionSearchRegionImageSource<TFixedImage, TMovingImage, TDisplacementImage>::SetOverlapSchedule(
  const double & schedule)
{
  // Check to make sure the PyramidSchedule has been set.
  if (m_PyramidSchedule.size() == 0)
  {
    itkExceptionMacro(<< "The PyramidSchedule must be set before calling this method.");
  }

  m_OverlapSchedule.SetSize(m_PyramidSchedule.rows(), m_PyramidSchedule.cols());
  m_OverlapSchedule.Fill(schedule);
}


template <typename TFixedImage, typename TMovingImage, typename TDisplacement>
void
MultiResolutionSearchRegionImageSource<TFixedImage, TMovingImage, TDisplacement>::GenerateOutputInformation()
{
  OutputImageType * output = this->GetOutput();

  if (m_FixedImage.IsNull())
  {
    itkExceptionMacro(<< "Fixed Image is not present.");
  }
  m_FixedImage->UpdateOutputInformation();
  output->SetDirection(m_FixedImage->GetDirection());

  // Set origin.  The first block is in the corner of the fixed image.
  typename FixedImageType::PointType  fixedOrigin = m_FixedImage->GetOrigin();
  typename OutputImageType::PointType origin;

  typename FixedImageType::SpacingType fixedSpacing = m_FixedImage->GetSpacing();

  RadiusType nullRadius;
  nullRadius.Fill(0);
  if (m_FixedBlockRadius == nullRadius)
  {
    itkExceptionMacro(<< "The FixedBlockRadius has not been set.");
  }

  if (m_OverlapSchedule.size() == 0)
  {
    itkExceptionMacro(<< "OverlapSchedule is not present.");
  }

  typename FixedImageType::IndexType fixedIndex = m_FixedImage->GetLargestPossibleRegion().GetIndex();
  for (unsigned int i = 0; i < ImageDimension; ++i)
  {
    origin[i] = fixedOrigin[i] + +(fixedIndex[i] + m_FixedBlockRadius[i]) * fixedSpacing[i];
  }
  output->SetOrigin(origin);

  typename OutputImageType::SpacingType spacing;
  for (unsigned int i = 0; i < ImageDimension; ++i)
  {
    spacing[i] = 2 * m_FixedBlockRadius[i] * fixedSpacing[i] * m_OverlapSchedule(m_CurrentLevel, i);
  }
  output->SetSpacing(spacing);

  OutputRegionType                     region;
  typename OutputRegionType::IndexType index;
  index.Fill(0);
  region.SetIndex(index);
  typename OutputRegionType::SizeType size;
  typename FixedImageType::SizeType   fixedSize = m_FixedImage->GetLargestPossibleRegion().GetSize();
  for (unsigned int i = 0; i < ImageDimension; ++i)
  {
    size[i] = static_cast<SizeValueType>(
      std::floor((fixedSize[i] - m_FixedBlockRadius[i] * 2 - 2) * fixedSpacing[i] / spacing[i]) - 2);
  }
  region.SetSize(size);
  output->SetLargestPossibleRegion(region);
}


template <typename TFixedImage, typename TMovingImage, typename TDisplacement>
void
MultiResolutionSearchRegionImageSource<TFixedImage, TMovingImage, TDisplacement>::BeforeThreadedGenerateData()
{
  if (this->m_CurrentLevel != 0)
  {
    // ! @todo these resampler should be replaced by resamplers for each
    // component that can specify a neumann boundary condition
    // ditto with FixedSearchRegionImageSource
    m_DisplacementResampler->SetInput(this->m_PreviousDisplacements);
    typename OutputImageType::Pointer output = this->GetOutput();
    m_DisplacementResampler->SetSize(output->GetRequestedRegion().GetSize());
    m_DisplacementResampler->SetOutputStartIndex(output->GetRequestedRegion().GetIndex());
    m_DisplacementResampler->SetOutputSpacing(output->GetSpacing());
    m_DisplacementResampler->SetOutputOrigin(output->GetOrigin());
    m_DisplacementResampler->SetOutputDirection(output->GetDirection());
    m_DisplacementResampler->UpdateLargestPossibleRegion();
  }
}

} // namespace BlockMatching
} // namespace itk

#endif
