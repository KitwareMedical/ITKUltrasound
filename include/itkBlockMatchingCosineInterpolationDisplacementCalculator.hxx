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
#ifndef itkBlockMatchingCosineInterpolationDisplacementCalculator_hxx
#define itkBlockMatchingCosineInterpolationDisplacementCalculator_hxx


#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk
{
namespace BlockMatching
{

template <typename TMetricImage, typename TDisplacementImage, typename TCoordRep>
CosineInterpolationDisplacementCalculator<TMetricImage, TDisplacementImage, TCoordRep>::
  CosineInterpolationDisplacementCalculator()
{}

template <typename TMetricImage, typename TDisplacementImage, typename TCoordRep>
void
CosineInterpolationDisplacementCalculator<TMetricImage, TDisplacementImage, TCoordRep>::SetMetricImagePixel(
  const PointType & centerPoint,
  const IndexType & displacementIndex,
  MetricImageType * metricImage)
{
  Superclass::SetMetricImagePixel(centerPoint, displacementIndex, metricImage);

  // Find index of the maximum value.
  PixelType max = NumericTraits<PixelType>::min();
  IndexType maxIndex;
  maxIndex.Fill(0);

  const typename MetricImageType::RegionType         region = metricImage->GetBufferedRegion();
  ImageRegionConstIteratorWithIndex<MetricImageType> it(metricImage, region);
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    if (it.Get() > max)
    {
      max = it.Get();
      maxIndex = it.GetIndex();
    }
  }

  TCoordRep y1 = max;
  TCoordRep y0;
  TCoordRep y2;
  TCoordRep omega;
  TCoordRep theta;

  PointType maxPoint;
  metricImage->TransformIndexToPhysicalPoint(maxIndex, maxPoint);

  SpacingType spacing = metricImage->GetSpacing();
  IndexType   tempIndex = maxIndex;
  for (unsigned int i = 0; i < ImageDimension; i++)
  {
    tempIndex[i] = maxIndex[i] - 1;
    if (!region.IsInside(tempIndex))
    {
      tempIndex[i] = maxIndex[i];
      continue;
    }
    y0 = metricImage->GetPixel(tempIndex);
    tempIndex[i] = maxIndex[i] + 1;
    if (!region.IsInside(tempIndex))
    {
      tempIndex[i] = maxIndex[i];
      continue;
    }
    y2 = metricImage->GetPixel(tempIndex);
    tempIndex[i] = maxIndex[i];

    omega = std::acos((y0 + y2) / (2 * y1));
    theta = std::atan((y0 - y2) / (2 * y1 * std::sin(omega)));
    // @todo is this rght ?
    maxPoint[i] += spacing[i] / ::itk::Math::pi * -1 * theta / omega;
  }

  this->m_DisplacementImage->SetPixel(displacementIndex, maxPoint - centerPoint);
}

} // end namespace BlockMatching
} // end namespace itk

#endif
