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
#ifndef itkBlockMatchingMaximumPixelDisplacementCalculator_hxx
#define itkBlockMatchingMaximumPixelDisplacementCalculator_hxx

#include "itkBlockMatchingMaximumPixelDisplacementCalculator.h"

#include "itkImageRegionConstIteratorWithIndex.h"

namespace itk
{
namespace BlockMatching
{

template <typename TMetricImage, typename TDisplacementImage>
void
MaximumPixelDisplacementCalculator<TMetricImage, TDisplacementImage>::SetMetricImagePixel(const PointType & point,
                                                                                          const IndexType & index,
                                                                                          MetricImageType * metricImage)
{
  Superclass::SetMetricImagePixel(point, index, metricImage);

  PixelType max = NumericTraits<PixelType>::min();
  IndexType maxIndex;
  maxIndex.Fill(0);

  itk::ImageRegionConstIteratorWithIndex<MetricImageType> it(metricImage, metricImage->GetBufferedRegion());
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    if (it.Get() > max)
    {
      max = it.Get();
      maxIndex = it.GetIndex();
    }
  }

  PointType maxPoint;
  metricImage->TransformIndexToPhysicalPoint(maxIndex, maxPoint);
  this->m_DisplacementImage->SetPixel(index, maxPoint - point);
}

} // end namespace BlockMatching
} // end namespace itk

#endif
