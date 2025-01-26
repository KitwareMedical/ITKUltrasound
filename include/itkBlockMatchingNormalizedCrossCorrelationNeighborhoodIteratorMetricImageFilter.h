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
#ifndef itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter_h
#define itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter_h

#include "itkBlockMatchingNormalizedCrossCorrelationMetricImageFilter.h"

namespace itk
{
namespace BlockMatching
{

/** \class NormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter
 *
 * \brief Create an image of the the normalized cross correlation with a kernel
 * calculated with a neighborhood iterator.
 *
 * \sa NormalizedCrossCorrelationMetricImageFilter
 *
 * \ingroup Ultrasound
 */
template <typename TFixedImage, typename TMovingImage, typename TMetricImage>
class ITK_TEMPLATE_EXPORT NormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter
  : public NormalizedCrossCorrelationMetricImageFilter<TFixedImage, TMovingImage, TMetricImage>
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(NormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter);

  /** Standard class type alias. */
  using Self = NormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter;
  using Superclass = NormalizedCrossCorrelationMetricImageFilter<TFixedImage, TMovingImage, TMetricImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(NormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter,
               NormalizedCrossCorrelationMetricImageFilter);

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int, TFixedImage::ImageDimension);

  /** Type of the fixed image. */
  using FixedImageType = typename Superclass::FixedImageType;

  /** Type of the moving image. */
  using MovingImageType = typename Superclass::MovingImageType;

  /** Type of the metric image. */
  using MetricImageType = typename Superclass::MetricImageType;
  using MetricImageRegionType = typename MetricImageType::RegionType;
  using MetricImagePixelType = typename MetricImageType::PixelType;

protected:
  NormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter() {}

  void
  BeforeThreadedGenerateData() override;

  void
  DynamicThreadedGenerateData(const MetricImageRegionType & outputRegion) override;

private:
};

} // end namespace BlockMatching
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter.hxx"
#endif

#endif
