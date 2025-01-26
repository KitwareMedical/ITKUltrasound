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
#ifndef itkBlockMatchingMaximumPixelDisplacementCalculator_h
#define itkBlockMatchingMaximumPixelDisplacementCalculator_h

#include "itkBlockMatchingMetricImageToDisplacementCalculator.h"

namespace itk
{
namespace BlockMatching
{

/** \class MaximumPixelDisplacementCalculator
 *
 * \brief The displacement is taken to be the pixel with the maximum metric
 * value.
 *
 * This is the simplest and fastest of the MetricImageToDisplacementCalculator's.
 *
 * \ingroup Ultrasound
 */
template <typename TMetricImage, typename TDisplacementImage>
class ITK_TEMPLATE_EXPORT MaximumPixelDisplacementCalculator
  : public MetricImageToDisplacementCalculator<TMetricImage, TDisplacementImage>
{
public:
  /** Standard class type alias. */
  using Self = MaximumPixelDisplacementCalculator;
  using Superclass = MetricImageToDisplacementCalculator<TMetricImage, TDisplacementImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(MaximumPixelDisplacementCalculator);

  using MetricImageType = typename Superclass::MetricImageType;
  using MetricImagePointerType = typename Superclass::MetricImagePointerType;
  using PixelType = typename MetricImageType::PixelType;

  using PointType = typename Superclass::PointType;

  using IndexType = typename Superclass::IndexType;

  void
  SetMetricImagePixel(const PointType & point, const IndexType & index, MetricImageType * image) override;

  void
  Compute() override
  {
    // We do this here instead of SetMetricImagePixel so it only has to be done
    // once.
    this->m_DisplacementImage->Modified();
  };

protected:
  MaximumPixelDisplacementCalculator() {}

private:
  MaximumPixelDisplacementCalculator(const Self &);
  void
  operator=(const Self &);
};

} // end namespace BlockMatching
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkBlockMatchingMaximumPixelDisplacementCalculator.hxx"
#endif

#endif
