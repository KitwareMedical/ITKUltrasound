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
#ifndef itkBlockMatchingParabolicInterpolationDisplacementCalculator_h
#define itkBlockMatchingParabolicInterpolationDisplacementCalculator_h

#include "itkBlockMatchingMetricImageToDisplacementCalculator.h"

namespace itk
{
namespace BlockMatching
{

/** \class ParabolicInterpolationDisplacementCalculator
 *
 * \brief The displacement around the maximum pixel is interpolated by a
 * parabola in each direction and the peak of the parabola is used.
 *
 * Cespedes et. al. Methods for estimation of subsample time delays of digitized
 * echo signals.  Ultrasonic Imaging 17. 142-171.  1995.
 *
 * \ingroup Ultrasound
 */
template <typename TMetricImage, typename TDisplacementImage, typename TCoordRep = double>
class ITK_TEMPLATE_EXPORT ParabolicInterpolationDisplacementCalculator
  : public MetricImageToDisplacementCalculator<TMetricImage, TDisplacementImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(ParabolicInterpolationDisplacementCalculator);

  /** Standard class type alias. */
  using Self = ParabolicInterpolationDisplacementCalculator;
  using Superclass = MetricImageToDisplacementCalculator<TMetricImage, TDisplacementImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(ParabolicInterpolationDisplacementCalculator);

  using MetricImageType = typename Superclass::MetricImageType;
  using MetricImagePointerType = typename Superclass::MetricImagePointerType;
  using PixelType = typename MetricImageType::PixelType;
  using SpacingType = typename MetricImageType::SpacingType;
  using RegionType = typename MetricImageType::RegionType;

  using MetricImageImageType = typename Superclass::MetricImageImageType;
  using MetricImageImagePointerType = typename Superclass::MetricImageImagePointerType;
  using MetricImageImageIteratorType = ImageRegionIterator<MetricImageImageType>;

  using CenterPointsImageType = typename Superclass::CenterPointsImageType;

  using DisplacementImageType = typename Superclass::DisplacementImageType;

  using PointType = typename Superclass::PointType;
  using IndexType = typename Superclass::IndexType;

  void
  Compute() override;

protected:
  ParabolicInterpolationDisplacementCalculator();

  /** Use a parabolic fit to find the subsample peak. */
  void
  ThreadedParabolicInterpolation(const RegionType & region);

private:
};

} // namespace BlockMatching
} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkBlockMatchingParabolicInterpolationDisplacementCalculator.hxx"
#endif

#endif
