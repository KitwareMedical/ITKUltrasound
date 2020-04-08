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
#ifndef itkBlockMatchingBlockAffineTransformMetricImageFilter_h
#define itkBlockMatchingBlockAffineTransformMetricImageFilter_h

#include "itkBlockMatchingMetricImageFilter.h"

#include "itkAffineTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkWindowedSincInterpolateImageFunction.h"
#include "itkZeroFluxNeumannBoundaryCondition.h"

namespace itk
{
namespace BlockMatching
{

/** \class BlockAffineTransformMetricImageFilter
 *
 * \brief Modifies the fixed image block according to the affine transform
 * implied by a strain image before passing to a delegate MetricImageFilter.
 *
 * For now we just apply the scaling implied by the strain image.  \todo : apply
 * the shearing implied by the strain image too?
 *
 * \sa MetricImageFilter
 *
 * \ingroup Ultrasound
 */
template <typename TFixedImage, typename TMovingImage, typename TMetricImage, typename TStrainValueType>
class ITK_TEMPLATE_EXPORT BlockAffineTransformMetricImageFilter
  : public MetricImageFilter<TFixedImage, TMovingImage, TMetricImage>
{
public:
  /** Standard class type alias. */
  using Self = BlockAffineTransformMetricImageFilter;
  using Superclass = MetricImageFilter<TFixedImage, TMovingImage, TMetricImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information (and related methods). */
  itkTypeMacro(BlockAffineTransformMetricImageFilter, MetricImageFilter);

  itkNewMacro(Self);

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Type of the fixed image. */
  using FixedImageType = typename Superclass::FixedImageType;
  using FixedImageConstPointerType = typename FixedImageType::ConstPointer;

  /** Type of the moving image. */
  using MovingImageType = typename Superclass::MovingImageType;
  using MovingImageRegionType = typename MovingImageType::RegionType;
  using MovingImageConstPointerType = typename MovingImageType::ConstPointer;

  /** Type of the metric image. */
  using MetricImageType = typename Superclass::MetricImageType;
  using MetricImagePointerType = typename MetricImageType::Pointer;
  using MetricImagePixelType = typename MetricImageType::PixelType;
  using MetricImageRegionType = typename Superclass::MetricImageRegionType;

  /** Type of the strain image. */
  using StrainImageType = Image<SymmetricSecondRankTensor<TStrainValueType, ImageDimension>, ImageDimension>;
  using StrainImagePointerType = typename StrainImageType::Pointer;
  using StrainInterpolatorType = LinearInterpolateImageFunction<StrainImageType, double>;

  /** Type of the transform. */
  using TransformType = AffineTransform<MetricImagePixelType, ImageDimension>;

  /** Type of the interpolator. */
  itkStaticConstMacro(SincWindowRadius, unsigned int, 4);
  using SincWindowType = Function::LanczosWindowFunction<SincWindowRadius>;
  using SincWindowBoundaryConditionType = ZeroFluxNeumannBoundaryCondition<FixedImageType>;
  using InterpolatorType = WindowedSincInterpolateImageFunction<FixedImageType,
                                                                SincWindowRadius,
                                                                SincWindowType,
                                                                SincWindowBoundaryConditionType,
                                                                double>;

  /** Set/Get the Internal MetricImageFilter that actually generates the metric
   * image. */
  itkSetObjectMacro(MetricImageFilter, Superclass);
  itkGetConstObjectMacro(MetricImageFilter, Superclass);

  /** Set/Get the strain image used to modify the fixed image block. */
  itkSetObjectMacro(StrainImage, StrainImageType);
  itkGetConstObjectMacro(StrainImage, StrainImageType);

protected:
  BlockAffineTransformMetricImageFilter();

  /** We need the entire input because we don't know where we will be resampling
   * from. */
  virtual void
  GenerateInputRequestedRegion() override;

  virtual void
  GenerateData() override;

  typename Superclass::Pointer m_MetricImageFilter;
  StrainImagePointerType       m_StrainImage;

  typename TransformType::Pointer          m_Transform;
  typename InterpolatorType::Pointer       m_Interpolator;
  typename StrainInterpolatorType::Pointer m_StrainInterpolator;

  typename FixedImageType::Pointer m_TransformedFixedImage;

private:
  BlockAffineTransformMetricImageFilter(const Self &);
  void
  operator=(const Self &);
};

} // end namespace BlockMatching
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkBlockMatchingBlockAffineTransformMetricImageFilter.hxx"
#endif

#endif
