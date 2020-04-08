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
#ifndef itkBlockMatchingNormalizedCrossCorrelationMetricImageFilter_h
#define itkBlockMatchingNormalizedCrossCorrelationMetricImageFilter_h

#include "itkBoxMeanImageFilter.h"
#include "itkBoxSigmaSqrtNMinusOneImageFilter.h"

#include "itkBlockMatchingMetricImageFilter.h"

namespace itk
{
namespace BlockMatching
{

/** \class NormalizedCrossCorrelationMetricImageFilter
 *
 * \brief A MetricImageFilter that calculates a metric image of
 * normalized cross correlation.
 *
 * r_{xy} = \frac{1}{n-1} \sum_{x,y}\frac{(x) - \overline{f})(y -
 * \overline{y})}{\sigma_f \sigma_y}
 *
 * This is an abstract base class that does the mean and standard deviation
 * calculation.  The cross correlation is left to inherited classes.
 *
 * \sa MetricImageFilter
 *
 * \ingroup Ultrasound
 */
template< class TFixedImage, class TMovingImage, class TMetricImage >
class ITK_TEMPLATE_EXPORT NormalizedCrossCorrelationMetricImageFilter :
  public MetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(NormalizedCrossCorrelationMetricImageFilter);

  /** Standard class type alias. */
  using Self = NormalizedCrossCorrelationMetricImageFilter;
  using Superclass = MetricImageFilter< TFixedImage, TMovingImage, TMetricImage >;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Run-time type information (and related methods). */
  itkTypeMacro(NormalizedCrossCorrelationMetricImageFilter, MetricImageFilter);

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

protected:
  NormalizedCrossCorrelationMetricImageFilter();

  /** The mean and pseudo-standarddeviation images are stored in the outputs so
    they fix in with the pipline architecture. */
  virtual void GenerateOutputInformation() override;

  /** All outputs generate the largest possible region. */
  virtual void EnlargeOutputRequestedRegion( DataObject * data ) override;

  /** Don't let the default mess with our output requested regions. */
  virtual void GenerateOutputRequestedRegion( DataObject * data ) override {};

  /** Generates helper images for the calculation.  These are only needed for
   * internal calculation, but they are put on the
   * outputs for use by subclasses and to use the pipeline memory management
   * system.
   *
   * Calculates an image of means for each block neighborhood.  Also calculates
   * an image of standard deviations (times sqrt(N-1)) for each block neighborhood. */
  virtual void GenerateHelperImages();

  using BoxMeanFilterType = BoxMeanImageFilter< MovingImageType, MetricImageType >;
  using BoxPseudoSigmaFilterType = BoxSigmaSqrtNMinusOneImageFilter< MovingImageType, MetricImageType >;

  typename BoxMeanFilterType::Pointer        m_BoxMeanFilter;
  typename BoxPseudoSigmaFilterType::Pointer m_BoxPseudoSigmaFilter;

private:
  using BoundaryConditionType = ConstantBoundaryCondition< MetricImageType >;
  BoundaryConditionType m_BoundaryCondition;
};

} // end namespace BlockMatching
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlockMatchingNormalizedCrossCorrelationMetricImageFilter.hxx"
#endif

#endif
