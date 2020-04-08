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
#ifndef itkBlockMatchingImageToImageMetricMetricImageFilter_h
#define itkBlockMatchingImageToImageMetricMetricImageFilter_h

#include "itkBlockMatchingMetricImageFilter.h"
#include "itkBlockMatchingImageToImageMetricMetricImageFilter.h"

namespace itk
{
namespace BlockMatching
{

/** \class ImageToImageMetricMetricImageFilter
 *
 * \brief A BlockMatching::MetricImageFilter that uses any of ITK's
 * ImageToImageMetric's to calculate the metric.
 *
 * \sa MetricImageFilter
 * \sa ImageToImageMetric
 *
 * \ingroup Ultrasound
 */
template< typename TFixedImage, typename TMovingImage, typename TMetricImage >
class ITK_TEMPLATE_EXPORT ImageToImageMetricMetricImageFilter :
  public MetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
{
public:
  /** Standard class type alias. */
  using Self = ImageToImageMetricMetricImageFilter;
  using Superclass = MetricImageFilter< TFixedImage, TMovingImage, TMetricImage >;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ImageToImageMetricMetricImageFilter, MetricImageFilter);

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      Superclass::ImageDimension);

  /** Type of the Moving image. */
  using MovingImageType = typename Superclass::MovingImageType;
  using MovingImageSpacingType = typename MovingImageType::SpacingType;
  using MovingImageRegionType = typename MovingImageType::RegionType;

  /** Type of the Metric image. */
  using MetricImageType = typename Superclass::MetricImageType;
  using MetricImageSpacingType = typename MetricImageType::SpacingType;
  using MetricImageRegionType = typename MetricImageType::RegionType;

  /** Set the spacing in the output metric image.  This is optional.  If not
   * set, the spacing in the moving image will be used. */
  void SetMetricImageSpacing( const MetricImageSpacingType& spacing );
  itkGetConstReferenceMacro( MetricImageSpacing, MetricImageSpacingType );
protected:
  ImageToImageMetricMetricImageFilter();
  virtual ~ImageToImageMetricMetricImageFilter() {}

  virtual void GenerateOutputInformation() override;

  MetricImageSpacingType m_MetricImageSpacing;
  bool                   m_MetricImageSpacingDefined;

private:
  ImageToImageMetricMetricImageFilter( const Self& ); // purposely not implemented
  void operator=( const Self& ); // purposely not implemented

};

} // end namespace BlockMatching
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlockMatchingImageToImageMetricMetricImageFilter.hxx"
#endif

#endif
