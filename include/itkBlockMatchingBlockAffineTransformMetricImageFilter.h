/*=========================================================================
 *
 *  Copyright Insight Software Consortium
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
template< typename TFixedImage, typename TMovingImage,
          typename TMetricImage, typename TStrainValueType >
class ITK_TEMPLATE_EXPORT BlockAffineTransformMetricImageFilter :
  public MetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
{
public:
  /** Standard class typedefs. */
  typedef BlockAffineTransformMetricImageFilter                        Self;
  typedef MetricImageFilter< TFixedImage, TMovingImage, TMetricImage > Superclass;
  typedef SmartPointer<Self>                                           Pointer;
  typedef SmartPointer<const Self>                                     ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(BlockAffineTransformMetricImageFilter, MetricImageFilter);

  itkNewMacro( Self );

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      Superclass::ImageDimension);

  /** Type of the fixed image. */
  typedef typename Superclass::FixedImageType   FixedImageType;
  typedef typename FixedImageType::ConstPointer FixedImageConstPointerType;

  /** Type of the moving image. */
  typedef typename Superclass::MovingImageType   MovingImageType;
  typedef typename MovingImageType::RegionType   MovingImageRegionType;
  typedef typename MovingImageType::ConstPointer MovingImageConstPointerType;

  /** Type of the metric image. */
  typedef typename Superclass::MetricImageType       MetricImageType;
  typedef typename MetricImageType::Pointer          MetricImagePointerType;
  typedef typename MetricImageType::PixelType        MetricImagePixelType;
  typedef typename Superclass::MetricImageRegionType MetricImageRegionType;

  /** Type of the strain image. */
  typedef Image< SymmetricSecondRankTensor< TStrainValueType, ImageDimension >, ImageDimension > StrainImageType;
  typedef typename StrainImageType::Pointer                                                      StrainImagePointerType;
  typedef LinearInterpolateImageFunction< StrainImageType, double >                              StrainInterpolatorType;

  /** Type of the transform. */
  typedef AffineTransform< MetricImagePixelType, ImageDimension > TransformType;

  /** Type of the interpolator. */
  itkStaticConstMacro(SincWindowRadius, unsigned int, 4);
  typedef Function::LanczosWindowFunction< SincWindowRadius > SincWindowType;
  typedef ZeroFluxNeumannBoundaryCondition< FixedImageType >  SincWindowBoundaryConditionType;
  typedef WindowedSincInterpolateImageFunction< FixedImageType, SincWindowRadius, SincWindowType,
                                                SincWindowBoundaryConditionType, double > InterpolatorType;

  /** Set/Get the Internal MetricImageFilter that actually generates the metric
   * image. */
  itkSetObjectMacro( MetricImageFilter, Superclass );
  itkGetObjectMacro( MetricImageFilter, Superclass );

  /** Set/Get the strain image used to modify the fixed image block. */
  itkSetObjectMacro( StrainImage, StrainImageType );
  itkGetObjectMacro( StrainImage, StrainImageType );

protected:
  BlockAffineTransformMetricImageFilter();

  /** We need the entire input because we don't know where we will be resampling
   * from. */
  virtual void GenerateInputRequestedRegion() ITK_OVERRIDE;

  virtual void GenerateData() ITK_OVERRIDE;

  typename Superclass::Pointer m_MetricImageFilter;
  StrainImagePointerType       m_StrainImage;

  typename TransformType::Pointer          m_Transform;
  typename InterpolatorType::Pointer       m_Interpolator;
  typename StrainInterpolatorType::Pointer m_StrainInterpolator;

  typename FixedImageType::Pointer         m_TransformedFixedImage;

private:
  BlockAffineTransformMetricImageFilter( const Self & );
  void operator=( const Self & );
};

} // end namespace BlockMatching
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlockMatchingBlockAffineTransformMetricImageFilter.hxx"
#endif

#endif
