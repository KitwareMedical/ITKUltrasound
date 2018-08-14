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
template < typename TMetricImage, typename TDisplacementImage, typename TCoordRep=double >
class ITK_TEMPLATE_EXPORT ParabolicInterpolationDisplacementCalculator:
  public MetricImageToDisplacementCalculator< TMetricImage, TDisplacementImage >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(ParabolicInterpolationDisplacementCalculator);

  /** Standard class typedefs. */
  typedef ParabolicInterpolationDisplacementCalculator                            Self;
  typedef MetricImageToDisplacementCalculator< TMetricImage, TDisplacementImage > Superclass;
  typedef SmartPointer< Self >                                                    Pointer;
  typedef SmartPointer< const Self >                                              ConstPointer;

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int, Superclass::ImageDimension);

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ParabolicInterpolationDisplacementCalculator, MetricImageToDisplacementCalculator );

  typedef typename Superclass::MetricImageType        MetricImageType;
  typedef typename Superclass::MetricImagePointerType MetricImagePointerType;
  typedef typename MetricImageType::PixelType         PixelType;
  typedef typename MetricImageType::SpacingType       SpacingType;
  typedef typename MetricImageType::RegionType        RegionType;

  typedef typename Superclass::MetricImageImageType MetricImageImageType;
  typedef typename Superclass::MetricImageImagePointerType
    MetricImageImagePointerType;
  typedef ImageRegionIterator< MetricImageImageType >
    MetricImageImageIteratorType;

  typedef typename Superclass::CenterPointsImageType CenterPointsImageType;

  typedef typename Superclass::DisplacementImageType DisplacementImageType;

  typedef typename Superclass::PointType PointType;
  typedef typename Superclass::IndexType IndexType;

  void Compute() override;

protected:
  ParabolicInterpolationDisplacementCalculator();

  /** Use a parabolic fit to find the subsample peak. */
  void ThreadedParabolicInterpolation( const RegionType& region );

private:
};

} // end namespace itk
} // end namespace BlockMatching

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlockMatchingParabolicInterpolationDisplacementCalculator.hxx"
#endif

#endif
