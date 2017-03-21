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
#ifndef itkBlockMatchingCosineInterpolationDisplacementCalculator_h
#define itkBlockMatchingCosineInterpolationDisplacementCalculator_h

#include "itkBlockMatchingMetricImageToDisplacementCalculator.h"

namespace itk
{
namespace BlockMatching
{

/** \class CosineInterpolationDisplacementCalculator
 *
 * \brief The displacement around the maximum pixel is interpolated by a cosine
 * in each direction and the peak of the cosine is used.
 *
 * Cespedes et. al. Methods for estimation of subsample time delays of digitized
 * echo signals.  Ultrasonic Imaging 17. 142-171.  1995.
 *
 * \ingroup Ultrasound
 */
template < class TMetricImage, class TDisplacementImage, class TCoordRep=double >
class ITK_TEMPLATE_EXPORT CosineInterpolationDisplacementCalculator:
  public MetricImageToDisplacementCalculator< TMetricImage, TDisplacementImage >
{
public:
  /** Standard class typedefs. */
  typedef CosineInterpolationDisplacementCalculator                               Self;
  typedef MetricImageToDisplacementCalculator< TMetricImage, TDisplacementImage > Superclass;
  typedef SmartPointer< Self >                                                    Pointer;
  typedef SmartPointer< const Self >                                              ConstPointer;

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      Superclass::ImageDimension);

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( CosineInterpolationDisplacementCalculator, MetricImageToDisplacementCalculator );

  typedef typename Superclass::MetricImageType        MetricImageType;
  typedef typename Superclass::MetricImagePointerType MetricImagePointerType;
  typedef typename MetricImageType::PixelType         PixelType;
  typedef typename MetricImageType::SpacingType       SpacingType;
  typedef typename Superclass::PointType              PointType;
  typedef typename Superclass::IndexType              IndexType;

  virtual void SetMetricImagePixel( const PointType & point, const IndexType& index, MetricImageType* image );

  virtual void Compute() {
    // We do this here instead of SetMetricImagePixel so it only has to be done
    // once.
    this->m_DisplacementImage->Modified();
  };

protected:
  CosineInterpolationDisplacementCalculator();

private:
  CosineInterpolationDisplacementCalculator( const Self & );
  void operator=( const Self & );
};

} // end namespace itk
} // end namespace BlockMatching

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlockMatchingCosineInterpolationDisplacementCalculator.hxx"
#endif

#endif
