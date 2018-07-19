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
#ifndef itkBlockMatchingStrainWindowDisplacementCalculator_h
#define itkBlockMatchingStrainWindowDisplacementCalculator_h

#include "itkBlockMatchingMetricImageToDisplacementCalculator.h"

#include "itkAbsImageFilter.h"
#include "itkBoxMeanImageFilter.h"
#include "itkImageToImageFilter.h"
#include "itkNthElementImageAdaptor.h"
#include "itkSymmetricSecondRankTensor.h"

namespace itk
{
namespace BlockMatching
{

/** \class StrainWindowDisplacementCalculator
 *
 * \brief Displacements outside a strain window are interpolated by surrounding
 * displacements.
 *
 * This class acts a filter to remove 'peak-hopping' artifacts in the
 * displacement.
 *
 * This class does not directly calculate displacement, but relies on another
 * MetricImageToDisplacementCalculator to first calculate the displacements.
 * Then, the resulting displacements are filtered.  Points whose strain values
 * lie outside the MaximumAbsStrain are
 * interpolated by surrounding displacements or extrapolated using the nearest
 * strain value if possible.
 *
 * The filter is templated over the MetricImage (for superclass, not used),
 * displacement image type, and strain image component value type.
 *
 * The filter can be applied iteratively to remove regions of high strain with
 * SetMaximumIterations().
 *
 * \ingroup Ultrasound
 */
template <typename TMetricImage, typename TDisplacementImage, typename TStrainValueType>
class ITK_TEMPLATE_EXPORT StrainWindowDisplacementCalculator :
  public MetricImageToDisplacementCalculator< TMetricImage, TDisplacementImage >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(StrainWindowDisplacementCalculator);

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int, TDisplacementImage::ImageDimension);

  /** Standard class typedefs. */
  typedef StrainWindowDisplacementCalculator                                    Self;
  typedef MetricImageToDisplacementCalculator<TMetricImage, TDisplacementImage> Superclass;
  typedef SmartPointer<Self>                                                    Pointer;
  typedef SmartPointer<const Self>                                              ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( StrainWindowDisplacementCalculator, MetricImageToDisplacementCalculator );

  typedef typename Superclass::DisplacementImageType DisplacementImageType;
  typedef typename DisplacementImageType::RegionType RegionType;

  typedef ImageToImageFilter<DisplacementImageType,
                                  Image<SymmetricSecondRankTensor<TStrainValueType, ImageDimension>, ImageDimension> >
  StrainImageFilterType;
  typedef typename StrainImageFilterType::OutputImagePixelType StrainTensorType;
  typedef typename StrainImageFilterType::OutputImageType      StrainImageType;

  /** @todo: Should this filter define ModifyGenerateInputRequestedRegion and ModifyEnlargeOutputRequestedRegion
   * so the LargestPossibleRegion is defined?  Behavior might not be as good on
   * the edges of the image, but then streaming is not possible.  Maybe it
   * should be optional? */

  typedef typename Superclass::PointType       PointType;
  typedef typename Superclass::IndexType       IndexType;
  typedef typename Superclass::MetricImageType MetricImageType;

  virtual void SetDisplacementImage( DisplacementImageType* image )
    {
    Superclass::SetDisplacementImage( image );
    m_DisplacementCalculator->SetDisplacementImage( image );
    }

  virtual void SetMetricImagePixel( const PointType & point,
    const IndexType& index,
    MetricImageType * image );

  void Compute() override;

  /** Set/Get the internal displacement calculator that is used to calculate the
   * displacements after regularization.  Defaults to a
   * MaximumPixelDisplacementCalcultor. */
  itkSetObjectMacro( DisplacementCalculator, Superclass );
  itkGetObjectMacro( DisplacementCalculator, Superclass );

  /** Set/Get the strain image filter that is used to calculate the strain
   * image. This defaults to an itk::StrainImageFilter. */
  itkSetObjectMacro( StrainImageFilter, StrainImageFilterType );
  itkGetObjectMacro( StrainImageFilter, StrainImageFilterType );

  /** Set/Get the maxmimum absolute strain allowed.  If any of strain components are below
   * the values in this tensor, the displacement is interpolated or extrapolated
   * by surrounding values. */
  virtual void SetMaximumAbsStrain( const StrainTensorType& max )
    {
    if( m_MaximumAbsStrain != max )
      {
      m_MaximumAbsStrain = max;
      this->Modified();
      }
    }
  itkGetConstReferenceMacro( MaximumAbsStrain, StrainTensorType );

  /** Set/Get the maximum iterations the algorithm is applied.  The algorithm
   * will be applied repeatedly until the number of strain pixels outside the
   * specified window converges or this value is reached.  Defaults to 1. */
  itkSetMacro( MaximumIterations, unsigned int );
  itkGetConstMacro( MaximumIterations, unsigned int );

  /** Get the current or last iteration of the algorithm.  Starts from 1.  See
   * also SetMaximumIterations(). */
  itkGetConstMacro( CurrentIteration, unsigned int );

protected:
  StrainWindowDisplacementCalculator();

  /** Determine the values outside the strain window, and set the internal mask
   * image to true if the value needs to be replaced.  Return the number of
   * values outside the strain window. */
  unsigned long long GenerateMask( const RegionType & region );

  /** Replace the displacements where the mask is true. */
  void ReplaceDisplacements( const RegionType & region );

  typename Superclass::Pointer m_DisplacementCalculator;

  typename StrainImageFilterType::Pointer m_StrainImageFilter;

  StrainTensorType m_MinimumStrain;
  StrainTensorType m_MaximumAbsStrain;

  unsigned int m_CurrentIteration;
  unsigned int m_MaximumIterations;

  typedef Image<bool, ImageDimension> MaskType;
  typename MaskType::Pointer m_Mask;

  typedef typename MetricImageType::PixelType                        MetricPixelType;
  typedef NthElementImageAdaptor< StrainImageType, MetricPixelType > NthElementAdaptorType;
  typename NthElementAdaptorType::Pointer                            m_NthElementAdaptor;

  typedef AbsImageFilter< NthElementAdaptorType, MetricImageType > AbsFilterType;
  typename AbsFilterType::Pointer                                  m_AbsFilter;

  // There are cases where you will get a large negative pixel, followed by a
  // normal pixel, followed by a large positive strain pixel, for instance.
  // This filter helps address that problem.
  typedef BoxMeanImageFilter< MetricImageType, MetricImageType > BoxMeanImageFilterType;
  typename BoxMeanImageFilterType::Pointer                       m_BoxMeanFilter;

private:
};

} // end namespace itk
} // end namespace BlockMatching

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlockMatchingStrainWindowDisplacementCalculator.hxx"
#endif

#endif
