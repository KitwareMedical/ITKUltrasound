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
#ifndef itkSpeckleReducingAnisotropicDiffusionFunction_h
#define itkSpeckleReducingAnisotropicDiffusionFunction_h

#include <cmath>

#include "itkScalarAnisotropicDiffusionFunction.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkDerivativeOperator.h"
#include "itkLaplacianOperator.h"
#include "itkMacro.h"

namespace itk
{
/** \class SpeckleReducingAnisotropicDiffusionFunction
 *
 * This class implements a 2-dimensional version of Yu-Acton
 * speckle reducing anisotropic diffusion (SRAD) for scalar-valued images.  See
 * itkAnisotropicDiffusionFunction for an overview of the anisotropic diffusion
 * framework and equation.
 *
 * \par References
 * Yongjian Yu and S. T. Acton, "Speckle reducing anisotropic diffusion,"
 * in IEEE Transactions on Image Processing, vol. 11, no. 11, pp. 1260-1270,
 * Nov. 2002, doi: 10.1109/TIP.2002.804276.
 *
 * \sa AnisotropicDiffusionFunction
 * \sa VectorAnisotropicDiffusionFunction
 * \sa VectorGradientAnisotropicDiffusionFunction
 * \sa CurvatureNDAnisotropicDiffusionFunction
 * \ingroup FiniteDifferenceFunctions
 * \ingroup ImageEnhancement
 * \ingroup ITKAnisotropicSmoothing
 * \ingroup ITKUltrasound
 */
template <typename TImage>
class ITK_TEMPLATE_EXPORT SpeckleReducingAnisotropicDiffusionFunction
  : public ScalarAnisotropicDiffusionFunction<TImage>
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(SpeckleReducingAnisotropicDiffusionFunction);

  /** Standard class type aliases. */
  using Self = SpeckleReducingAnisotropicDiffusionFunction;
  using Superclass = ScalarAnisotropicDiffusionFunction<TImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods) */
  itkOverrideGetNameOfClassMacro(SpeckleReducingAnisotropicDiffusionFunction);

  /** Inherit some parameters from the superclass type. */
  using ImageType = typename Superclass::ImageType;
  using PixelType = typename Superclass::PixelType;
  using PixelRealType = typename Superclass::PixelRealType;
  using TimeStepType = typename Superclass::TimeStepType;
  using RadiusType = typename Superclass::RadiusType;
  using NeighborhoodType = typename Superclass::NeighborhoodType;
  using FloatOffsetType = typename Superclass::FloatOffsetType;

  using SizeType = typename ImageType::SizeType;
  using RegionType = typename ImageType::RegionType;
  using MetricImageType = ImageType;
  using MetricImagePointerType = typename ImageType::Pointer;
  using NeighborhoodSizeValueType = SizeValueType;

  /** Inherit some parameters from the superclass type. */
  static constexpr unsigned int ImageDimension = Superclass::ImageDimension;

  /** Compute the equation value. */
  PixelType
  ComputeUpdate(const NeighborhoodType & it,
                void *                   globalData,
                const FloatOffsetType &  offset = FloatOffsetType(0.0)) override;

  /** This method is called prior to each iteration of the solver. */
  void
  InitializeIteration() override;

  PixelType
  ComputeDiffusionCoefficient(PixelType intensity, PixelType averageGradient, PixelType laplacian);

  /** Set/Get the rho parameter used in calculating the
   *  speckle scale function q_0(t). */
  void
  SetRho(const PixelType & rho)
  {
    m_rho = rho;
  }
  const PixelType &
  GetRho() const
  {
    return m_rho;
  }

  /** Set/Get the speckle coefficient of variation parameter q0
   * used in calculating the speckle scale function q_0(t). */
  void
  SetQ0(const PixelType & q0)
  {
    m_q0 = q0;
  }
  const PixelType &
  GetQ0() const
  {
    return m_q0;
  }

  void
  SetElapsedIterationCount(const IdentifierType & iteration)
  {
    m_ElapsedIterationCount = iteration;
  }
  const IdentifierType &
  GetElapsedIterationCount()
  {
    return m_ElapsedIterationCount;
  }

protected:
  SpeckleReducingAnisotropicDiffusionFunction();
  ~SpeckleReducingAnisotropicDiffusionFunction() override = default;

  /** Inner product function. */
  NeighborhoodInnerProduct<ImageType> m_InnerProduct;

  /** Slices for computing neighborhood gradient and Laplacian.
   *  Abstracts pixel stride to allow operator to step in an
   *  arbitrary pixel direction. */
  std::slice x_slice[ImageDimension];

  /** Derivative operator. */
  DerivativeOperator<PixelType, Self::ImageDimension> m_DerivativeOperator;

  /** Laplacian operator. */
  LaplacianOperator<PixelType, Self::ImageDimension> m_LaplacianOperator;
  // DerivativeOperator<PixelType, Self::ImageDimension> m_LaplacianOperator;

  /** Constant for speckle scale function q0(t) */
  PixelType m_rho;
  /** Speckle coefficient of variation in the given image.
   *  Unity (1) assumed unless user specifies otherwise. */
  PixelType m_q0;

  /** Speckle scale function at time step t */
  PixelType m_q0_t;

  NeighborhoodSizeValueType m_Center;
  NeighborhoodSizeValueType m_Stride[ImageDimension];

  /** Elapsed iterations, with n=0 being the first iteration.
   *  Set from the parent image filter. */
  SizeValueType m_ElapsedIterationCount;

  /** Iterator for local diffusion coefficient computation */
  NeighborhoodType m_DiffusionNeighborhood;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkSpeckleReducingAnisotropicDiffusionFunction.hxx"
#endif

#endif // itkSpeckleReducingAnisotropicDiffusionFunction_h
