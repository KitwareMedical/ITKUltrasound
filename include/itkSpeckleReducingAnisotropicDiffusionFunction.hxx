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
#ifndef itkSpeckleReducingAnisotropicDiffusionFunction_hxx
#define itkSpeckleReducingAnisotropicDiffusionFunction_hxx


#include <cmath>
#include "itkMath.h"
#include "itkSize.h"
#include "itkNumericTraits.h"

namespace itk
{

template <typename TImage>
SpeckleReducingAnisotropicDiffusionFunction<TImage>::SpeckleReducingAnisotropicDiffusionFunction()
  : m_rho(0.2)
  , m_q0(1.0)
  , m_ElapsedIterationCount(0)
{
  unsigned int i;
  RadiusType   r;

  // Set up neighborhood iterator parameters such that we use a
  // 3x3 window to calculate the diffusion coefficient at each pixel.
  for (i = 0; i < ImageDimension; ++i)
  {
    r[i] = 1;
  }
  this->SetRadius(r);

  // Dummy neighborhood used to set up the slices.
  Neighborhood<PixelType, ImageDimension> it;
  it.SetRadius(r);

  // Slice the neighborhood
  m_Center = it.Size() / 2;

  for (i = 0; i < ImageDimension; ++i)
  {
    m_Stride[i] = it.GetStride(i);
  }

  // Set slice to access the local neighborhood for computing coefficients.
  // We use 3x3 gradient and laplacian operators.
  for (i = 0; i < ImageDimension; ++i)
  {
    // Contains positional information to slice a given neighborhood along axis i
    // for derivative and laplacian computations.
    x_slice[i] = std::slice(m_Center - m_Stride[i], 3, m_Stride[i]);
  }

  // Allocate the derivative and laplacian operators.
  m_DerivativeOperator.SetDirection(0); // Not relevant, will be applied in a slice-based
                                        // fashion.
  m_DerivativeOperator.SetOrder(1);
  m_DerivativeOperator.CreateDirectional();
  m_LaplacianOperator.CreateToRadius(1);
};

template <typename TImage>
void
SpeckleReducingAnisotropicDiffusionFunction<TImage>::InitializeIteration()
{
  m_LaplacianOperator.SetDerivativeScalings(this->m_ScaleCoefficients);

  // Estimate q0(t) according to Yu,Acton (37)
  // where t := (n+1) * delta_t
  m_q0_t = static_cast<PixelType>(m_q0 * std::exp(-1.0f * m_rho * this->GetTimeStep() * (m_ElapsedIterationCount + 1)));
}

template <typename TImage>
typename SpeckleReducingAnisotropicDiffusionFunction<TImage>::PixelType
SpeckleReducingAnisotropicDiffusionFunction<TImage>::ComputeDiffusionCoefficient(PixelType intensity,
                                                                                 PixelType averageGradient,
                                                                                 PixelType laplacian)
{
  // Epsilon to avoid divide-by-zero errors
  const float eps = itk::Math::eps;

  // Compute instantaneous coefficient of variation.
  // See (35) in Yu, Acton.
  PixelType numerator =
    0.5 * itk::Math::sqr(averageGradient / (intensity + eps)) - itk::Math::sqr(0.25 * laplacian / (intensity + eps));
  // We only use the real part of q(t) so we discard any values n < 0
  // that would result in q(t) = sqrt(n/d) = a + b i where b != 0
  if (numerator < 0)
  {
    numerator = 0;
  }

  PixelType denominator = itk::Math::sqr(1 + 0.25 * laplacian / (intensity + eps));

  // Instantaneous coefficient of variation
  PixelType q_t = std::sqrt(numerator / denominator);

  // Compute diffusion coefficient.
  // See (33) in Yu, Acton.
  PixelRealType C = 1 / (1 + (itk::Math::sqr(q_t) - itk::Math::sqr(m_q0_t)) /
                               (itk::Math::sqr(m_q0_t) * (1 + itk::Math::sqr(m_q0_t)) + eps));
  return C;
}

template <typename TImage>
typename SpeckleReducingAnisotropicDiffusionFunction<TImage>::PixelType
SpeckleReducingAnisotropicDiffusionFunction<TImage>::ComputeUpdate(const NeighborhoodType & it,
                                                                   void *,
                                                                   const FloatOffsetType &)
{
  // Ignore edge pixels. We assume the user has mirror-padded the input image.
  if (!it.InBounds())
    return 0;

  unsigned int i;

  const size_t               NUM_C_TERMS = 3; // c[i,j], c[i,j+1], c[i+1,j]
  std::vector<PixelRealType> c;               // pixel-wise diffusion coefficients

  PixelRealType dx;
  PixelRealType dx_avg;  // local derivative
  PixelRealType dx2_avg; // local laplacian
  PixelRealType delta;   // calculated intensity update

  // Create an iterator to walk over a local 2x2 region to generate
  // diffusion coefficient statistics at each pixel
  RadiusType radius = RadiusType();
  radius[0] = 1;
  radius[1] = 1;
  RadiusType regionSize = RadiusType();
  regionSize[0] = 2;
  regionSize[1] = 2;
  RegionType region;
  region.SetSize(regionSize);
  region.SetIndex(it.GetIndex());
  NeighborhoodType coeffIt(radius, it.GetImagePointer(), region);

  while (!coeffIt.IsAtEnd() && c.size() < NUM_C_TERMS)
  {
    // Try to get the diffusion coefficient from the cache
    // Calculate the centralized derivatives for each dimension.
    dx = 0;
    dx_avg = 0;
    for (i = 0; i < ImageDimension; ++i)
    {
      dx = m_InnerProduct(x_slice[i], coeffIt, m_DerivativeOperator);
      dx *= this->m_ScaleCoefficients[i];
      dx_avg += itk::Math::sqr(dx);
    }
    dx_avg = std::sqrt(dx_avg);

    // Laplacian
    dx2_avg = m_InnerProduct(coeffIt, m_LaplacianOperator);

    c.push_back(ComputeDiffusionCoefficient(coeffIt.GetPixel(m_Center), dx_avg, dx2_avg));
    ++coeffIt;
  }
  // Check that we got values for c[i,j]; c[i,j+1]; and c[i+1,j]
  itkAssertInDebugAndIgnoreInReleaseMacro(c.size() == NUM_C_TERMS);

  // Compute divergence (see (58) in Yu, Acton).
  delta = NumericTraits<PixelRealType>::ZeroValue();
  delta += c[2] * static_cast<PixelRealType>(it.GetPixel(m_Center + m_Stride[1]) - it.GetPixel(m_Center));
  delta += c[0] * static_cast<PixelRealType>(it.GetPixel(m_Center - m_Stride[1]) - it.GetPixel(m_Center));
  delta += c[1] * static_cast<PixelRealType>(it.GetPixel(m_Center + m_Stride[0]) - it.GetPixel(m_Center));
  delta += c[0] * static_cast<PixelRealType>(it.GetPixel(m_Center - m_Stride[0]) - it.GetPixel(m_Center));
  for (auto scale : this->m_ScaleCoefficients)
  {
    delta *= scale;
  }

  PixelType pixel_delta = static_cast<PixelType>(delta / 4);
  return pixel_delta;
}
} // end namespace itk

#endif
