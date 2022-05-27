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
#ifndef itkSpeckleReducingAnisotropicDiffusionImageFilter_h
#define itkSpeckleReducingAnisotropicDiffusionImageFilter_h

#include "itkAnisotropicDiffusionImageFilter.h"
#include "itkSpeckleReducingAnisotropicDiffusionFunction.h"

namespace itk
{
/** \class SpeckleReducingAnisotropicDiffusionImageFilter
 * \brief This filter performs Speckle Reducing Anisotropic Diffusion (SRAD)
 * on a scalar itk::Image using the updated Yu-Acton gradient magnitude based
 * equation.
 *
 * For detailed information on anisotropic diffusion, see
 * itkAnisotropicDiffusionFunction and
 * SpeckleReducingAnisotropicDiffusionFunction.
 *
 * \par Inputs and Outputs
 * The input to this filter should be a 2D padded scalar itk::Image.
 * The output image will be a diffused copy of the input.

 * \par Parameters
 * Please see the description of parameters given in
 * itkAnisotropicDiffusionImageFilter.
 *
 * \sa AnisotropicDiffusionImageFilter
 * \sa AnisotropicDiffusionFunction
 * \sa SpeckleReducingAnisotropicDiffusionFunction
 * \ingroup ImageEnhancement
 * \ingroup ImageFilters
 * \ingroup ITKAnisotropicSmoothing
 * \ingroup ITKUltrasound
 */
template <typename TInputImage, typename TOutputImage = TInputImage>
class ITK_TEMPLATE_EXPORT SpeckleReducingAnisotropicDiffusionImageFilter
  : public AnisotropicDiffusionImageFilter<TInputImage, TOutputImage>
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(SpeckleReducingAnisotropicDiffusionImageFilter);

  /** Standard class type aliases. */
  using Self = SpeckleReducingAnisotropicDiffusionImageFilter;
  using Superclass = AnisotropicDiffusionImageFilter<TInputImage, TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Standard method for creation through object factory. */
  itkNewMacro(Self);

  /** Run-time class information. */
  itkTypeMacro(SpeckleReducingAnisotropicDiffusionImageFilter, AnisotropicDiffusionImageFilter);

  /** Extract information from the superclass. */
  using UpdateBufferType = typename Superclass::UpdateBufferType;

  /** Extract information from the superclass. */
  static constexpr unsigned int ImageDimension = Superclass::ImageDimension;

#ifdef ITK_USE_CONCEPT_CHECKING
  // Begin concept checking
  itkConceptMacro(UpdateBufferHasNumericTraitsCheck, (Concept::HasNumericTraits<typename UpdateBufferType::PixelType>));
  // End concept checking
#endif

  itkSetMacro(Rho, double);
  itkGetConstMacro(Rho, double);

  itkSetMacro(Q0, double);
  itkGetConstMacro(Q0, double);

protected:
  SpeckleReducingAnisotropicDiffusionImageFilter()
    : m_Q0(1.0)
    , m_Rho(0.2)
  {
    typename SpeckleReducingAnisotropicDiffusionFunction<UpdateBufferType>::Pointer p =
      SpeckleReducingAnisotropicDiffusionFunction<UpdateBufferType>::New();
    this->SetDifferenceFunction(p);
  }

  ~SpeckleReducingAnisotropicDiffusionImageFilter() override = default;

  /** Prepare for the iteration process. */
  void
  InitializeIteration() override
  {
    auto * f = dynamic_cast<SpeckleReducingAnisotropicDiffusionFunction<UpdateBufferType> *>(
      this->GetDifferenceFunction().GetPointer());
    if (!f)
    {
      throw ExceptionObject(
        __FILE__, __LINE__, "Speckle reducing anisotropic diffusion function is not set.", ITK_LOCATION);
    }

    f->SetRho(this->m_Rho);
    f->SetQ0(this->m_Q0);
    f->SetElapsedIterationCount(this->m_ElapsedIterations);

    Superclass::InitializeIteration();
  }

  /** Constant for speckle scale function q0(t) */
  double m_Rho;
  /** Speckle coefficient of variation in the given image.
   *  Unity (1) assumed unless user specifies otherwise. */
  double m_Q0;
};
} // namespace itk

#endif
