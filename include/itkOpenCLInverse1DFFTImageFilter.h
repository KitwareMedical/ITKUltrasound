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
#if !defined(itkOpenCLInverse1DFFTImageFilter_h) && defined(ITKUltrasound_USE_clFFT)
#  define itkOpenCLInverse1DFFTImageFilter_h

#  include "itkInverse1DFFTImageFilter.h"

#  define __CL_ENABLE_EXCEPTIONS
#  include "CL/cl.hpp"
#  include "clFFT.h"

namespace itk
{
/** \class OpenCLInverse1DFFTImageFilter
 * \brief Do FFT along only one dimension using clFFT library as a backend.
 *
 * The size of the image in the transformed direction
 * must be a multiple of powers of 2, 3, 5, and 7.
 *
 * There is considerable overhead to generate the FFT plan, which occurs
 * whenever the input image size changes.  Therefore, the throughput benefit
 * will only be realized for large images or many small images of
 * the same size.
 *
 * \ingroup FourierTransform
 * \ingroup Ultrasound
 */

template <typename TInputImage,
          typename TOutputImage =
            Image<typename NumericTraits<typename TInputImage::PixelType>::ValueType, TInputImage::ImageDimension>>
class ITK_TEMPLATE_EXPORT OpenCLInverse1DFFTImageFilter : public Inverse1DFFTImageFilter<TInputImage, TOutputImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(OpenCLInverse1DFFTImageFilter);
  using TPixel = typename NumericTraits<typename TInputImage::PixelType>::ValueType;

  using Self = OpenCLInverse1DFFTImageFilter;
  using Superclass = Inverse1DFFTImageFilter<TInputImage, TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Standard class type alias.*/
  using InputImageType = typename Superclass::InputImageType;
  using OutputImageType = typename Superclass::OutputImageType;

  using OpenCLComplexType = std::complex<TPixel>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(OpenCLInverse1DFFTImageFilter, FFT1DComplexConjugateToRealImageFilter);

  SizeValueType
  GetSizeGreatestPrimeFactor() const override
  {
    return 7; // clFFT supports prime factors 2, 3, 5 and 7
  }

protected:
  OpenCLInverse1DFFTImageFilter();
  ~OpenCLInverse1DFFTImageFilter() override
  {
    if (m_PlanComputed)
    {
      clfftDestroyPlan(&this->m_Plan);
      delete[] this->m_InputBuffer;
      delete[] this->m_OutputBuffer;
    }
    delete m_clQueue;
    delete m_clContext;
  }

  virtual void
  GenerateData(); // generates output from input

  ///** Method to check if an array dimension is legal for current OpenCL FFT */
  bool
  Legaldim(int n);

private:
  bool                m_PlanComputed = false;
  clfftPlanHandle     m_Plan = 0;
  unsigned int        m_LastImageSize = 0;
  OpenCLComplexType * m_InputBuffer = nullptr;
  OpenCLComplexType * m_OutputBuffer = nullptr;
  cl::Context *       m_clContext = nullptr;
  cl::CommandQueue *  m_clQueue = nullptr;
};

} // namespace itk

#  ifndef ITK_MANUAL_INSTANTIATION
#    include "itkOpenCLInverse1DFFTImageFilter.hxx"
#  endif

#endif // itkOpenCLInverse1DFFTImageFilter_h
