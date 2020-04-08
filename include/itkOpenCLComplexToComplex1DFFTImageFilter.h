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
#ifndef itkOpenCLComplexToComplex1DFFTImageFilter_h
#define itkOpenCLComplexToComplex1DFFTImageFilter_h

#include "itkComplexToComplex1DFFTImageFilter.h"

#define __CL_ENABLE_EXCEPTIONS
#include "CL/cl.hpp"
//#include "clFFT.h"

namespace itk
{
/** /class OpenCLComplexToComplex1DFFTImageFilter
 * /brief only do FFT along one dimension using OpenCL_FFT as a backend.
 *
 * The size of the image in the transformed direction must be a power of 2.
 *
 * There is considerable overhead to generate the FFT plan, which occurs
 * whenever the input image size changes.  Therefore, the throughput benefit
 * will only be realized for large images or many small images of
 * the same size.
 *
 * \ingroup
 */

template <typename TInputImage, typename TOutputImage = TInputImage>
class ITK_EXPORT OpenCLComplexToComplex1DFFTImageFilter :
    public ComplexToComplex1DFFTImageFilter<TInputImage, TOutputImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(OpenCLComplexToComplex1DFFTImageFilter);
  using TPixel = typename NumericTraits< typename TInputImage::PixelType >::ValueType;

  typedef OpenCLComplexToComplex1DFFTImageFilter Self;
  typedef ComplexToComplex1DFFTImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Standard class typedefs.*/
  typedef typename Superclass::InputImageType InputImageType;
  typedef typename Superclass::OutputImageType OutputImageType;

  struct OpenCLComplexType
  {
    TPixel real;
    TPixel imag;
  };

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(OpenCLComplexToComplex1DFFTImageFilter,
               FFT1DComplexToComplexImageFilter);


protected:
  OpenCLComplexToComplex1DFFTImageFilter();
  virtual ~OpenCLComplexToComplex1DFFTImageFilter()
  {
    if(m_PlanComputed)
      {
      //clFFT_DestroyPlan(this->m_Plan);
      delete [] this->m_InputBuffer;
      delete [] this->m_OutputBuffer;
      }
    delete m_clQueue;
    delete m_clContext;
  }

  virtual void GenerateData();  // generates output from input

  ///** Method to check if an array dimension is legal for current OpenCL FFT */
  bool Legaldim(int n); 
private:
  bool m_PlanComputed;
  //clFFT_Plan m_Plan;
  unsigned int m_LastImageSize;
  OpenCLComplexType *m_InputBuffer;
  OpenCLComplexType *m_OutputBuffer;
  cl::Context *m_clContext;
  cl::CommandQueue *m_clQueue;

};

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkOpenCLComplexToComplex1DFFTImageFilter.hxx"
#endif

#endif //itkOpenCLComplexToComplex1DFFTImageFilter_h
