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
#ifndef itkOpenCLForward1DFFTImageFilter_hxx
#define itkOpenCLForward1DFFTImageFilter_hxx

#include "itkForward1DFFTImageFilter.hxx"
#include "itkOpenCLForward1DFFTImageFilter.h"
#include "itkclFFTInitializer.h"

#include <vector>

#include "itkIndent.h"
#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkMetaDataObject.h"

namespace itk
{

template <typename TInputImage, typename TOutputImage>
OpenCLForward1DFFTImageFilter<TInputImage, TOutputImage>::OpenCLForward1DFFTImageFilter()
{
  try
  {
    auto initObject = clFFFInitialization();
    m_clContext = new cl::Context(CL_DEVICE_TYPE_ALL);
    std::vector<cl::Device> devices = m_clContext->getInfo<CL_CONTEXT_DEVICES>();
    if (devices.size() < 1)
    {
      itkExceptionMacro("No OpenCL devices found.");
    }
    // @todo: code to select the fastest device, or the device that is
    // CL_DEVICE_TYPE_ACCELERATOR
    this->m_clQueue = new cl::CommandQueue(*m_clContext, devices[0]);
  }
  catch (const cl::Error & e)
  {
    itkExceptionMacro("Error in OpenCL: " << e.what() << "(" << e.err() << ")");
  }
}

template <typename TInputImage, typename TOutputImage>
bool
OpenCLForward1DFFTImageFilter<TInputImage, TOutputImage>::Legaldim(int n)
{
  return clFFTFactorization(n);
}

template <typename TInputImage, typename TOutputImage>
void
OpenCLForward1DFFTImageFilter<TInputImage, TOutputImage>::GenerateData()
{
  // get pointers to the input and output
  typename InputImageType::ConstPointer inputPtr = this->GetInput();
  typename OutputImageType::Pointer     outputPtr = this->GetOutput();

  if (!inputPtr || !outputPtr)
  {
    return;
  }

  // allocate output buffer memory
  outputPtr->SetBufferedRegion(outputPtr->GetRequestedRegion());
  outputPtr->Allocate();

  const typename InputImageType::SizeType &  inputSize = inputPtr->GetRequestedRegion().GetSize();
  const typename OutputImageType::SizeType & outputSize = outputPtr->GetRequestedRegion().GetSize();

  unsigned int vec_size = inputSize[this->GetDirection()];
  if (!this->Legaldim(vec_size))
  {
    ExceptionObject exception(__FILE__, __LINE__);
    exception.SetDescription("Illegal dimension for FFT: " + std::to_string(vec_size));
    exception.SetLocation(ITK_LOCATION);
    throw exception;
  }

  cl_int batchSize = 1;
  for (unsigned int i = 0; i < TInputImage::ImageDimension; i++)
  {
    batchSize *= outputSize[i];
  }
  unsigned int totalSize = batchSize;
  batchSize /= outputSize[this->GetDirection()];
  cl_command_queue queue = (*m_clQueue)();

  if (this->m_PlanComputed) // if we've already computed a plan
  {
    // if the image sizes aren't the same,
    // we have to compute the plan again
    if (this->m_LastImageSize != totalSize)
    {
      delete[] this->m_InputBuffer;
      delete[] this->m_OutputBuffer;
      clfftDestroyPlan(&this->m_Plan);
      this->m_PlanComputed = false;
    }
  }
  if (!this->m_PlanComputed)
  {
    try
    {
      this->m_InputBuffer = new OpenCLComplexType[totalSize];
      this->m_OutputBuffer = new OpenCLComplexType[totalSize];
    }
    catch (std::bad_alloc &)
    {
      itkExceptionMacro("Problem allocating memory for internal computations");
    }

    this->m_LastImageSize = totalSize;
    const size_t n[3] = { inputSize[this->GetDirection()], 1, 1 };
    clfftStatus  error_code = clfftCreateDefaultPlan(&this->m_Plan, (*m_clContext)(), CLFFT_1D, n);
    if (!this->m_Plan || error_code)
    {
      itkExceptionMacro("Could not create OpenCL FFT Plan.");
    }

    error_code = clfftSetResultLocation(this->m_Plan, CLFFT_INPLACE);
    error_code = clfftSetPlanBatchSize(this->m_Plan, batchSize);
    if (std::is_same<TPixel, double>::value) // float by default
    {
      error_code = clfftSetPlanPrecision(this->m_Plan, CLFFT_DOUBLE);
    }
    clfftBakePlan(this->m_Plan, 1, &queue, nullptr, nullptr);
    this->m_PlanComputed = true;
  }

  using InputIteratorType = itk::ImageLinearConstIteratorWithIndex<InputImageType>;
  using OutputIteratorType = itk::ImageLinearIteratorWithIndex<OutputImageType>;
  InputIteratorType  inputIt(inputPtr, inputPtr->GetRequestedRegion());
  OutputIteratorType outputIt(outputPtr, outputPtr->GetRequestedRegion());

  inputIt.SetDirection(this->GetDirection());
  outputIt.SetDirection(this->GetDirection());

  OpenCLComplexType * inputBufferIt = this->m_InputBuffer;
  // for every fft line
  for (inputIt.GoToBegin(); !inputIt.IsAtEnd(); inputIt.NextLine())
  {
    // copy the input line into our buffer
    inputIt.GoToBeginOfLine();
    while (!inputIt.IsAtEndOfLine())
    {
      inputBufferIt->real(inputIt.Get());
      ++inputIt;
      ++inputBufferIt;
    }
  }

  try
  {
    // do the transform
    cl::Buffer clDataBuffer(
      *m_clContext, CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR, totalSize * sizeof(TPixel) * 2, m_InputBuffer);
    cl_mem      clPointer = clDataBuffer();
    clfftStatus err =
      clfftEnqueueTransform(this->m_Plan, CLFFT_FORWARD, 1, &queue, 0, nullptr, nullptr, &clPointer, nullptr, nullptr);
    if (err)
    {
      itkExceptionMacro("Error in clfftEnqueueTransform(" << err << ")");
    }

    cl_int err2 =
      m_clQueue->enqueueReadBuffer(clDataBuffer, CL_TRUE, 0, totalSize * sizeof(TPixel) * 2, m_OutputBuffer);
  }
  catch (const cl::Error & e)
  {
    itkExceptionMacro("Error in OpenCL: " << e.what() << "(" << e.err() << ")");
  }

  OpenCLComplexType * outputBufferIt = this->m_OutputBuffer;
  for (outputIt.GoToBegin(); !outputIt.IsAtEnd(); outputIt.NextLine())
  {
    outputIt.GoToBeginOfLine();
    while (!outputIt.IsAtEndOfLine())
    {
      outputIt.Set(*reinterpret_cast<typename OutputIteratorType::PixelType *>(outputBufferIt));
      ++outputIt;
      ++outputBufferIt;
    }
  }
}

} // namespace itk

#endif // itkOpenCLForward1DFFTImageFilter_hxx
