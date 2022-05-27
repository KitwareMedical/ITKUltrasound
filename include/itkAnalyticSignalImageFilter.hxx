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
#ifndef itkAnalyticSignalImageFilter_hxx
#define itkAnalyticSignalImageFilter_hxx


#include "itkVnlForward1DFFTImageFilter.h"
#include "itkVnlComplexToComplex1DFFTImageFilter.h"

#if defined(ITK_USE_FFTWD) || defined(ITK_USE_FFTWF)
#  include "itkFFTWForward1DFFTImageFilter.h"
#  include "itkFFTWComplexToComplex1DFFTImageFilter.h"
#endif

#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkMetaDataObject.h"

namespace itk
{

template <typename TInputImage, typename TOutputImage>
AnalyticSignalImageFilter<TInputImage, TOutputImage>::AnalyticSignalImageFilter()
{
  m_FFTRealToComplexFilter = FFTRealToComplexType::New();
  m_FFTComplexToComplexFilter = FFTComplexToComplexType::New();
  m_FFTComplexToComplexFilter->SetTransformDirection(FFTComplexToComplexType::INVERSE);

  this->SetDirection(0);
}


template <typename TInputImage, typename TOutputImage>
void
AnalyticSignalImageFilter<TInputImage, TOutputImage>::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // get pointers to the inputs
  InputImageType *  input = const_cast<InputImageType *>(this->GetInput());
  OutputImageType * output = this->GetOutput();

  if (!input || !output)
  {
    return;
  }

  // we need to compute the input requested region (size and start index)
  using OutputSizeType = const typename OutputImageType::SizeType &;
  OutputSizeType outputRequestedRegionSize = output->GetRequestedRegion().GetSize();
  using OutputIndexType = const typename OutputImageType::IndexType &;
  OutputIndexType outputRequestedRegionStartIndex = output->GetRequestedRegion().GetIndex();

  //// the regions other than the fft direction are fine
  typename InputImageType::SizeType  inputRequestedRegionSize = outputRequestedRegionSize;
  typename InputImageType::IndexType inputRequestedRegionStartIndex = outputRequestedRegionStartIndex;

  // we but need all of the input in the fft direction
  const unsigned int                        direction = this->GetDirection();
  const typename InputImageType::SizeType & inputLargeSize = input->GetLargestPossibleRegion().GetSize();
  inputRequestedRegionSize[direction] = inputLargeSize[direction];
  const typename InputImageType::IndexType & inputLargeIndex = input->GetLargestPossibleRegion().GetIndex();
  inputRequestedRegionStartIndex[direction] = inputLargeIndex[direction];

  typename InputImageType::RegionType inputRequestedRegion;
  inputRequestedRegion.SetSize(inputRequestedRegionSize);
  inputRequestedRegion.SetIndex(inputRequestedRegionStartIndex);

  input->SetRequestedRegion(inputRequestedRegion);
}


template <typename TInputImage, typename TOutputImage>
void
AnalyticSignalImageFilter<TInputImage, TOutputImage>::EnlargeOutputRequestedRegion(DataObject * out)
{
  OutputImageType * output = dynamic_cast<OutputImageType *>(out);

  // we need to enlarge the region in the fft direction to the
  // largest possible in that direction
  using ConstOutputSizeType = const typename OutputImageType::SizeType &;
  ConstOutputSizeType requestedSize = output->GetRequestedRegion().GetSize();
  ConstOutputSizeType outputLargeSize = output->GetLargestPossibleRegion().GetSize();
  using ConstOutputIndexType = const typename OutputImageType::IndexType &;
  ConstOutputIndexType requestedIndex = output->GetRequestedRegion().GetIndex();
  ConstOutputIndexType outputLargeIndex = output->GetLargestPossibleRegion().GetIndex();

  typename OutputImageType::SizeType  enlargedSize = requestedSize;
  typename OutputImageType::IndexType enlargedIndex = requestedIndex;
  const unsigned int                  direction = this->GetDirection();
  enlargedSize[direction] = outputLargeSize[direction];
  enlargedIndex[direction] = outputLargeIndex[direction];

  typename OutputImageType::RegionType enlargedRegion;
  enlargedRegion.SetSize(enlargedSize);
  enlargedRegion.SetIndex(enlargedIndex);
  output->SetRequestedRegion(enlargedRegion);
}


template <typename TInputImage, typename TOutputImage>
void
AnalyticSignalImageFilter<TInputImage, TOutputImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  const unsigned int direction = this->GetDirection();
  os << indent << "Direction: " << direction << std::endl;

  os << indent << "FFTRealToComplexFilter: " << std::endl;
  m_FFTRealToComplexFilter->Print(os, indent);
  if (m_FrequencyFilter.IsNotNull())
  {
    os << indent << "FrequencyFilter: " << std::endl;
    m_FrequencyFilter->Print(os, indent);
  }
  os << indent << "FFTComplexToComplexFilter: " << std::endl;
  m_FFTComplexToComplexFilter->Print(os, indent);
}


template <typename TInputImage, typename TOutputImage>
void
AnalyticSignalImageFilter<TInputImage, TOutputImage>::GenerateData()
{
  this->AllocateOutputs();

  OutputImageType * output = this->GetOutput();

  m_FFTRealToComplexFilter->SetInput(this->GetInput());
  if (m_FrequencyFilter.IsNotNull())
  {
    m_FrequencyFilter->SetInput(m_FFTRealToComplexFilter->GetOutput());
    m_FrequencyFilter->GetOutput()->SetRequestedRegion(output->GetRequestedRegion());
    m_FrequencyFilter->GetOutput()->SetLargestPossibleRegion(output->GetLargestPossibleRegion());
    m_FrequencyFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    m_FrequencyFilter->Update();
  }
  else
  {
    m_FFTRealToComplexFilter->GetOutput()->SetRequestedRegion(output->GetRequestedRegion());
    m_FFTRealToComplexFilter->GetOutput()->SetLargestPossibleRegion(output->GetLargestPossibleRegion());
    m_FFTRealToComplexFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
    m_FFTRealToComplexFilter->Update();
  }

  // get pointers to the input and output
  const typename FFTRealToComplexType::OutputImageType * input;
  if (m_FrequencyFilter.IsNotNull())
  {
    input = m_FrequencyFilter->GetOutput();
  }
  else
  {
    input = m_FFTRealToComplexFilter->GetOutput();
  }

  const typename FFTRealToComplexType::OutputImageType::SizeType & inputSize = input->GetRequestedRegion().GetSize();
  const unsigned int                                               direction = this->GetDirection();
  const unsigned int                                               size = inputSize[direction];
  unsigned int                                                     dubSize;
  bool                                                             even;
  if (size % 2 == 0)
  {
    even = true;
    dubSize = size / 2 - 1;
  }
  else
  {
    even = false;
    dubSize = (size + 1) / 2 - 1;
  }

  MultiThreaderBase * multiThreader = this->GetMultiThreader();
  multiThreader->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
  multiThreader->template ParallelizeImageRegionRestrictDirection<ImageDimension>(
    direction,
    output->GetRequestedRegion(),
    [this, dubSize, even, input](const typename OutputImageType::RegionType & lambdaRegion) {
      OutputImageType *  output = this->GetOutput();
      const unsigned int direction = this->GetDirection();

      using InputIteratorType = ImageLinearConstIteratorWithIndex<typename FFTRealToComplexType::OutputImageType>;
      using OutputIteratorType = ImageLinearIteratorWithIndex<OutputImageType>;
      InputIteratorType  inputIt(input, lambdaRegion);
      OutputIteratorType outputIt(output, lambdaRegion);
      inputIt.SetDirection(direction);
      outputIt.SetDirection(direction);

      // for every fft line
      for (inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd(); outputIt.NextLine(), inputIt.NextLine())
      {
        inputIt.GoToBeginOfLine();
        outputIt.GoToBeginOfLine();

        // DC
        outputIt.Set(inputIt.Get());
        ++inputIt;
        ++outputIt;
        for (unsigned int i = 0; i < dubSize; i++)
        {
          outputIt.Set(inputIt.Get() * static_cast<typename TInputImage::PixelType>(2));
          ++outputIt;
          ++inputIt;
        }
        if (even)
        {
          outputIt.Set(inputIt.Get());
          ++inputIt;
          ++outputIt;
        }
        while (!outputIt.IsAtEndOfLine())
        {
          outputIt.Set(static_cast<typename TInputImage::PixelType>(0));
          ++outputIt;
        }
      }
    },
    this);

  // Trippy, eh?
  m_FFTComplexToComplexFilter->SetInput(output);
  m_FFTComplexToComplexFilter->GetOutput()->SetRequestedRegion(output->GetRequestedRegion());
  m_FFTComplexToComplexFilter->GetOutput()->SetLargestPossibleRegion(output->GetLargestPossibleRegion());
  m_FFTComplexToComplexFilter->SetNumberOfWorkUnits(this->GetNumberOfWorkUnits());
  m_FFTComplexToComplexFilter->Update();
  this->GraftOutput(m_FFTComplexToComplexFilter->GetOutput());
}

} // end namespace itk

#endif // itkAnalyticSignalImageFilter_hxx
