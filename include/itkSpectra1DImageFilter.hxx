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
#ifndef itkSpectra1DImageFilter_hxx
#define itkSpectra1DImageFilter_hxx


#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageScanlineIterator.h"
#include "itkImageScanlineConstIterator.h"
#include "itkMetaDataObject.h"

#include "itkSpectra1DSupportWindowImageFilter.h"

namespace itk
{

template <typename TInputImage, typename TSupportWindowImage, typename TOutputImage>
Spectra1DImageFilter<TInputImage, TSupportWindowImage, TOutputImage>::Spectra1DImageFilter()
{
  this->AddRequiredInputName("SupportWindowImage");
  this->DynamicMultiThreadingOff();
}


template <typename TInputImage, typename TSupportWindowImage, typename TOutputImage>
void
Spectra1DImageFilter<TInputImage, TSupportWindowImage, TOutputImage>::GenerateOutputInformation()
{
  Superclass::GenerateOutputInformation();
  OutputImageType *              output = this->GetOutput();
  const SupportWindowImageType * supportWindowImage = this->GetSupportWindowImage();

  // Support windows inform output image region and metadata
  output->SetSpacing(supportWindowImage->GetSpacing());
  output->SetLargestPossibleRegion(supportWindowImage->GetLargestPossibleRegion());

  const MetaDataDictionary & dict = supportWindowImage->GetMetaDataDictionary();
  FFT1DSizeType              fft1DSize = 32;
  ExposeMetaData<FFT1DSizeType>(dict, "FFT1DSize", fft1DSize);

  // Number of frequency bins represented by each vector pixel.
  // Divide by two for Hermitian symmetry. Divide by two for Welch's method
  // with 50% overlap. Subtract one for discarding DC component.
  const FFT1DSizeType spectraComponents = fft1DSize / 2 / 2 - 1;

  output->SetVectorLength(spectraComponents);
}


template <typename TInputImage, typename TSupportWindowImage, typename TOutputImage>
void
Spectra1DImageFilter<TInputImage, TSupportWindowImage, TOutputImage>::BeforeThreadedGenerateData()
{
  const SupportWindowImageType * supportWindowImage = this->GetSupportWindowImage();
  const MetaDataDictionary &     dict = supportWindowImage->GetMetaDataDictionary();
  FFT1DSizeType                  fft1DSize = 32;
  ExposeMetaData<FFT1DSizeType>(dict, "FFT1DSize", fft1DSize);
  // Divide by two for Hermitian symmetry. Divide by two for Welch's method
  // with 50% overlap. Subtract one for discarding DC component.
  const FFT1DSizeType spectraComponents = fft1DSize / 2 / 2 - 1;

  const ThreadIdType numberOfWorkUnits = this->GetNumberOfWorkUnits();
  this->m_PerThreadDataContainer.resize(numberOfWorkUnits);
  for (ThreadIdType threadId = 0; threadId < numberOfWorkUnits; ++threadId)
  {
    PerThreadData & perThreadData = this->m_PerThreadDataContainer[threadId];
    perThreadData.ComplexVector.set_size(fft1DSize / 2);
    perThreadData.SpectraVector.resize(spectraComponents);
    perThreadData.LineImageRegionSize.Fill(1);
    perThreadData.LineImageRegionSize[0] = fft1DSize;
  }
}


template <typename TInputImage, typename TSupportWindowImage, typename TOutputImage>
void
Spectra1DImageFilter<TInputImage, TSupportWindowImage, TOutputImage>::VerifyInputInformation() const
{}


template <typename TInputImage, typename TSupportWindowImage, typename TOutputImage>
void
Spectra1DImageFilter<TInputImage, TSupportWindowImage, TOutputImage>::AddLineWindow(FFT1DSizeType       length,
                                                                                    LineWindowMapType & lineWindowMap)
{
  if (lineWindowMap.count(length) == 1)
  {
    return;
  }
  // Currently using a Hamming Window
  SpectraVectorType window(length);
  ScalarType        sum = NumericTraits<ScalarType>::ZeroValue();
  for (FFT1DSizeType sample = 0; sample < length; ++sample)
  {
    window[sample] = 0.54 + 0.46 * std::cos((Math::twopi * sample) / (length - 1));
    sum += window[sample];
  }
  for (FFT1DSizeType sample = 0; sample < length; ++sample)
  {
    window[sample] /= sum;
  }
  lineWindowMap[length] = window;
}


template <typename TInputImage, typename TSupportWindowImage, typename TOutputImage>
void
Spectra1DImageFilter<TInputImage, TSupportWindowImage, TOutputImage>::ComputeSpectra(const IndexType & lineIndex,
                                                                                     ThreadIdType      threadId,
                                                                                     SpectraLineType & spectraLine)
{
  const InputImageType * input = this->GetInput();
  PerThreadData &        perThreadData = this->m_PerThreadDataContainer[threadId];

  const FFT1DSizeType fftSize = static_cast<FFT1DSizeType>(perThreadData.ComplexVector.size());

  const typename InputImageType::RegionType lineRegion(lineIndex, perThreadData.LineImageRegionSize);
  InputImageIteratorType                    inputIt(input, lineRegion);
  inputIt.GoToBegin();
  perThreadData.ComplexVector.fill(0);
  typename ComplexVectorType::iterator       complexVectorIt = perThreadData.ComplexVector.begin();
  const typename ComplexVectorType::iterator complexVectorEnd = perThreadData.ComplexVector.end();
  typename SpectraVectorType::const_iterator windowIt = perThreadData.LineWindowMap[fftSize].begin();
  typename ComplexVectorType::const_iterator complexVectorConstIt = perThreadData.ComplexVector.begin();
  typename SpectraVectorType::iterator       spectraVectorIt = perThreadData.SpectraVector.begin();
  const size_t                               highFreq = perThreadData.SpectraVector.size();
  for (size_t freq = 0; freq < highFreq; ++freq)
  {
    spectraVectorIt[freq] = 0.0f;
  }
  const double overlap = 0.5;
  IndexType    segmentIndex(lineIndex);
  const double spectralScale = 1.0 / (fftSize * fftSize);
  // 3 segments for Welch's method with 50% overlap (D=M/2)
  for (unsigned int segment = 0; segment < 3; ++segment)
  {
    segmentIndex[0] =
      static_cast<IndexValueType>(lineIndex[0] + segment * perThreadData.LineImageRegionSize[0] * overlap / 3.0);
    inputIt.SetIndex(segmentIndex);
    complexVectorIt = perThreadData.ComplexVector.begin();
    windowIt = perThreadData.LineWindowMap[fftSize].begin();
    while (complexVectorIt != complexVectorEnd)
    {
      *complexVectorIt = inputIt.Value() * *windowIt;
      ++inputIt;
      ++complexVectorIt;
      ++windowIt;
    }
    FFT1DType fft1D(fftSize);
    fft1D.bwd_transform(perThreadData.ComplexVector);
    complexVectorConstIt = perThreadData.ComplexVector.begin();
    spectraVectorIt = perThreadData.SpectraVector.begin();
    // drop DC component
    ++complexVectorConstIt;

    // Each spectral component = (Re^2 + Im^2) / 3 / (fftSize)^2
    for (size_t freq = 0; freq < highFreq; ++freq)
    {
      spectraVectorIt[freq] +=
        std::real(*complexVectorConstIt * std::conj(*complexVectorConstIt)) / 3.0 * spectralScale;
      ++complexVectorConstIt;
    }
  }

  spectraLine.first = lineIndex;
  spectraLine.second = perThreadData.SpectraVector;
}


template <typename TInputImage, typename TSupportWindowImage, typename TOutputImage>
void
Spectra1DImageFilter<TInputImage, TSupportWindowImage, TOutputImage>::ThreadedGenerateData(
  const OutputImageRegionType & outputRegionForThread,
  ThreadIdType                  threadId)
{
  OutputImageType *              output = this->GetOutput();
  const SupportWindowImageType * supportWindowImage = this->GetSupportWindowImage();

  using OutputIteratorType = ImageLinearIteratorWithIndex<OutputImageType>;
  OutputIteratorType outputIt(output, outputRegionForThread);
  outputIt.SetDirection(1);

  const MetaDataDictionary & dict = supportWindowImage->GetMetaDataDictionary();
  PerThreadData &            perThreadData = this->m_PerThreadDataContainer[threadId];
  this->AddLineWindow(perThreadData.ComplexVector.size(), perThreadData.LineWindowMap);

  SpectraLinesContainerType spectraLines;

  using SupportWindowIteratorType = ImageLinearConstIteratorWithIndex<SupportWindowImageType>;
  SupportWindowIteratorType supportWindowIt(supportWindowImage, outputRegionForThread);
  supportWindowIt.SetDirection(1);

  SpectraLineType spectraLine;
  for (outputIt.GoToBegin(), supportWindowIt.GoToBegin(); !outputIt.IsAtEnd();
       outputIt.NextLine(), supportWindowIt.NextLine())
  {
    spectraLines.clear();
    while (!outputIt.IsAtEndOfLine())
    {
      // Compute the per line spectra.
      const SupportWindowType & supportWindow = supportWindowIt.Value();
      if (spectraLines.size() == 0) // first window in this lateral direction
      {
        const typename SupportWindowType::const_iterator windowLineEnd = supportWindow.end();
        for (typename SupportWindowType::const_iterator windowLine = supportWindow.begin(); windowLine != windowLineEnd;
             ++windowLine)
        {
          const IndexType & lineIndex = *windowLine;
          this->ComputeSpectra(lineIndex, threadId, spectraLine);
          spectraLines.push_back(spectraLine);
        }
      }
      else // subsequent window along a line
      {
        const IndexValueType desiredFirstLine = supportWindow.front()[1];
        while (spectraLines.front().first[1] < desiredFirstLine)
        {
          spectraLines.pop_front();
        }
        const typename SupportWindowType::const_iterator   windowLineEnd = supportWindow.end();
        typename SpectraLinesContainerType::iterator       spectraLinesIt = spectraLines.begin();
        const typename SpectraLinesContainerType::iterator spectraLinesEnd = spectraLines.end();
        for (typename SupportWindowType::const_iterator windowLine = supportWindow.begin(); windowLine != windowLineEnd;
             ++windowLine)
        {
          const IndexType & lineIndex = *windowLine;
          if (spectraLinesIt == spectraLinesEnd) // past the end of the previously processed lines
          {
            this->ComputeSpectra(lineIndex, threadId, spectraLine);
            spectraLines.push_back(spectraLine);
          }
          else if (lineIndex[1] == (spectraLinesIt->first)[1]) // one of the same lines that was previously computed
          {
            if (lineIndex[0] != (spectraLinesIt->first)[0])
            {
              this->ComputeSpectra(lineIndex, threadId, spectraLine);
              *spectraLinesIt = spectraLine;
            }
            ++spectraLinesIt;
          }
          else
          {
            itkExceptionMacro("Unexpected line");
          }
        }
      }

      // lateral window and sum
      const size_t spectraLinesCount = spectraLines.size();
      this->AddLineWindow(spectraLinesCount, perThreadData.LineWindowMap);
      typename OutputImageType::PixelType outputPixel;
      const FFT1DSizeType                 spectralComponents = perThreadData.SpectraVector.size();
      outputPixel.SetSize(spectralComponents);
      outputPixel.Fill(NumericTraits<ScalarType>::ZeroValue());
      typename SpectraVectorType::const_iterator   windowIt = perThreadData.LineWindowMap[spectraLinesCount].begin();
      typename SpectraLinesContainerType::iterator linesIt = spectraLines.begin();
      for (size_t line = 0; line < spectraLinesCount; ++line)
      {
        typename SpectraVectorType::const_iterator spectraIt = linesIt->second.begin();
        for (FFT1DSizeType sample = 0; sample < spectralComponents; ++sample)
        {
          outputPixel[sample] += *windowIt * *spectraIt;
          ++spectraIt;
        }
        ++windowIt;
        ++linesIt;
      }
      outputIt.Set(outputPixel);

      ++outputIt;
      ++supportWindowIt;
    }
  }

  // Optionally normalize for system noise via reference spectra image input
  const OutputImageType * referenceSpectra = this->GetReferenceSpectraImage();
  if (referenceSpectra != nullptr)
  {
    using ReferenceSpectraIteratorType = ImageScanlineConstIterator<OutputImageType>;
    ReferenceSpectraIteratorType referenceSpectraIt(referenceSpectra, outputRegionForThread);

    using PopulatedOutputIteratorType = ImageScanlineIterator<OutputImageType>;
    PopulatedOutputIteratorType populatedOutputIt(output, outputRegionForThread);

    const unsigned int numberOfComponents = referenceSpectra->GetNumberOfComponentsPerPixel();
    if (numberOfComponents != output->GetNumberOfComponentsPerPixel())
    {
      itkExceptionMacro("ReferenceSpectraImage has " << numberOfComponents << " while the output image has "
                                                     << output->GetNumberOfComponentsPerPixel() << " components");
    }

    for (referenceSpectraIt.GoToBegin(), populatedOutputIt.GoToBegin(); !populatedOutputIt.IsAtEnd();)
    {
      while (!populatedOutputIt.IsAtEndOfLine())
      {
        using PixelType = typename OutputImageType::PixelType;
        PixelType       outputPixel = populatedOutputIt.Get();
        const PixelType referencePixel = referenceSpectraIt.Get();
        for (unsigned int component = 0; component < numberOfComponents; ++component)
        {
          if (Math::FloatAlmostEqual(referencePixel[component], NumericTraits<typename PixelType::ValueType>::Zero))
          {
            outputPixel[component] = NumericTraits<typename PixelType::ValueType>::Zero;
          }
          else
          {
            outputPixel[component] /= referencePixel[component];
          }
        }
        populatedOutputIt.Set(outputPixel);

        ++populatedOutputIt;
        ++referenceSpectraIt;
      }
      populatedOutputIt.NextLine();
      referenceSpectraIt.NextLine();
    }
  }
}


} // end namespace itk

#endif
