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
#ifndef itkSpectra1DImageFilter_h
#define itkSpectra1DImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkDefaultConvertPixelTraits.h"
#include "itkImageRegionConstIterator.h"

#include "vnl/algo/vnl_fft_base.h"
#include "vnl/algo/vnl_fft_1d.h"

#include <utility>

#include <unordered_map>

#include "itkSpectra1DSupportWindowImageFilter.h"

namespace itk
{

/** \class Spectra1DImageFilter
 * \brief Generate an image of local spectra.
 *
 * This image filter takes in the input image and the support window image
 * that has indexes of the local lines used to compute the local spectra.
 * The output vector image contains spectral content information for each
 * support window in the given input.
 *
 * Each pixel in the output image represents the frequency content for the
 * given input support window, which encompasses a subset of the given beam line
 * averaged with adjacent local beam lines. Pixel vector values correspond to frequency
 * bins on the range (0,nyquist] with DC content discarded.
 *
 * FFT is computed in each support window using a Hamming window. A reference spectra
 * image may be provided to compensate for system noise.
 *
 * This filter expects that beam input lies along the zeroth dimension and lateral lines
 * lie along the first dimension. Images not matching this description may be permuted
 * with itk::PermuteAxesImageFilter prior to running the filter.
 *
 * \ingroup Ultrasound
 *
 * \sa Spectra1DSupportWindowImageFilter
 * \sa PermuteAxesImageFilter
 */
template <typename TInputImage, typename TSupportWindowImage, typename TOutputImage>
class ITK_TEMPLATE_EXPORT Spectra1DImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(Spectra1DImageFilter);

  itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension);

  using InputImageType = TInputImage;
  using SupportWindowImageType = TSupportWindowImage;
  using OutputImageType = TOutputImage;

  using ScalarType = typename DefaultConvertPixelTraits<typename OutputImageType::PixelType>::ComponentType;

  /** Standard class type alias. */
  using Self = Spectra1DImageFilter;
  using Superclass = ImageToImageFilter<InputImageType, OutputImageType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  itkOverrideGetNameOfClassMacro(Spectra1DImageFilter);
  itkNewMacro(Self);

  /** Set/get the input image describing support windows for local spectra
   * computation. This may be generated with itk::Spectra1DSupportWindowImageFilter. */
  itkSetInputMacro(SupportWindowImage, SupportWindowImageType);
  itkGetInputMacro(SupportWindowImage, SupportWindowImageType);

  /** Set/get an optional reference spectra image use to normalize the
   * output, such as from a phantom image.*/
  itkSetInputMacro(ReferenceSpectraImage, OutputImageType);
  itkGetInputMacro(ReferenceSpectraImage, OutputImageType);

protected:
  Spectra1DImageFilter();
  ~Spectra1DImageFilter() override = default;

  using OutputImageRegionType = typename OutputImageType::RegionType;

  void
  GenerateOutputInformation() override;
  void
  BeforeThreadedGenerateData() override;
  void
  ThreadedGenerateData(const OutputImageRegionType & outputRegionForThread, ThreadIdType threadId) override;
  void
  VerifyInputInformation() const override;

private:
  using ComplexType = std::complex<ScalarType>;
  using ComplexVectorType = vnl_vector<ComplexType>;
  using SpectraVectorType = std::vector<ScalarType>;
  using IndexType = typename InputImageType::IndexType;
  using SpectraLineType = std::pair<IndexType, SpectraVectorType>;
  using SpectraLinesContainerType = std::list<SpectraLineType>;
  using SupportWindowType = typename SupportWindowImageType::PixelType;
  using InputImageIteratorType = ImageRegionConstIterator<InputImageType>;
  using FFT1DType = vnl_fft_1d<ScalarType>;

  using Spectra1DSupportWindowFilterType = Spectra1DSupportWindowImageFilter<InputImageType>;
  using FFT1DSizeType = typename Spectra1DSupportWindowFilterType::FFT1DSizeType;

  using LineWindowMapType = std::unordered_map<FFT1DSizeType, SpectraVectorType>;

  struct PerThreadData
  {
    ComplexVectorType                 ComplexVector;
    SpectraVectorType                 SpectraVector;
    typename InputImageType::SizeType LineImageRegionSize;
    LineWindowMapType                 LineWindowMap;
  };
  using PerThreadDataContainerType = std::vector<PerThreadData>;
  PerThreadDataContainerType m_PerThreadDataContainer;

  typename OutputImageType::Pointer m_ReferenceSpectraImage;

  void
  ComputeSpectra(const IndexType & lineIndex, ThreadIdType threadId, SpectraLineType & spectraLine);
  static void
  AddLineWindow(FFT1DSizeType length, LineWindowMapType & lineWindowMap);
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkSpectra1DImageFilter.hxx"
#endif

#endif // itkSpectra1DImageFilter_h
