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
#ifndef itkBackscatterImageFilter_h
#define itkBackscatterImageFilter_h

#include <vector>

#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include "itkMacro.h"
#include "itkNumericTraits.h"

namespace itk
{
/** \class BackscatterImageFilter
 * \brief Computes the estimated backscatter coefficient
 *
 * BackscatterImageFilter receives an input vector image representing
 * RF spectra. One image direction represents the direction of an
 * RF waveform emitted from an ultrasound probe. Remaining image
 * directions may represent directions in physical space, such as
 * where the probe contains multiple elements, or directions in
 * time as multiple samples are captured in a sweep with the probe.
 * Each pixel in the input image is a vector representing frequency
 * components at bins based on the sampling frequency.
 *
 * BackscatterImageFilter generates three scalar metric output images with
 * intensities representing backscatter estimates.
 *
 * Output 0 (primary): Average of intensities of selected frequency components.
 * Output 1: Negative slope of a line fit to the selected frequency components.
 *   The slope is negated because it is expected to be mostly negative.
 * Output 2: Intercept of a line fit to the selected frequency components.
 *
 * The slope and intercept correspond to the classic Lizzi-Feleppa parameters.
 *
 * Lizzi, Frederic L., S. Kaisar Alam, Samuel Mikaelian, Paul Lee, and Ernest J. Feleppa.
 * "On the statistics of ultrasonic spectral parameters."
 * Ultrasound in medicine & biology 32, no. 11 (2006): 1671-1685.
 *
 * Lizzi, Frederic L., Michael Greenebaum, Ernest J. Feleppa, Marek Elbaum, and D. Jackson Coleman.
 * "Theoretical framework for spectrum analysis in ultrasonic tissue characterization."
 * The Journal of the Acoustical Society of America 73, no. 4 (1983): 1366-1373.
 *
 * They usually used the mid-band fit instead of the average.
 *
 * \sa MaskedImageToHistogramFilter
 * \sa Spectra1DImageFilter
 * \sa VariableLengthVector
 * \sa VectorImage
 *
 * \ingroup ITKImageStatistics
 * \ingroup Ultrasound
 */
template <typename TInputImage, typename TOutputImage = Image<float, TInputImage::ImageDimension>>
class ITK_TEMPLATE_EXPORT BackscatterImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(BackscatterImageFilter);

  /** Standard class type aliases. */
  using Self = BackscatterImageFilter;
  using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(BackscatterImageFilter, ImageToImageFilter);

  /** Image type alias support */
  static constexpr unsigned int ImageDimension = TInputImage::ImageDimension;

  using InputImageType = TInputImage;
  using InputImagePointer = typename TInputImage::Pointer;

  using OutputImageType = TOutputImage;
  using OutputImagePointer = typename TOutputImage::Pointer;

  using InputRegionType = typename TInputImage::RegionType;
  using InputSizeType = typename TInputImage::SizeType;
  using InputIndexType = typename TInputImage::IndexType;
  using InputPixelType = typename TInputImage::PixelType;

  using OutputRegionType = typename TOutputImage::RegionType;
  using OutputPixelType = typename TOutputImage::PixelType;

  /** RF sampling frequency */
  itkSetMacro(SamplingFrequencyMHz, float);
  itkGetConstMacro(SamplingFrequencyMHz, float);

  /** Low end of RF frequency band to use in backscatter analysis.
   * Must be a positive value. */
  itkSetMacro(FrequencyBandStartMHz, float);
  itkGetConstMacro(FrequencyBandStartMHz, float);

  /* High end of RF frequency band to use in backscatter analysis.
   * Must be a positive value.*/
  itkSetMacro(FrequencyBandEndMHz, float);
  itkGetConstMacro(FrequencyBandEndMHz, float);

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

protected:
  BackscatterImageFilter();
  ~BackscatterImageFilter() override = default;

  void
  VerifyPreconditions() const override;

  /** Allocate buffers and other initializations before threaded execution. */
  void
  BeforeThreadedGenerateData() override;

  void
  DynamicThreadedGenerateData(const OutputRegionType & regionForThread) override;

  /** Compute backscatter for a pixel in the RF spectra vector image. */
  void
  ComputeBackscatter(const InputIndexType & index,
                     OutputPixelType &      average,
                     OutputPixelType &      slope,
                     OutputPixelType &      intercept) const;

private:
  float m_SamplingFrequencyMHz = 0.0f;
  float m_FrequencyBandStartMHz = 0.0f;
  float m_FrequencyBandEndMHz = 0.0f;
  float m_FrequencyDelta = 0.0f;

  // Frequency band to consider for backscatter
  unsigned int m_StartComponent = 0;
  unsigned int m_EndComponent = 0;
  unsigned int m_ConsideredComponents = 1;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkBackscatterImageFilter.hxx"
#endif

#endif
