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
#ifndef itkAttenuationImageFilter_h
#define itkAttenuationImageFilter_h

#include <vector>

#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include "itkImageRegionSplitterDirection.h"
#include "itkMacro.h"

#include "itkNumericTraits.h"

namespace itk
{
/** \class AttenuationImageFilter
 * \brief Computes the estimated attentuation in dB/(MHz*cm)
 *
 * Attenuation is a measure of how an RF signal fades in strength
 * as it passes through a physical region. In ultrasound analysis
 * signal attenuation tends to be roughly similar over areas
 * of similar material composition, such as different types of
 * tissue within an image.
 *
 * AttenuationImageFilter receives an input vector image representing
 * RF spectra. One image direction represents the direction of an
 * RF waveform emitted from an ultrasound probe. Remaining image
 * directions may represent directions in physical space, such as
 * where the probe contains multiple elements, or directions in
 * time as multiple samples are captured in a sweep with the probe.
 * Each pixel in the input image is a vector representing frequency
 * components at bins based on the sampling frequency.
 * The filter also receives a mandatory mask input indicating the
 * region over which attenuations should be estimated.
 *
 * AttenuationImageFilter generates a scalar output image with
 * pixel intensities representing attenuation estimates between
 * the given pixel's location in the input image and either a pixel
 * that is a fixed distance away or a pixel at the end of the given
 * mask region, depending on filter settings. Pixels outside of
 * the mask after padding erosion will have a value of zero.
 *
 * \sa MaskedImageToHistogramFilter
 * \sa Spectra1DImageFilter
 * \sa VariableLengthVector
 * \sa VectorImage
 *
 * \ingroup ITKImageStatistics
 * \ingroup Ultrasound
 */
template <typename TInputImage,
          typename TOutputImage = itk::Image<float, TInputImage::ImageDimension>,
          typename TMaskImage = itk::Image<unsigned char, TInputImage::ImageDimension>>
class ITK_TEMPLATE_EXPORT AttenuationImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(AttenuationImageFilter);

  /** Standard class type aliases. */
  using Self = AttenuationImageFilter;
  using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(AttenuationImageFilter, ImageToImageFilter);

  /** Image type alias support */
  static constexpr unsigned int ImageDimension = TInputImage::ImageDimension;

  using InputImageType = TInputImage;
  using InputImagePointer = typename TInputImage::Pointer;

  using OutputImageType = TOutputImage;
  using OutputImagePointer = typename TOutputImage::Pointer;

  using MaskImageType = TMaskImage;
  using MaskImagePointer = typename TMaskImage::Pointer;
  using MaskPixelType = typename TMaskImage::PixelType;

  using InputRegionType = typename TInputImage::RegionType;
  using InputSizeType = typename TInputImage::SizeType;
  using InputIndexType = typename TInputImage::IndexType;
  using InputPixelType = typename TInputImage::PixelType;

  using OutputRegionType = typename TOutputImage::RegionType;
  using OutputPixelType = typename TOutputImage::PixelType;

  /** Input mask image represents input region for analysis */
  itkSetInputMacro(InputMaskImage, TMaskImage);
  itkGetInputMacro(InputMaskImage, TMaskImage);

  /** Output mask image represents output region for analysis
   *  after padding is applied */
  itkGetConstMacro(OutputMaskImage, TMaskImage *);

  /** Label value indicating which mask pixel values should be included in analysis.
   *  A value of zero indicates that any nonzero pixel should be included.
   */
  itkSetMacro(LabelValue, MaskPixelType);
  itkGetConstMacro(LabelValue, MaskPixelType);

  /** Fix the pixel distance between voxels for estimating
   *  attenuation in a scan line.
   *  If set to zero then the last continuous pixel
   *  in the inclusion will always be chosen as the
   *  second pixel for attenuation calculation. */
  itkSetMacro(FixedEstimationDepth, unsigned int);
  itkGetConstMacro(FixedEstimationDepth, unsigned int);

  /** Set/get fixed estimation depth in physical space.
   *  Assumes input image spacing is in millimenters. */
  void
  SetFixedEstimationDepthMM(const float distanceMM);
  float
  GetFixedEstimationDepthMM() const;

  /** RF scanline direction */
  itkSetMacro(Direction, unsigned int);
  itkGetConstMacro(Direction, unsigned int);

  /** RF sampling frequency */
  itkSetMacro(SamplingFrequencyMHz, float);
  itkGetConstMacro(SamplingFrequencyMHz, float);

  /** Low end of RF frequency band. Must be a positive value. */
  itkSetMacro(FrequencyBandStartMHz, float);
  itkGetConstMacro(FrequencyBandStartMHz, float);

  /* High end of RF frequency band. Must be a positive value.*/
  itkSetMacro(FrequencyBandEndMHz, float);
  itkGetConstMacro(FrequencyBandEndMHz, float);

  /** Optionally discard negative attenuation estimates
   *  so that they are not considered in statistic computations.
   *  Negative attenuation implies that a signal strengthened
   *  while passing through tissue and may result from
   *  sampling error or external interference. */
  itkSetMacro(ConsiderNegativeAttenuations, bool);
  itkGetConstMacro(ConsiderNegativeAttenuations, bool);

  /** Skip attenuation estimation for a fixed number
   *  of pixels at the start of an inclusion region.
   *  Applied before spatial padding.
   *  Can help with uncertainty at borders of mask.
   */
  itkSetMacro(PadUpperBounds, unsigned int);
  itkGetConstMacro(PadUpperBounds, unsigned int);
  /** Skip attenuation estimation at the end of an inclusion region. */
  itkSetMacro(PadLowerBounds, unsigned int);
  itkGetConstMacro(PadLowerBounds, unsigned int);

  /** Set padding at start of inclusion region based on
   *  physical distance. Assumes input image is in millimeters. */
  void
  SetPadLowerBoundsMM(const float distanceMM);
  float
  GetPadLowerBoundsMM() const;

  /** Set padding at end of inclusion region based on
   *  physical distance. Assumes input image is in millimeters. */
  void
  SetPadUpperBoundsMM(const float distanceMM);
  float
  GetPadUpperBoundsMM() const;

  // Alias for setting direction of RF waveform in data collection
  void
  SetScanDirection(unsigned int direction)
  {
    this->SetDirection(direction);
  };
  unsigned int
  GetScanDirection()
  {
    return this->GetDirection();
  };

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

protected:
  AttenuationImageFilter();
  ~AttenuationImageFilter() override = default;


  void
  VerifyPreconditions() const override;

  /** Allocate buffers and other initializations before threaded execution. */
  void
  BeforeThreadedGenerateData() override;

  void
  ThreadedGenerateData(const OutputRegionType & regionForThread, ThreadIdType) override;

  const ImageRegionSplitterBase *
  GetImageRegionSplitter() const override;

  /** Compute attenuation between two pixels in the RF spectra vector image.
   *  Assumes that image spacing is in MM. */
  OutputPixelType
  ComputeAttenuation(const InputIndexType & end, const InputIndexType & start) const;

  /** Transform spatial distance along an RF scan line to continuous pixel distance
   *  in the input image.
   *  Assumes input distance and image spacing are in millimeters.
   *  Rounded to nearest integer. */
  float
  TransformPhysicalToPixelScanLineDistance(float distanceMM) const;

  /** Transform pixel distance along an RF scan line to spatial distance
   *  in the input image.
   *  Assumes input image spacing is in millimeters. Output is in millimeters. */
  float
  TransformPixelToPhysicalScanLineDistance(unsigned int distance) const;

  /** Check whether given pixel index is included in the mask */
  bool
  ThreadedIsIncluded(InputIndexType index) const;

private:
  MaskPixelType m_LabelValue = 0;

  unsigned int m_Direction = 0;

  float m_ScanStepMM = 1.0f;

  unsigned int m_FixedEstimationDepth = 0;

  std::vector<float> m_DistanceWeights;

  float m_SamplingFrequencyMHz = 0.0f;

  float m_FrequencyBandStartMHz = 0.0f;
  float m_FrequencyBandEndMHz = 0.0f;
  float m_FrequencyDelta = 0.0f;

  // Frequency band to consider for attenuation
  unsigned int m_StartComponent = 0;
  unsigned int m_EndComponent = 0;
  unsigned int m_ConsideredComponents = 1;

  bool m_ConsiderNegativeAttenuations = false;

  unsigned int m_PadUpperBounds = 0;
  unsigned int m_PadLowerBounds = 0;

  /** Region splitter to ensure scanline is intact in threaded regions */
  ImageRegionSplitterDirection::Pointer m_RegionSplitter = ImageRegionSplitterDirection::New();

  /** Cache mask image reference before threaded execution to reduce calls to GetMaskImage() */
  mutable const MaskImageType * m_ThreadedInputMaskImage;

  /** Output mask image may be eroded via m_PadUpperBounds and m_PadLowerBounds
   *  along scan line direction */
  MaskImagePointer m_OutputMaskImage = MaskImageType::New();
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkAttenuationImageFilter.hxx"
#endif

#endif
