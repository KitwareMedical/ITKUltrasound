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
#ifndef itkBlockMatchingNormalizedCrossCorrelationFFTMetricImageFilter_h
#define itkBlockMatchingNormalizedCrossCorrelationFFTMetricImageFilter_h

#include "itkBlockMatchingNormalizedCrossCorrelationMetricImageFilter.h"

#include "itkComplexConjugateImageFilter.h"
#include "itkHalfHermitianToRealInverseFFTImageFilter.h"
#include "itkConstantPadImageFilter.h"
#include "itkRealToHalfHermitianForwardFFTImageFilter.h"
#include "itkFFTShiftImageFilter.h"
#include "itkRegionFromReferenceImageFilter.h"

#include "itkMultiplyImageFilter.h"

namespace itk
{
namespace BlockMatching
{

/** \class NormalizedCrossCorrelationFFTMetricImageFilter
 *
 * \brief Create an image of the the normalized cross correlation with a kernel
 * calculated with FFT based correlation (multiply forward FFT's and take IFFT).
 *
 * \sa NormalizedCrossCorrelationMetricImageFilter
 *
 * \ingroup Ultrasound
 */
template <class TFixedImage, class TMovingImage, class TMetricImage>
class ITK_TEMPLATE_EXPORT NormalizedCrossCorrelationFFTMetricImageFilter
  : public NormalizedCrossCorrelationMetricImageFilter<TFixedImage, TMovingImage, TMetricImage>
{
public:
  /** Standard class type alias. */
  using Self = NormalizedCrossCorrelationFFTMetricImageFilter;
  using Superclass = NormalizedCrossCorrelationMetricImageFilter<TFixedImage, TMovingImage, TMetricImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(NormalizedCrossCorrelationFFTMetricImageFilter);

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int, TFixedImage::ImageDimension);

  /** Type of the fixed image. */
  using FixedImageType = typename Superclass::FixedImageType;
  using FixedImageConstPointerType = typename FixedImageType::ConstPointer;

  /** Type of the moving image. */
  using MovingImageType = typename Superclass::MovingImageType;
  using MovingImageConstPointerType = typename MovingImageType::ConstPointer;

  /** Type of the metric image. */
  using MetricImageType = typename Superclass::MetricImageType;
  using MetricImagePointerType = typename MetricImageType::Pointer;
  using MetricImageConstPointerType = typename MetricImageType::ConstPointer;
  using MetricImageRegionType = typename MetricImageType::RegionType;
  using MetricImagePixelType = typename MetricImageType::PixelType;

  /**
   * Set/Get the greatest prime factor allowed on the size of the padded image.
   * The filter increase the size of the image to reach a size with the greatest
   * prime factor smaller or equal to the specified value. The default value is
   * 13, which is the greatest prime number for which the FFT are precomputed
   * in FFTW, and thus gives very good performance.
   * A greatest prime factor of 2 produce a size which is a power of 2, and thus
   * is suitable for vnl base fft filters.
   * A greatest prime factor of 1 or less - typically 0 - disable the extra padding.
e  */
  itkGetConstMacro(SizeGreatestPrimeFactor, SizeValueType);
  itkSetMacro(SizeGreatestPrimeFactor, SizeValueType);

protected:
  NormalizedCrossCorrelationFFTMetricImageFilter();

  virtual void
  GenerateData() override;

  using PadFilterType = ConstantPadImageFilter<MetricImageType, MetricImageType>;
  using FFTShiftFilterType = FFTShiftImageFilter<MetricImageType, MetricImageType>;
  using FFTFilterType = RealToHalfHermitianForwardFFTImageFilter<MetricImageType>;
  using ComplexImageType = typename FFTFilterType::OutputImageType;
  using IFFTFilterType = HalfHermitianToRealInverseFFTImageFilter<ComplexImageType>;
  using ComplexConjugateFilterType = ComplexConjugateImageFilter<ComplexImageType, ComplexImageType>;
  using MultiplyFilterType = MultiplyImageFilter<ComplexImageType, ComplexImageType, ComplexImageType>;
  using CropFilterType = RegionFromReferenceImageFilter<MetricImageType, MetricImageType>;

  typename PadFilterType::Pointer              m_KernelPadFilter;
  typename PadFilterType::Pointer              m_MovingPadFilter;
  typename FFTShiftFilterType::Pointer         m_FFTShiftFilter;
  typename FFTFilterType::Pointer              m_KernelFFTFilter;
  typename FFTFilterType::Pointer              m_MovingFFTFilter;
  typename ComplexConjugateFilterType::Pointer m_ComplexConjugateImageFilter;
  typename MultiplyFilterType::Pointer         m_MultiplyFilter;
  typename IFFTFilterType::Pointer             m_IFFTFilter;
  typename CropFilterType::Pointer             m_CropFilter;

  SizeValueType m_SizeGreatestPrimeFactor;

private:
  NormalizedCrossCorrelationFFTMetricImageFilter(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented
};

} // end namespace BlockMatching
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkBlockMatchingNormalizedCrossCorrelationFFTMetricImageFilter.hxx"
#endif

#endif
