#ifndef __itkBlockMatchingNormalizedCrossCorrelationFFTMetricImageFilter_h
#define __itkBlockMatchingNormalizedCrossCorrelationFFTMetricImageFilter_h

#include "itkBlockMatchingNormalizedCrossCorrelationMetricImageFilter.h"

#include "itkComplexConjugateImageFilter.h"
#include "itkHalfHermitianToRealInverseFFTImageFilter.h"
#include "itkFFTPadImageFilter.h"
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
template< class TFixedImage, class TMovingImage, class TMetricImage >
class ITK_EXPORT NormalizedCrossCorrelationFFTMetricImageFilter :
  public NormalizedCrossCorrelationMetricImageFilter< TFixedImage,
                                                      TMovingImage, TMetricImage >
{
public:
  /** Standard class typedefs. */
  typedef NormalizedCrossCorrelationFFTMetricImageFilter Self;
  typedef NormalizedCrossCorrelationMetricImageFilter< TFixedImage,
                                                       TMovingImage, TMetricImage > Superclass;
  typedef SmartPointer<Self>       Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(NormalizedCrossCorrelationFFTMetricImageFilter,
               NormalizedCrossCorrelationMetricImageFilter);

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TFixedImage::ImageDimension);

  /** Type of the fixed image. */
  typedef typename Superclass::FixedImageType   FixedImageType;
  typedef typename FixedImageType::ConstPointer FixedImageConstPointerType;

  /** Type of the moving image. */
  typedef typename Superclass::MovingImageType   MovingImageType;
  typedef typename MovingImageType::ConstPointer MovingImageConstPointerType;

  /** Type of the metric image. */
  typedef typename Superclass::MetricImageType   MetricImageType;
  typedef typename MetricImageType::Pointer      MetricImagePointerType;
  typedef typename MetricImageType::ConstPointer MetricImageConstPointerType;
  typedef typename MetricImageType::RegionType   MetricImageRegionType;
  typedef typename MetricImageType::PixelType    MetricImagePixelType;

  /**
   * Set/Get the greatest prime factor allowed on the size of the padded image.
   * The filter increase the size of the image to reach a size with the greatest
   * prime factor smaller or equal to the specified value. The default value is
   * 13, which is the greatest prime number for which the FFT are precomputed
   * in FFTW, and thus gives very good performance.
   * A greatest prime factor of 2 produce a size which is a power of 2, and thus
   * is suitable for vnl base fft filters.
   * A greatest prime factor of 1 or less - typically 0 - disable the extra padding.
   *
   * Warning: this parameter is not used (and useful) only when ITK is built with
   * FFTW support.
   */
  itkGetConstMacro(GreatestPrimeFactor, int);
  itkSetMacro(GreatestPrimeFactor, int);

protected:
  NormalizedCrossCorrelationFFTMetricImageFilter();

  virtual void GenerateData();

  typedef itk::FFTPadImageFilter< MetricImageType, MetricImageType >
  PadFilterType;

  typedef itk::FFTShiftImageFilter< MetricImageType, MetricImageType >
  FFTShiftFilterType;

  typedef itk::RealToHalfHermitianForwardFFTImageFilter< MetricImageType > FFTFilterType;
  typedef typename FFTFilterType::OutputImageType ComplexImageType;
  typedef itk::HalfHermitianToRealInverseFFTImageFilter< ComplexImageType > IFFTFilterType;

  typedef itk::ComplexConjugateImageFilter< ComplexImageType,
                                            ComplexImageType > ComplexConjugateFilterType;

  typedef itk::MultiplyImageFilter< ComplexImageType, ComplexImageType,
                                    ComplexImageType > MultiplyFilterType;

  typedef itk::RegionFromReferenceImageFilter< MetricImageType,
                                               MetricImageType > CropFilterType;

  typename PadFilterType::Pointer m_PadFilter;

  typename FFTShiftFilterType::Pointer  m_FFTShiftFilter;

  typename FFTFilterType::Pointer m_KernelFFTFilter;
  typename FFTFilterType::Pointer m_MovingFFTFilter;

  typename ComplexConjugateFilterType::Pointer m_ComplexConjugateImageFilter;

  typename MultiplyFilterType::Pointer  m_MultiplyFilter;

  typename IFFTFilterType::Pointer m_IFFTFilter;

  typename CropFilterType::Pointer m_CropFilter;

  int m_GreatestPrimeFactor;

private:
  NormalizedCrossCorrelationFFTMetricImageFilter( const Self& ); // purposely not implemented
  void operator=( const Self& ); // purposely not implemented
};

} // end namespace BlockMatching
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlockMatchingNormalizedCrossCorrelationFFTMetricImageFilter.hxx"
#endif

#endif
