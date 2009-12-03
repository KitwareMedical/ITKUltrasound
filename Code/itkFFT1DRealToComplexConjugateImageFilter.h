#ifndef __itkFFT1DRealToComplexConjugateImageFilter_h
#define __itkFFT1DRealToComplexConjugateImageFilter_h

#include <complex>

#include "itkImage.h"
#include "itkImageToImageFilter.h"

namespace itk
{
/** @class itkFFT1DRealToComplexConjugateImageFilter
 * @brief Perform the Fast Fourier Transform, in the forward direction, with
 * real inputs, but only along one dimension.
 *
 * @ingroup FourierTransform
 */
template <class TPixel, unsigned int Dimension = 3>
class ITK_EXPORT FFT1DRealToComplexConjugateImageFilter:
  public ImageToImageFilter< TPixel, Dimension >
{
public:
  /** Standard class typedefs. */ 
  typedef Image<TPixel,VDimension>                    TInputImageType;
  typedef Image< std::complex< TPixel > , VDimension> TOutputImageType;

  typedef FFT1DRealToComplexConjugateImageFilter Self;
  typedef ImageToImageFilter< TInputImageType, TOutputImageType > Superclass;
  typedef SmartPointer<Self>                                      Pointer;
  typedef SmartPointer<const Self>                                ConstPointer;

  itkTypeMacro( FFT1DRealToComplexConjugateImageFilter, ImageToImageFilter );

  /** Customized object creation methods that support configuration-based 
    * selection of FFT implementation.
    *
    * Default implementation is VnlFFT1D.
    */
  static Pointer New(void);

protected:
  FFT1DRealToComplexConjugateImageFilter() {}
  virtual ~FFT1DRealToComplexConjugateImageFilter() {}

  virtual void GenerateInputRequestedRegion();
  virtual void GenerateInputRequestedRegion(); 
  virtual void EnlargeOutputRequestedRegion(DataObject *output); 
  virtual bool FullMatrix() = 0; // must be implemented in child

private:
  FFT1DRealToComplexConjugateImageFilter( const Self& );
  void operator=( const Self& );
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#ifndef __itkVnlFFT1DRealToComplexConjugateImageFilter_h
#ifndef __itkVnlFFT1DRealToComplexConjugateImageFilter_txx
#ifndef __itkFFTW1DRealToComplexConjugateImageFilter_h
#ifndef __itkFFTW1DRealToComplexConjugateImageFilter_txx
#include "itkFFT1DRealToComplexConjugateImageFilter.txx"
#endif
#endif
#endif
#endif
#endif

#endif // __itkFFT1DRealToComplexConjugateImageFilter_h
