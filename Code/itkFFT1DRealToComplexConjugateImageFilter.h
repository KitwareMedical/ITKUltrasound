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
template <class TPixel, unsigned int VDimension = 3>
class ITK_EXPORT FFT1DRealToComplexConjugateImageFilter:
  public ImageToImageFilter< Image< TPixel, VDimension >,
			     Image< std::complex< TPixel >, VDimension > >
{
public:
  /** Standard class typedefs. */ 
  typedef Image<TPixel,VDimension>                    TInputImageType;
  typedef Image< std::complex< TPixel > , VDimension> TOutputImageType;
  typedef typename TOutputImageType::RegionType OutputImageRegionType;

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

  /** Get the direction in which the filter is to be applied. */
  itkGetMacro(Direction, unsigned int);

  /** Set the direction in which the filter is to be applied. */
  itkSetMacro(Direction, unsigned int);

protected:
  FFT1DRealToComplexConjugateImageFilter(): m_Direction(0) {}
  virtual ~FFT1DRealToComplexConjugateImageFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  virtual void GenerateInputRequestedRegion(); 
  virtual void EnlargeOutputRequestedRegion(DataObject *output); 

  /** Split the output's RequestedRegion into "num" pieces, returning
   * region "i" as "splitRegion". This method is called "num" times. The
   * regions must not overlap. The method returns the number of pieces that
   * the routine is capable of splitting the output RequestedRegion,
   * i.e. return value is less than or equal to "num". 
   *
   * This method is reimplemented from itk::ImageSource to ensure that the
   * region does not get split along the direction of the transform. */
  virtual
  int SplitRequestedRegion(int i, int num, OutputImageRegionType& splitRegion);

  /** Direction in which the filter is to be applied
   * this should be in the range [0,ImageDimension-1]. */
  unsigned int m_Direction;

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
#ifndef __itkOpenCL1DRealToComplexConjugateImageFilter_h
#ifndef __itkOpenCL1DRealToComplexConjugateImageFilter_txx
#include "itkFFT1DRealToComplexConjugateImageFilter.txx"
#endif
#endif
#endif
#endif
#endif
#endif
#endif

#endif // __itkFFT1DRealToComplexConjugateImageFilter_h
