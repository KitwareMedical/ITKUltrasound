#ifndef __itkFFT1DRealToComplexConjugateImageFilter_h
#define __itkFFT1DRealToComplexConjugateImageFilter_h

#include <complex>

#include "itkImageToImageFilter.h"
#include "itkImageRegionSplitterDirection.h"

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
  typedef Image<TPixel,VDimension>                    InputImageType;
  typedef Image< std::complex< TPixel > , VDimension> OutputImageType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  typedef FFT1DRealToComplexConjugateImageFilter Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType > Superclass;
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
  FFT1DRealToComplexConjugateImageFilter();
  virtual ~FFT1DRealToComplexConjugateImageFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  virtual void GenerateInputRequestedRegion();
  virtual void EnlargeOutputRequestedRegion(DataObject *output);

  /** Direction in which the filter is to be applied
   * this should be in the range [0,ImageDimension-1]. */
  unsigned int m_Direction;

  virtual void BeforeThreadedGenerateData();

  /** Override to return a splitter that does not split along the direction we
   * are performing the transform. */
  virtual const ImageRegionSplitterBase* GetImageRegionSplitter(void) const;

private:
  FFT1DRealToComplexConjugateImageFilter( const Self& );
  void operator=( const Self& );

  ImageRegionSplitterDirection::Pointer m_ImageRegionSplitter;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#ifndef __itkVnlFFT1DRealToComplexConjugateImageFilter_h
#ifndef __itkVnlFFT1DRealToComplexConjugateImageFilter_hxx
#ifndef __itkFFTW1DRealToComplexConjugateImageFilter_h
#ifndef __itkFFTW1DRealToComplexConjugateImageFilter_hxx
#include "itkFFT1DRealToComplexConjugateImageFilter.hxx"
#endif
#endif
#endif
#endif
#endif

#endif // __itkFFT1DRealToComplexConjugateImageFilter_h
