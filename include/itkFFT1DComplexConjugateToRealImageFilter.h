#ifndef __itkFFT1DComplexConjugateToRealImageFilter_h
#define __itkFFT1DComplexConjugateToRealImageFilter_h

#include <complex>

#include "itkImageToImageFilter.h"
#include "itkImageRegionSplitterDirection.h"

namespace itk
{
/** \class itkFFT1DComplexConjugateToRealImageFilter
 * \brief Perform the Fast Fourier Transform, in the reverse direction, with
 * real output, but only along one dimension.
 *
 * \ingroup FourierTransform
 * \ingroup FFT1D
 */
template <class TPixel, unsigned int VDimension = 3>
class ITK_EXPORT FFT1DComplexConjugateToRealImageFilter:
  public ImageToImageFilter< Image< std::complex< TPixel >, VDimension >,
			     Image< TPixel, VDimension > >
{
public:
  /** Standard class typedefs. */ 
  typedef Image< std::complex< TPixel > , VDimension> InputImageType;
  typedef Image<TPixel,VDimension>                    OutputImageType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  typedef FFT1DComplexConjugateToRealImageFilter		  Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType > Superclass;
  typedef SmartPointer<Self>                                      Pointer;
  typedef SmartPointer<const Self>                                ConstPointer;

  itkStaticConstMacro(ImageDimension, unsigned int,
                      InputImageType::ImageDimension );

  itkTypeMacro( FFT1DComplexConjugateToRealImageFilter, ImageToImageFilter );

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
  FFT1DComplexConjugateToRealImageFilter();
  virtual ~FFT1DComplexConjugateToRealImageFilter() {}
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
  FFT1DComplexConjugateToRealImageFilter( const Self& );
  void operator=( const Self& );

  ImageRegionSplitterDirection::Pointer m_ImageRegionSplitter;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#ifndef __itkVnlFFT1DComplexConjugateToRealImageFilter_h
#ifndef __itkVnlFFT1DComplexConjugateToRealImageFilter_hxx
#ifndef __itkFFTW1DComplexConjugateToRealImageFilter_h
#ifndef __itkFFTW1DComplexConjugateToRealImageFilter_hxx
#include "itkFFT1DComplexConjugateToRealImageFilter.hxx"
#endif
#endif
#endif
#endif
#endif

#endif // __itkFFT1DComplexConjugateToRealImageFilter_h
