#ifndef __itkAnalyticSignalImageFilter_h
#define __itkAnalyticSignalImageFilter_h

#include <complex>

#include "itkFFT1DComplexToComplexImageFilter.h"
#include "itkFFT1DRealToComplexConjugateImageFilter.h"
#include "itkImage.h"
#include "itkImageToImageFilter.h"

namespace itk
{
/** @class itkAnalyticSignalImageFilter
 * @brief Generates the analytic signal from one direction of an image.
 * 
 * This filter generates the complex valued analytic signal along one direction
 * of an image.  This input is a real valued image, and the output is a complex
 * image.
 *
 * The analytic signal is given by
 *
 * f_a(x) = f(x) - i f_H(x)
 *
 * Where i is the square root of one and f_H(x) is the Hibert transform of f(x).
 *
 * Since the Hilbert transform in the Fourier domain is
 *
 * F_H(k) = F(k) i sign(k),
 *
 * f_a(x) is calculated by
 *
 * f_a(x) = F^{-1}( F(k) 2 U(k) )
 *
 * where U(k) is the unit step function.
 *
 * @ingroup FourierTransform
 */
template <class TPixel, unsigned int VDimension = 3>
class ITK_EXPORT AnalyticSignalImageFilter:
  public ImageToImageFilter< Image< TPixel, VDimension >,
			     Image< std::complex< TPixel >, VDimension > >
{
public:
  /** Standard class typedefs. */ 
  typedef Image<TPixel,VDimension>                    TInputImageType;
  typedef Image< std::complex< TPixel > , VDimension> TOutputImageType;
  typedef typename TOutputImageType::RegionType OutputImageRegionType;

  typedef AnalyticSignalImageFilter Self;
  typedef ImageToImageFilter< TInputImageType, TOutputImageType > Superclass;
  typedef SmartPointer<Self>                                      Pointer;
  typedef SmartPointer<const Self>                                ConstPointer;

  itkTypeMacro( AnalyticSignalImageFilter, ImageToImageFilter );
  itkNewMacro( Self );

  /** Get the direction in which the filter is to be applied. */
  virtual unsigned int GetDirection() const
    {
    return this->m_FFTRealToComplexFilter->GetDirection();
    }

  /** Set the direction in which the filter is to be applied. */
  virtual void SetDirection( const unsigned int direction )
    {
    if( this->m_FFTRealToComplexFilter->GetDirection() != direction )
      {
      this->m_FFTRealToComplexFilter->SetDirection( direction );
      this->m_FFTComplexToComplexFilter->SetDirection( direction );
      this->Modified();
      }
    }

protected:
  AnalyticSignalImageFilter();
  virtual ~AnalyticSignalImageFilter() {}
  void PrintSelf(std::ostream& os, Indent indent) const;

  // These behave like their analogs in FFT1DRealToComplexConjugateImageFilter.
  virtual void GenerateInputRequestedRegion(); 
  virtual void EnlargeOutputRequestedRegion(DataObject *output); 

  /** Split the output's RequestedRegion into "num" pieces, returning
   * region "i" as "splitRegion". This method is called "num" times. The
   * regions must not overlap. The method returns the number of pieces that
   * the routine is capable of splitting the output RequestedRegion,
   * i.e. return value is less than or equal to "num". 
   *
   * This method is reimplemented from itk::ImageSource to ensure that the
   * region does not get split along the direction of the filter. */
  virtual
  int SplitRequestedRegion(int i, int num, OutputImageRegionType& splitRegion);

  virtual void BeforeThreadedGenerateData ();
  virtual void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread, int threadId );
  virtual void AfterThreadedGenerateData ();

  typedef FFT1DRealToComplexConjugateImageFilter< TPixel, VDimension > FFTRealToComplexType;
  typedef FFT1DComplexToComplexImageFilter< TPixel, VDimension > FFTComplexToComplexType;

  typename FFTRealToComplexType::Pointer m_FFTRealToComplexFilter;
  typename FFTComplexToComplexType::Pointer m_FFTComplexToComplexFilter;

private:
  AnalyticSignalImageFilter( const Self& );
  void operator=( const Self& );
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAnalyticSignalImageFilter.txx"
#endif

#endif // __itkAnalyticSignalImageFilter_h
