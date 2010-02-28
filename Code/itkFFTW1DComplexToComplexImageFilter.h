#ifndef __itkFFTW1DComplexToComplexImageFilter_h
#define __itkFFTW1DComplexToComplexImageFilter_h

#include "itkFFT1DComplexToComplexImageFilter.h"
#include "itkFFTWCommonExtended.h"


namespace itk
{
/** /class FFTW1DComplexToComplexImageFilter
 * /brief only do FFT along one dimension using FFTW as a backend.
 *
 * \ingroup
 */

template <class TPixel, unsigned int Dimension = 3>
class ITK_EXPORT FFTW1DComplexToComplexImageFilter :
    public FFT1DComplexToComplexImageFilter<TPixel,Dimension>
{
public:
  typedef FFTW1DComplexToComplexImageFilter Self;
  typedef FFT1DComplexToComplexImageFilter<TPixel,Dimension> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Standard class typedefs.*/
  typedef typename Superclass::TInputImageType TInputImageType;
  typedef typename Superclass::TOutputImageType TOutputImageType;

  /**
   * the proxy type is a wrapper for the fftw API
   * since the proxy is only defined over double and float,
   * trying to use any other pixel type is inoperative, as
   * is trying to use double if only the float FFTW1D version is
   * configured in, or float if only double is configured.
   */
  typedef typename fftw::ComplexToComplexProxy<TPixel> FFTW1DProxyType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(FFTW1DComplexToComplexImageFilter,
               FFT1DComplexToComplexImageFilter);


protected:
  FFTW1DComplexToComplexImageFilter() : m_PlanComputed(false),
                                            m_LastImageSize(0),
                                            m_InputBuffer(0),
                                            m_OutputBuffer(0)
  {
  }
  virtual ~FFTW1DComplexToComplexImageFilter()
  {
    if(m_PlanComputed)
      {
      FFTW1DProxyType::DestroyPlan(this->m_Plan);
      delete [] this->m_InputBuffer;
      delete [] this->m_OutputBuffer;
      }
  }

  virtual void GenerateData();  // generates output from input

private:
  FFTW1DComplexToComplexImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  bool m_PlanComputed;
  typename FFTW1DProxyType::PlanType m_Plan;
  unsigned int m_LastImageSize;
  typename FFTW1DProxyType::ComplexType *m_InputBuffer;
  typename FFTW1DProxyType::ComplexType	*m_OutputBuffer;

};

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFFTW1DComplexToComplexImageFilter.txx"
#endif

#endif //__itkFFTW1DComplexToComplexImageFilter_h
