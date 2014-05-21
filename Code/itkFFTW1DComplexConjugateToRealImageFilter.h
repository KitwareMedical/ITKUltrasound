#ifndef __itkFFTW1DComplexConjugateToRealImageFilter_h
#define __itkFFTW1DComplexConjugateToRealImageFilter_h

#include "itkFFT1DComplexConjugateToRealImageFilter.h"
#include "itkFFTWCommonExtended.h"

#include <vector>

namespace itk
{
/** \class FFTW1DComplexConjugateToRealImageFilter
 * \brief only do FFT along one dimension using FFTW as a backend.
 *
 * \ingroup
 */

template <class TPixel, unsigned int Dimension = 3>
class ITK_EXPORT FFTW1DComplexConjugateToRealImageFilter :
    public FFT1DComplexConjugateToRealImageFilter<TPixel,Dimension>
{
public:
  typedef FFTW1DComplexConjugateToRealImageFilter Self;
  typedef FFT1DComplexConjugateToRealImageFilter<TPixel,Dimension> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Standard class typedefs.*/
  typedef typename Superclass::InputImageType InputImageType;
  typedef typename Superclass::OutputImageType OutputImageType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  /**
   * the proxy type is a wrapper for the fftw API
   * since the proxy is only defined over double and float,
   * trying to use any other pixel type is inoperative, as
   * is trying to use double if only the float FFTW1D version is
   * configured in, or float if only double is configured.
   */
  typedef typename fftw::ComplexToComplexProxy<TPixel> FFTW1DProxyType;
  typedef typename std::vector< typename FFTW1DProxyType::PlanType > PlanArrayType;
  typedef typename std::vector< typename FFTW1DProxyType::ComplexType* > PlanBufferPointerType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(FFTW1DComplexConjugateToRealImageFilter,
               FFT1DComplexConjugateToRealImageFilter);


protected:
  FFTW1DComplexConjugateToRealImageFilter(): m_PlanComputed( false ),
    m_LastImageSize( 0 )
  {}
  virtual ~FFTW1DComplexConjugateToRealImageFilter()
  {
  if ( m_PlanComputed )
    {
    this->DestroyPlans();
    }
  }

  virtual void BeforeThreadedGenerateData();
  virtual void ThreadedGenerateData( const OutputImageRegionType&, ThreadIdType threadID );  // generates output from input

private:
  FFTW1DComplexConjugateToRealImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** Destroy FFTW Plans and associated buffers. */
  void DestroyPlans();

  bool m_PlanComputed;
  PlanArrayType m_PlanArray;
  unsigned int m_LastImageSize;
  PlanBufferPointerType m_InputBufferArray;
  PlanBufferPointerType m_OutputBufferArray;
};

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFFTW1DComplexConjugateToRealImageFilter.txx"
#endif

#endif //__itkFFTW1DComplexConjugateToRealImageFilter_h
