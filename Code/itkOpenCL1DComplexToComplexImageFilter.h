#ifndef __itkOpenCL1DComplexToComplexImageFilter_h
#define __itkOpenCL1DComplexToComplexImageFilter_h

#include "itkFFT1DComplexToComplexImageFilter.h"

#define __CL_ENABLE_EXCEPTIONS
#include "CL/cl.hpp"
#include "clFFT.h"

namespace itk
{
/** /class OpenCL1DComplexToComplexImageFilter
 * /brief only do FFT along one dimension using OpenCL_FFT as a backend.
 *
 * The size of the image in the transformed direction must be a power of 2.
 *
 * There is considerable overhead to generate the FFT plan, which occurs
 * whenever the input image size changes.  Therefore, the throughput benefit
 * will only be realized for large images or many small images of
 * the same size.
 *
 * \ingroup
 */

template <class TPixel, unsigned int Dimension = 3>
class ITK_EXPORT OpenCL1DComplexToComplexImageFilter :
    public FFT1DComplexToComplexImageFilter<TPixel,Dimension>
{
public:
  typedef OpenCL1DComplexToComplexImageFilter Self;
  typedef FFT1DComplexToComplexImageFilter<TPixel,Dimension> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Standard class typedefs.*/
  typedef typename Superclass::InputImageType InputImageType;
  typedef typename Superclass::OutputImageType OutputImageType;

  struct OpenCLComplexType
  {
    TPixel real;
    TPixel imag;
  };

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(OpenCL1DComplexToComplexImageFilter,
               FFT1DComplexToComplexImageFilter);


protected:
  OpenCL1DComplexToComplexImageFilter();
  virtual ~OpenCL1DComplexToComplexImageFilter()
  {
    if(m_PlanComputed)
      {
      clFFT_DestroyPlan(this->m_Plan);
      delete [] this->m_InputBuffer;
      delete [] this->m_OutputBuffer;
      }
    delete m_clQueue;
    delete m_clContext;
  }

  virtual void GenerateData();  // generates output from input

  ///** Method to check if an array dimension is legal for current OpenCL FFT */
  bool Legaldim(int n); 
private:
  OpenCL1DComplexToComplexImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  bool m_PlanComputed;
  clFFT_Plan m_Plan;
  unsigned int m_LastImageSize;
  OpenCLComplexType *m_InputBuffer;
  OpenCLComplexType *m_OutputBuffer;
  cl::Context *m_clContext;
  cl::CommandQueue *m_clQueue;

};

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkOpenCL1DComplexToComplexImageFilter.hxx"
#endif

#endif //__itkOpenCL1DComplexToComplexImageFilter_h
