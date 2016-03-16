/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkFFTWComplexToComplexFFT1DImageFilter_h
#define itkFFTWComplexToComplexFFT1DImageFilter_h

#include "itkComplexToComplexFFT1DImageFilter.h"
#include "itkFFTWCommonExtended.h"

#include <vector>


namespace itk
{

/** \class FFTWComplexToComplexFFT1DImageFilter
 * \brief only do FFT along one dimension using FFTW as a backend.
 *
 * \ingroup Ultrasound
 */
template< typename TInputImage, typename TOutputImage >
class FFTWComplexToComplexFFT1DImageFilter :
    public ComplexToComplexFFT1DImageFilter< TInputImage, TOutputImage >
{
public:
  typedef FFTWComplexToComplexFFT1DImageFilter                             Self;
  typedef ComplexToComplexFFT1DImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                                          Pointer;
  typedef SmartPointer< const Self >                                    ConstPointer;

  /** Standard class typedefs.*/
  typedef typename Superclass::InputImageType  InputImageType;
  typedef typename Superclass::OutputImageType OutputImageType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  /**
   * the proxy type is a wrapper for the fftw API
   * since the proxy is only defined over double and float,
   * trying to use any other pixel type is inoperative, as
   * is trying to use double if only the float FFTW1D version is
   * configured in, or float if only double is configured.
   */
  typedef typename fftw::ComplexToComplexProxy< typename TInputImage::PixelType::value_type > FFTW1DProxyType;
  typedef typename std::vector< typename FFTW1DProxyType::PlanType >     PlanArrayType;
  typedef typename std::vector< typename FFTW1DProxyType::ComplexType* > PlanBufferPointerType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( FFTWComplexToComplexFFT1DImageFilter, ComplexToComplexFFT1DImageFilter );


protected:
  FFTWComplexToComplexFFT1DImageFilter(): m_PlanComputed( false ),
    m_LastImageSize( 0 )
  {}
  virtual ~FFTWComplexToComplexFFT1DImageFilter()
  {
  if ( m_PlanComputed )
    {
    this->DestroyPlans();
    }
  }

  virtual void BeforeThreadedGenerateData() ITK_OVERRIDE;
  virtual void ThreadedGenerateData( const OutputImageRegionType&, ThreadIdType threadID ) ITK_OVERRIDE;

private:
  FFTWComplexToComplexFFT1DImageFilter(const Self&); //purposely not implemented
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
#include "itkFFTWComplexToComplexFFT1DImageFilter.hxx"
#endif

#endif //itkFFTWComplexToComplexFFT1DImageFilter_h
