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
#ifndef itkFFTWInverse1DFFTImageFilter_h
#define itkFFTWInverse1DFFTImageFilter_h

#include "itkInverse1DFFTImageFilter.h"
#include "itkFFTWCommonExtended.h"

#include <vector>

namespace itk
{
/** \class FFTWInverse1DFFTImageFilter
 * \brief only do FFT along one dimension using FFTW as a backend.
 *
 * \ingroup Ultrasound
 */

template< typename TInputImage, typename TOutputImage=Image< typename NumericTraits< typename TInputImage::PixelType >::ValueType, TInputImage::ImageDimension > >
class ITK_TEMPLATE_EXPORT FFTWInverse1DFFTImageFilter:
  public Inverse1DFFTImageFilter< TInputImage, TOutputImage >
{
public:
  typedef FFTWInverse1DFFTImageFilter                          Self;
  typedef Inverse1DFFTImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                                 Pointer;
  typedef SmartPointer< const Self >                           ConstPointer;

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
  typedef typename fftw::ComplexToComplexProxy< typename TOutputImage::PixelType > FFTW1DProxyType;
  typedef typename std::vector< typename FFTW1DProxyType::PlanType >               PlanArrayType;
  typedef typename std::vector< typename FFTW1DProxyType::ComplexType* >           PlanBufferPointerType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(FFTWInverse1DFFTImageFilter, Inverse1DFFTImageFilter);


protected:
  FFTWInverse1DFFTImageFilter(): m_PlanComputed( false ),
    m_LastImageSize( 0 )
  {}
  virtual ~FFTWInverse1DFFTImageFilter()
  {
  if ( m_PlanComputed )
    {
    this->DestroyPlans();
    }
  }

  virtual void BeforeThreadedGenerateData() ITK_OVERRIDE;
  virtual void ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread, ThreadIdType threadID ) ITK_OVERRIDE;

private:
  FFTWInverse1DFFTImageFilter(const Self&) ITK_DELETED_FUNCTION;
  void operator=(const Self&) ITK_DELETED_FUNCTION;

  /** Destroy FFTW Plans and associated buffers. */
  void DestroyPlans();

  bool                  m_PlanComputed;
  PlanArrayType         m_PlanArray;
  unsigned int          m_LastImageSize;
  PlanBufferPointerType m_InputBufferArray;
  PlanBufferPointerType m_OutputBufferArray;
};

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFFTWInverse1DFFTImageFilter.hxx"
#endif

#endif //itkFFTWInverse1DFFTImageFilter_h
