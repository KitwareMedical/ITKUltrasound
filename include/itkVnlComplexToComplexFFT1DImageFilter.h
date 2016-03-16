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
#ifndef itkVnlComplexToComplexFFT1DImageFilter_h
#define itkVnlComplexToComplexFFT1DImageFilter_h

#include "itkComplexToComplexFFT1DImageFilter.h"
#include <complex>

namespace itk
{

/** \class VnlComplexToComplexFFT1DImageFilter
 *
 * \brief Perform the FFT along one dimension of an image using Vnl as a
 * backend.
 *
 * \ingroup Ultrasound
 */
template< typename TInputImage, typename TOutputImage >
class VnlComplexToComplexFFT1DImageFilter:
    public ComplexToComplexFFT1DImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef VnlComplexToComplexFFT1DImageFilter                           Self;
  typedef ComplexToComplexFFT1DImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                                          Pointer;
  typedef SmartPointer< const Self >                                    ConstPointer;

  typedef typename Superclass::InputImageType                           InputImageType;
  typedef typename Superclass::OutputImageType                          OutputImageType;
  typedef typename OutputImageType::RegionType                          OutputImageRegionType;

  typedef typename Superclass::TransformDirectionType                   TransformDirectionType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( VnlComplexToComplexFFT1DImageFilter, ComplexToComplexFFT1DImageFilter );

protected:
  VnlComplexToComplexFFT1DImageFilter() {}
  virtual ~VnlComplexToComplexFFT1DImageFilter() {}

  virtual void ThreadedGenerateData( const OutputImageRegionType&, ThreadIdType threadID ) ITK_OVERRIDE;

private:
  VnlComplexToComplexFFT1DImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVnlComplexToComplexFFT1DImageFilter.hxx"
#endif

#endif
