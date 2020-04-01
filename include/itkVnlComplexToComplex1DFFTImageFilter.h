/*=========================================================================
 *
 *  Copyright NumFOCUS
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
#ifndef itkVnlComplexToComplex1DFFTImageFilter_h
#define itkVnlComplexToComplex1DFFTImageFilter_h

#include "itkComplexToComplex1DFFTImageFilter.h"
#include <complex>

namespace itk
{

/** \class VnlComplexToComplex1DFFTImageFilter
 *
 * \brief Perform the FFT along one dimension of an image using Vnl as a
 * backend.
 *
 * \ingroup Ultrasound
 */
template< typename TInputImage, typename TOutputImage >
class ITK_TEMPLATE_EXPORT VnlComplexToComplex1DFFTImageFilter:
    public ComplexToComplex1DFFTImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef VnlComplexToComplex1DFFTImageFilter                           Self;
  typedef ComplexToComplex1DFFTImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                                          Pointer;
  typedef SmartPointer< const Self >                                    ConstPointer;

  typedef typename Superclass::InputImageType                           InputImageType;
  typedef typename Superclass::OutputImageType                          OutputImageType;
  typedef typename OutputImageType::RegionType                          OutputImageRegionType;

  typedef typename Superclass::TransformDirectionType                   TransformDirectionType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( VnlComplexToComplex1DFFTImageFilter, ComplexToComplex1DFFTImageFilter );

protected:
  VnlComplexToComplex1DFFTImageFilter() {}
  virtual ~VnlComplexToComplex1DFFTImageFilter() {}

  void GenerateData() override;

private:
  VnlComplexToComplex1DFFTImageFilter(const Self&) ITK_DELETED_FUNCTION;
  void operator=(const Self&) ITK_DELETED_FUNCTION;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVnlComplexToComplex1DFFTImageFilter.hxx"
#endif

#endif
