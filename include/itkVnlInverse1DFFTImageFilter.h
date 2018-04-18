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
#ifndef itkVnlInverse1DFFTImageFilter_h
#define itkVnlInverse1DFFTImageFilter_h

#include "itkInverse1DFFTImageFilter.h"
#include <complex>

namespace itk
{

/** \class VnlInverse1DFFTImageFilter
 *
 * \brief Perform the FFT along one dimension of an image using Vnl as a
 * backend.
 *
 * \ingroup Ultrasound
 */
template< typename TInputImage, typename TOutputImage=Image< typename NumericTraits< typename TInputImage::PixelType >::ValueType, TInputImage::ImageDimension > >
class ITK_TEMPLATE_EXPORT VnlInverse1DFFTImageFilter:
  public Inverse1DFFTImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef VnlInverse1DFFTImageFilter                           Self;
  typedef Inverse1DFFTImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                                 Pointer;
  typedef SmartPointer< const Self >                           ConstPointer;

  typedef typename Superclass::InputImageType                  InputImageType;
  typedef typename Superclass::OutputImageType                 OutputImageType;
  typedef typename OutputImageType::RegionType                 OutputImageRegionType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( VnlInverse1DFFTImageFilter, Inverse1DFFTImageFilter );

protected:
  void ThreadedGenerateData( const OutputImageRegionType&, ThreadIdType threadID ) ITK_OVERRIDE;

  VnlInverse1DFFTImageFilter() { }
  virtual ~VnlInverse1DFFTImageFilter() { }

private:
  VnlInverse1DFFTImageFilter(const Self&) ITK_DELETE_FUNCTION;
  void operator=(const Self&) ITK_DELETE_FUNCTION;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVnlInverse1DFFTImageFilter.hxx"
#endif

#endif
