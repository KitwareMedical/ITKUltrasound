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
#ifndef itkVnlForward1DFFTImageFilter_h
#define itkVnlForward1DFFTImageFilter_h

#include "itkForward1DFFTImageFilter.h"
#include <complex>

namespace itk
{

/** \class VnlForward1DFFTImageFilter
 *
 * \brief Perform the FFT along one dimension of an image using Vnl as a
 * backend.
 *
 * \ingroup Ultrasound
 */
template< typename TInputImage, typename TOutputImage=Image< std::complex< typename TInputImage::PixelType >, TInputImage::ImageDimension > >
class ITK_TEMPLATE_EXPORT VnlForward1DFFTImageFilter :
    public Forward1DFFTImageFilter< TInputImage, TOutputImage >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(VnlForward1DFFTImageFilter);

  /** Standard class typedefs. */
  typedef VnlForward1DFFTImageFilter                           Self;
  typedef Forward1DFFTImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                                 Pointer;
  typedef SmartPointer< const Self >                           ConstPointer;

  typedef typename Superclass::InputImageType                  InputImageType;
  typedef typename Superclass::OutputImageType                 OutputImageType;
  typedef typename OutputImageType::RegionType                 OutputImageRegionType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( VnlForward1DFFTImageFilter, Forward1DFFTImageFilter );

protected:
  virtual void ThreadedGenerateData( const OutputImageRegionType&, ThreadIdType threadID ) ITK_OVERRIDE;

  VnlForward1DFFTImageFilter() { }
  ~VnlForward1DFFTImageFilter() { }

private:
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVnlForward1DFFTImageFilter.hxx"
#endif

#endif
