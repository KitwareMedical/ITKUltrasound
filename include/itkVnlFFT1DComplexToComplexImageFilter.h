/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVnlFFT1DComplexToComplexImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009-01-27 19:30:16 $
  Version:   $Revision: 1.8 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVnlFFT1DComplexToComplexImageFilter_h
#define __itkVnlFFT1DComplexToComplexImageFilter_h

#include "itkFFT1DComplexToComplexImageFilter.h"
#include <complex>

namespace itk
{

/** \class VnlFFT1DComplexToComplexImageFilter
 *
 * \brief Perform the FFT along one dimension of an image using Vnl as a
 * backend.
 *
 * \ingroup FFT1D
 */
template <class TPixel, unsigned int VDimension = 3>
class VnlFFT1DComplexToComplexImageFilter :
    public FFT1DComplexToComplexImageFilter<TPixel,VDimension>
{
public:
  /** Standard class typedefs. */
  typedef VnlFFT1DComplexToComplexImageFilter                 Self;
  typedef FFT1DComplexToComplexImageFilter<TPixel,VDimension> Superclass;
  typedef SmartPointer<Self>                                      Pointer;
  typedef SmartPointer<const Self>                                ConstPointer;

  typedef typename Superclass::InputImageType  InputImageType;
  typedef typename Superclass::OutputImageType OutputImageType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  typedef typename Superclass::TransformDirectionType TransformDirectionType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(VnlFFT1DComplexToComplexImageFilter,
               FFT1DComplexToComplexImageFilter);

protected:
  virtual void ThreadedGenerateData( const OutputImageRegionType&, ThreadIdType threadID );  // generates output from input

  VnlFFT1DComplexToComplexImageFilter() { }
  ~VnlFFT1DComplexToComplexImageFilter() { }

private:
  VnlFFT1DComplexToComplexImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVnlFFT1DComplexToComplexImageFilter.hxx"
#endif

#endif
