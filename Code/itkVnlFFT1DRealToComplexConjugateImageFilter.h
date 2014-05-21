/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVnlFFT1DRealToComplexConjugateImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009-01-27 19:30:16 $
  Version:   $Revision: 1.8 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVnlFFT1DRealToComplexConjugateImageFilter_h
#define __itkVnlFFT1DRealToComplexConjugateImageFilter_h

#include "itkFFT1DRealToComplexConjugateImageFilter.h"
#include <complex>

namespace itk
{

/** \class VnlFFT1DRealToComplexConjugateImageFilter
 *
 * \brief Perform the FFT along one dimension of an image using Vnl as a
 * backend.
 */
template <class TPixel, unsigned int VDimension = 3>
class VnlFFT1DRealToComplexConjugateImageFilter :
    public FFT1DRealToComplexConjugateImageFilter<TPixel,VDimension>
{
public:
  /** Standard class typedefs. */
  typedef VnlFFT1DRealToComplexConjugateImageFilter                 Self;
  typedef FFT1DRealToComplexConjugateImageFilter<TPixel,VDimension> Superclass;
  typedef SmartPointer<Self>                                        Pointer;
  typedef SmartPointer<const Self>                                  ConstPointer;

  typedef typename Superclass::InputImageType  InputImageType;
  typedef typename Superclass::OutputImageType OutputImageType;
  typedef typename OutputImageType::RegionType OutputImageRegionType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(VnlFFT1DRealToComplexConjugateImageFilter,
               FFT1DRealToComplexConjugateImageFilter);

protected:
  virtual void ThreadedGenerateData( const OutputImageRegionType&, ThreadIdType threadID );

  VnlFFT1DRealToComplexConjugateImageFilter() { }
  ~VnlFFT1DRealToComplexConjugateImageFilter() { }

  ///** Method to check if an array dimension is legal for PFA FFT */
  bool Legaldim(int n);

private:
  VnlFFT1DRealToComplexConjugateImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkVnlFFT1DRealToComplexConjugateImageFilter.txx"
#endif

#endif
