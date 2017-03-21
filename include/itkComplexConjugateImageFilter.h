/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkComplexConjugateImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2009-04-01 14:36:10 $
  Version:   $Revision: 1.8 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkComplexConjugateImageFilter_h
#define __itkComplexConjugateImageFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include "vnl/vnl_math.h"

#include <complex>

namespace itk
{
  
/** \class ComplexConjugateImageFilter
 * \brief Computes pixel-wise the complex conjugate of a complex image.
 * 
 * \ingroup IntensityImageFilters  Multithreaded
 *
 * \ingroup Ultrasound
 */
namespace Function {  
  
template< class TInput, class TOutput>
class ComplexConjugate
{
public:
  ComplexConjugate() {}
  ~ComplexConjugate() {}
  bool operator!=( const ComplexConjugate & ) const
    {
    return false;
    }
  bool operator==( const ComplexConjugate & other ) const
    {
    return !(*this != other);
    }
  inline TOutput operator()( const TInput & A ) const
    {
    return (TOutput)( std::conj( A ) );
    }
}; 
}

template <class TInputImage, class TOutputImage>
class ITK_EXPORT ComplexConjugateImageFilter :
    public
UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                        Function::ComplexConjugate< 
  typename TInputImage::PixelType, 
  typename TOutputImage::PixelType>   >
{
public:
  /** Standard class typedefs. */
  typedef ComplexConjugateImageFilter  Self;
  typedef UnaryFunctorImageFilter<
    TInputImage,TOutputImage, 
    Function::ComplexConjugate< typename TInputImage::PixelType, 
                                  typename TOutputImage::PixelType> >
                                         Superclass;
  typedef SmartPointer<Self>             Pointer;
  typedef SmartPointer<const Self>       ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(ComplexConjugateImageFilter, 
               UnaryFunctorImageFilter);

  typedef typename TInputImage::PixelType                     InputPixelType;
  typedef typename TOutputImage::PixelType                    OutputPixelType;
  typedef typename NumericTraits< InputPixelType >::ValueType InputPixelValueType;

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputConvertibleToOutputCheck,
    (Concept::Convertible<InputPixelValueType, OutputPixelType>));
  /** End concept checking */
#endif

protected:
  ComplexConjugateImageFilter() {}
  virtual ~ComplexConjugateImageFilter() {}

private:
  ComplexConjugateImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk


#endif
