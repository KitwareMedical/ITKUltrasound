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
#ifndef itkComplexConjugateImageFilter_h
#define itkComplexConjugateImageFilter_h

#include "itkUnaryFunctorImageFilter.h"

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
namespace Function
{
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
} // end namespace Function


template <class TInputImage, class TOutputImage>
class ITK_TEMPLATE_EXPORT ComplexConjugateImageFilter :
    public
UnaryFunctorImageFilter<TInputImage,TOutputImage,
                        Function::ComplexConjugate<
  typename TInputImage::PixelType,
  typename TOutputImage::PixelType>   >
{
public:
  /** Standard class typedefs. */
  typedef ComplexConjugateImageFilter    Self;
  typedef UnaryFunctorImageFilter< TInputImage,TOutputImage,
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
