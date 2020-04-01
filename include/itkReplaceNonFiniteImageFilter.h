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
#ifndef itkReplaceNonFiniteImageFilter_h
#define itkReplaceNonFiniteImageFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include "itkMath.h"

namespace itk
{
namespace Functor
{
/**
 * \class ReplaceNonFinite
 * \brief Replaces a non-finite value with a finite one
 * \ingroup Ultrasound
 */
template< typename TInput, typename TOutput >
class ReplaceNonFinite
{
public:
  typedef TInput  InputType;
  typedef TOutput OutputType;

  ReplaceNonFinite(): m_ReplacementValue( NumericTraits< OutputType >::ZeroValue() ) {}
  ~ReplaceNonFinite() {}

  const OutputType & GetReplacementValue() const { return this->m_ReplacementValue; }
  void SetReplacementValue( const OutputType & value ) { this->m_ReplacementValue = value; }

  bool operator!=(const ReplaceNonFinite & other) const
  {
    return (other.GetReplacementValue() != this->m_ReplacementValue);
  }

  bool operator==(const ReplaceNonFinite & other) const
  {
    return !( *this != other );
  }

  inline OutputType operator()(const InputType & input) const
  {
    if( Math::isfinite( input ) )
      {
      return static_cast< OutputType >( input );
      }
    else
      {
      return this->m_ReplacementValue;
      }
  }

private:
  OutputType m_ReplacementValue;
};
}


/** \class ReplaceNonFiniteImageFilter
 * \brief Replace -inf, inf, and NaN values in an image.
 *
 * -inf, inf, and NaN values are replaced by the ReplacementValue.
 *
 * \ingroup Ultrasound
 */
template< typename TInputImage, typename TOutputImage = TInputImage >
class ReplaceNonFiniteImageFilter:
  public UnaryFunctorImageFilter< TInputImage, TOutputImage,
    Functor::ReplaceNonFinite< typename TInputImage::PixelType, typename TOutputImage::PixelType > >
{
public:
  /** Standard class typedefs. */
  typedef ReplaceNonFiniteImageFilter Self;
  typedef UnaryFunctorImageFilter< TInputImage, TOutputImage,
                                   Functor::ReplaceNonFinite< typename TInputImage::PixelType,
                                                  typename TOutputImage::PixelType > >  Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(ReplaceNonFiniteImageFilter, UnaryFunctorImageFilter);

protected:
  ReplaceNonFiniteImageFilter() {}
  virtual ~ReplaceNonFiniteImageFilter() {}

private:
  ReplaceNonFiniteImageFilter(const Self &) ITK_DELETED_FUNCTION;
  void operator=(const Self &) ITK_DELETED_FUNCTION;
};

} // end namespace itk

#endif // itkReplaceNonFiniteImageFilter_h
