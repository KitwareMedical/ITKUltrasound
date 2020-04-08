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
#ifndef itkFrequencyDomain1DFilterFunction_h
#define itkFrequencyDomain1DFilterFunction_h

#include "itkObject.h"
#include "itkObjectFactory.h"

namespace itk
{
/** \class FrequencyDomain1DFilterFunction
 * \brief
 * Class to implment filter functions for FrequencyDomain1DImageFilter
 *
 * Supports caching of precomputed function values (SetUseCache) for applying
 * the function to multiple signals of the same length.
 * For the caching to work properly make sure this->Modified gets called in subclasses
 * whenever a parmater is changed that changes the function values.
 *
 * \ingroup FourierTransform
 * \ingroup Ultrasound
 */
class FrequencyDomain1DFilterFunction:
  public Object
{
public:

  /** Standard class type alias. */
  using Self = FrequencyDomain1DFilterFunction;
  using Superclass = Object;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  itkTypeMacro( FrequencyDomain1DFilterFunction, Object);
  itkNewMacro( Self );


  /** discrete frequency index such as retunr by FFTW, i.e.
   * 0 is DC component
   * i / SignalSize * 2*pi is the frequency,
   * i.e., first half positive frequencies and
   * second half negative frequencies in reverse order.
   */
  double EvaluateIndex(SizeValueType &i) const
    {
    if( m_UseCache )
      {
      //TODO: Check for out of bounds?
      return m_Cache[i];
      }
    else
      {
      return this->EvaluateFrequency( this->GetFrequency(i) );
      }
    }

  void SetSignalSize( const SizeValueType &size)
    {
    if( this->m_SignalSize != size )
      {
      this->m_SignalSize = size;
      if( this->m_UseCache )
        {
        this->m_Cache.resize( size );
        }
      this->Modified();
      }
    }

  SizeValueType GetSignalSize() const
    {
    return m_SignalSize;
    }

   itkSetMacro( UseCache, bool );
   itkGetMacro( UseCache, bool );

  /**
   * Override this function to implement a specific filter.
   * The input ranges from -1 to 1
   * Default is identity function.
   */
  virtual double EvaluateFrequency( double itkNotUsed( frequency ) ) const
    {
    return 1.0;
    }

  virtual void Modified( ) const override
    {
    //Force a chache update
    const_cast< FrequencyDomain1DFilterFunction *>(this)->UpdateCache();
    Superclass::Modified();
    }

protected:

  virtual void PrintSelf( std::ostream& os, Indent indent ) const override
    {
    Superclass::PrintSelf( os, indent );

    os << indent << "SignalSize: " << m_SignalSize << std::endl;
    os << indent << "UseCache: " << m_UseCache << std::endl;
    os << indent << "CacheSize: " << m_Cache.size() << std::endl;
    }

    FrequencyDomain1DFilterFunction()
      {
      m_UseCache = true;
      m_SignalSize = 0;
      }

private:

  double GetFrequency(SizeValueType &i) const
    {
    double f = (2.0 * i ) / m_SignalSize;
    if( f > 1.0 )
      {
      f = f - 2.0;
      }
    return f;
    }

  void UpdateCache()
    {
    if( this->m_UseCache )
      {
      for( SizeValueType i=0; i < m_Cache.size(); i++)
        {
        this->m_Cache[i] = this->EvaluateFrequency( this->GetFrequency(i) );
        }
      }
    }

  bool                  m_UseCache;
  std::vector< double > m_Cache;
  SizeValueType         m_SignalSize;
};

}
#endif // itkFrequencyDomain1DFilterFunction.h
