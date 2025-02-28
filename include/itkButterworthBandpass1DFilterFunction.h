/*=========================================================================
 *
 *  Copyright NumFOCUS
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         https://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkButterworthBandpass1DFilterFunction_h
#define itkButterworthBandpass1DFilterFunction_h

#include "itkFrequencyDomain1DFilterFunction.h"
#include "itkMath.h"

namespace itk
{
class ButterworthBandpass1DFilterFunction : public FrequencyDomain1DFilterFunction
{
public:
  /** Standard class type alias. */

  using Self = ButterworthBandpass1DFilterFunction;
  using Superclass = FrequencyDomain1DFilterFunction;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  itkOverrideGetNameOfClassMacro(ButterworthBandpass1DFilterFunction);
  itkNewMacro(Self);

  double
  EvaluateFrequency(double f) const override
  {
    // highpass
    double x = 1.0;
    if (m_LowerFrequency > 0.0)
    {
      x = 1.0 - 1.0 / (1.0 + std::pow(f / m_LowerFrequency, 2 * m_Order));
    }
    // lowpass
    double y = 1.0;
    if (m_UpperFrequency > 0.0)
    {
      y = 1.0 / (1.0 + std::pow(f / m_UpperFrequency, 2 * m_Order));
    }
    return x * y;
  };

  itkSetMacro(UpperFrequency, double);
  itkGetMacro(UpperFrequency, double);

  itkSetMacro(LowerFrequency, double);
  itkGetMacro(LowerFrequency, double);

  itkSetMacro(Order, int);
  itkGetMacro(Order, int);


protected:

    ButterworthBandpass1DFilterFunction()
  {
    m_Order = 1;
    m_LowerFrequency = 0.0;
    m_UpperFrequency = 1.0;
  }

  virtual void
  PrintSelf(std::ostream & os, Indent indent) const override
  {
    Superclass::PrintSelf(os, indent);

    os << indent << "LowerFrequency: " << m_LowerFrequency << std::endl;
    os << indent << "UpperFrequency: " << m_UpperFrequency << std::endl;
    os << indent << "Order: " << m_Order << std::endl;
  }

private:
  int    m_Order;
  double m_LowerFrequency;
  double m_UpperFrequency;
};

} // namespace itk
#endif // ButterworthBandpass1DFilterFunction.h
