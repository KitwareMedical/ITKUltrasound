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
#ifndef itkFrequencyDomain1DFilterFunction_h
#define itkFrequencyDomain1DFilterFunction_h

#include "itkLightObject.h"

namespace itk
{
/** \class FrequencyDomain1DImageFilter
 * \brief
 * Class to implment filter functions for FrequencyDomain1DImageFilter
 * 
 * \ingroup FourierTransform
 * \ingroup Ultrasound
 */
class FrequencyDomain1DFilterFunction:
  public LightObject 
{
public:

  /** Standard class typedefs. */

  typedef FrequencyDomain1DFilterFunction      Self;
  typedef LightObject                          Superclass;
  typedef SmartPointer< Self >                 Pointer;
  typedef SmartPointer< const Self >           ConstPointer;

  itkTypeMacro( FrequencyDomain1DFilterFunction, LightObject);
  itkNewMacro( Self );


  /** Implement default identity function */
  virtual double Evaluate(double frequency){
    return frequency;
  };

};

}

#endif // itkFrequencyDomain1DFilterFunction.h
