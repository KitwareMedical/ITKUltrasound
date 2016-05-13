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
#ifndef itkSpecialCoordinatesImageToVTKStructuredGridFilter_h
#define itkSpecialCoordinatesImageToVTKStructuredGridFilter_h

#include "itkProcessObject.h"
#include "vtkStructuredGrid.h"

namespace itk
{

/** \class itkSpecialCoordinatesImageToVTKStructuredGridFilter
 *
 * \brief Convert an itk::SpecialCoordinatesImage to a vtkStructuredGrid.
 *
 * \ingroup Ultrasound
 */
template< typename TInputImage >
class SpecialCoordinatesImageToVTKStructuredGridFilter: public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef SpecialCoordinatesImageToVTKStructuredGridFilter Self;
  typedef ProcessObject                                    Superclass;
  typedef SmartPointer<Self>                               Pointer;
  typedef SmartPointer<const Self>                         ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SpecialCoordinatesImageToVTKStructuredGridFilter, ProcessObject);

};
} // end namespace itk

#endif
