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
#ifndef itkSpecialCoordinatesImageToVTKStructuredGridFilter_hxx
#define itkSpecialCoordinatesImageToVTKStructuredGridFilter_hxx

#include "itkSpecialCoordinatesImageToVTKStructuredGridFilter.h"

namespace itk
{

template< typename TInputImage >
SpecialCoordinatesImageToVTKStructuredGridFilter< TInputImage >
::SpecialCoordinatesImageToVTKStructuredGridFilter()
{
  this->SetNumberOfRequiredInputs( 1 );
  m_StructuredGrid = vtkSmartPointer< vtkStructuredGrid >::New();
}


template< typename TInputImage >
SpecialCoordinatesImageToVTKStructuredGridFilter< TInputImage >
::~SpecialCoordinatesImageToVTKStructuredGridFilter()
{
}


template< typename TInputImage >
void
SpecialCoordinatesImageToVTKStructuredGridFilter< TInputImage >
::SetInput(const InputImageType *input)
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput( 0,
                                    const_cast< InputImageType * >( input ) );
}


template< typename TInputImage >
void
SpecialCoordinatesImageToVTKStructuredGridFilter< TInputImage >
::GenerateData()
{
}

} // end namespace itk

#endif
