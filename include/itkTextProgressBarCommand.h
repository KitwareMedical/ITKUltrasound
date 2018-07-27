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
#ifndef itkTextProgressCommand_h
#define itkTextProgressCommand_h

#include "itkCommand.h"

#include <string>

namespace itk
{

/** \class TextProgressBarCommand
 *
 * \brief A simple command that outputs a text progress bar the associated filter.
 *
 * \ingroup Ultrasound
 * */
class TextProgressBarCommand:
  public Command
{
public:
  typedef TextProgressBarCommand Self;
  typedef Command                Superclass;
  typedef SmartPointer< Self >   Pointer;

  itkNewMacro( Self );

protected:
  TextProgressBarCommand();

  void Execute(itk::Object *caller, const itk::EventObject & event) ITK_OVERRIDE;

  void Execute(const itk::Object * object, const itk::EventObject & event) ITK_OVERRIDE;

  std::string m_Progress;
};

} // end namespace itk

#endif
