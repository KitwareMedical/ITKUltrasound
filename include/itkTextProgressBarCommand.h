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
#ifndef itkTextProgressBarCommand_h
#define itkTextProgressBarCommand_h

#include "itkCommand.h"

#include "UltrasoundExport.h"

#include <string>

namespace itk
{

/** \class TextProgressBarCommand
 *
 * \brief A simple command that outputs a text progress bar the associated filter.
 *
 * \ingroup Ultrasound
 * */
class Ultrasound_EXPORT TextProgressBarCommand:
  public Command
{
public:
  using Self = TextProgressBarCommand;
  using Superclass = Command;
  using Pointer = SmartPointer< Self >;

  itkNewMacro( Self );

protected:
  TextProgressBarCommand();

  void Execute(itk::Object *caller, const itk::EventObject & event) override;

  void Execute(const itk::Object * object, const itk::EventObject & event) override;

  std::string m_Progress;
};

} // end namespace itk

#endif
