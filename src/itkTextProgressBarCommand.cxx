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
#include "itkTextProgressBarCommand.h"

#include "itkProcessObject.h"

#include <iostream>

namespace itk
{

TextProgressBarCommand
::TextProgressBarCommand():
  m_Progress(   "[>                                                  ]" )
{
}


void
TextProgressBarCommand
::Execute(itk::Object *caller, const itk::EventObject & event)
{
  this->Execute( (const itk::Object *)caller, event);
}


void
TextProgressBarCommand
::Execute(const itk::Object * object, const itk::EventObject & event)
{
  if( !(itk::ProgressEvent().CheckEvent( &event ) ) )
    {
    return;
    }

  itk::ProcessObject::ConstPointer process = dynamic_cast< const itk::ProcessObject * > ( object );
  if( !process )
    {
    return;
    }
  double progress = process->GetProgress();

  const unsigned int position = static_cast< unsigned int >( progress * 50 );

  unsigned int i;
  for( i = 0; i < position; i++ )
    {
    m_Progress[i+1] = '=';
    }

  if( position == 50 )
    {
    m_Progress[position+1] = '=';
    }
  else
    {
    m_Progress[position+1] = '>';
    }

  for( i = position + 1; i < 51; i++ )
    {
    m_Progress[i+1] = ' ';
    }

  std::cout << process->GetNameOfClass() << ": ";
  //std::cout << m_Progress << std::endl;
  // If you are in a console, this is nicer because you don't have repeats.  But
  // in ctest -V nothing shows up.
  //std::cout << m_Identifier << ": ";
  std::cout << m_Progress << '\r' << std::flush;

  if( position == 50 )
    {
    std::cout << std::endl;
    }
}

} // end namespace itk
