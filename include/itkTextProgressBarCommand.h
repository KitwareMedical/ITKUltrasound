#ifndef __itkTextProgressCommand_h
#define __itkTextProgressCommand_h

#include "itkCommand.h"

#include <string>

namespace itk
{

/** A simple command that outputs a text progress bar the associated filter. */
class TextProgressBarCommand:
  public Command
{
public:
  typedef TextProgressBarCommand  Self;
  typedef Command              Superclass;
  typedef itk::SmartPointer< Self > Pointer;
  itkNewMacro( Self );

protected:
  TextProgressBarCommand();

  void Execute(itk::Object *caller, const itk::EventObject & event);

  void Execute(const itk::Object * object, const itk::EventObject & event);

  std::string m_Progress;
};

} // end namespace itk

#endif
