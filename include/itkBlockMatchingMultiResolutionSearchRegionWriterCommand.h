#ifndef __itkBlockMatchingMultiResolutionSearchRegionWriterCommand_h
#define __itkBlockMatchingMultiResolutionSearchRegionWriterCommand_h

#include "itkBlockMatchingMultiResolutionIterationCommand.h"

#include "itkImageFileWriter.h"

namespace itk
{
namespace BlockMatching
{

/** \class MultiResolutionSearchRegionWriterCommand
 *
 * Write the search region images size images for each level. */
template < class TMultiResolutionMethod >
class MultiResolutionSearchRegionWriterCommand :
  public MultiResolutionIterationCommand< TMultiResolutionMethod >
{
public:
  typedef MultiResolutionSearchRegionWriterCommand                Self;
  typedef MultiResolutionIterationCommand<TMultiResolutionMethod> Superclass;
  typedef SmartPointer<Self>                                      Pointer;

  itkNewMacro( Self );

  virtual void Execute(itk::Object *caller, const itk::EventObject & event)
    {
    Execute( (const itk::Object *)caller, event);
    }

  virtual void Execute(const itk::Object * object, const itk::EventObject & event);

  typedef TMultiResolutionMethod MultiResolutionMethodType;

  typedef typename MultiResolutionMethodType::SearchRegionImageSourceType SearchRegionImageSourceType;
  typedef typename SearchRegionImageSourceType::OutputImageType           SearchRegionImageType;
  typedef typename SearchRegionImageType::PixelType                       SearchRegionType;

  typedef Image< unsigned short, 2 > SearchRegionImageComponentType;
  typedef SearchRegionImageComponentType::Pointer SearchRegionImageComponentPointer;

  typedef ImageFileWriter< SearchRegionImageComponentType > SearchRegionComponentWriterType;

  itkSetStringMacro( OutputFilePrefix );
  itkGetConstMacro( OutputFilePrefix, std::string );

protected:
  MultiResolutionSearchRegionWriterCommand();
  ~MultiResolutionSearchRegionWriterCommand();

  SearchRegionImageComponentPointer        m_SearchRegionImageComponent;
  SearchRegionComponentWriterType::Pointer m_SearchRegionComponentWriter;

  std::string m_OutputFilePrefix;

private:
  MultiResolutionSearchRegionWriterCommand( const Self& );
  void operator=( const Self & );
};

} // end namespace BlockMatching
} // end namespace itk

#include "itkBlockMatchingMultiResolutionSearchRegionWriterCommand.txx"

#endif
