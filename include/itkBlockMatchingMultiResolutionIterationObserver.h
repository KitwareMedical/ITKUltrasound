#ifndef __itkBlockMatchingMultiResolutionIterationObserver_h
#define __itkBlockMatchingMultiResolutionIterationObserver_h

#include "itkBlockMatchingMultiResolutionIterationCommand.h"

#include "itkStrainImageFilter.h"
#include "itkSplitComponentsImageFilter.h"

#include "itkImageFileWriter.h"

#include <fstream>
#include <string>

namespace itk
{
namespace BlockMatching
{

template< class TMultiResolutionMethod >
/** \class MultiresolutionIterationObserver
 * \brief Save status images at different resolutions. */
class MultiResolutionIterationObserver :
  public MultiResolutionIterationCommand< TMultiResolutionMethod >
{
public:
  typedef MultiResolutionIterationObserver                        Self;
  typedef MultiResolutionIterationCommand<TMultiResolutionMethod> Superclass;
  typedef SmartPointer<Self>                                      Pointer;

  itkNewMacro( Self );

  virtual void Execute(itk::Object *caller, const itk::EventObject & event)
    {
    Execute( (const itk::Object *)caller, event);
    }

  virtual void Execute(const itk::Object * object, const itk::EventObject & event);

  typedef TMultiResolutionMethod
    MultiResolutionMethodType;
  typedef typename MultiResolutionMethodType::FixedImageType
    FixedImageType;
  typedef typename MultiResolutionMethodType::MovingImageType
    MovingImageType;
  typedef typename MultiResolutionMethodType::DisplacementImageType
    DisplacementImageType;

  itkSetStringMacro( OutputFilePrefix );
  itkGetConstMacro( OutputFilePrefix, std::string );

protected:
  MultiResolutionIterationObserver();
  virtual ~MultiResolutionIterationObserver();

  // search region radius
  std::string m_OutputFilePrefix;

  std::ofstream m_CSVFile;

  typedef ImageFileWriter< FixedImageType > FixedImageWriterType;
  typename FixedImageWriterType::Pointer m_FixedImageWriter;

  typedef ImageFileWriter< MovingImageType > MovingImageWriterType;
  typename MovingImageWriterType::Pointer m_MovingImageWriter;

  typedef ImageFileWriter< DisplacementImageType > DisplacementWriterType;
  typename DisplacementWriterType::Pointer m_DisplacementWriter;

  typedef typename MultiResolutionMethodType::MetricImageType MetricImageType;
  typedef typename MetricImageType::PixelType                 MetricPixelType;
  typedef itk::SplitComponentsImageFilter< DisplacementImageType, MetricImageType,
          MetricImageType::ImageDimension > DisplacementComponentsFilterType;
  typename DisplacementComponentsFilterType::Pointer
    m_DisplacementComponentsFilter;
  typedef typename itk::ImageFileWriter< MetricImageType >
    DisplacementComponentsWriterType;
  typename DisplacementComponentsWriterType::Pointer
    m_DisplacementComponentsWriter;

  typedef itk::StrainImageFilter< DisplacementImageType, MetricPixelType,
          MetricPixelType > StrainFilterType;
  typename StrainFilterType::Pointer m_StrainFilter;
  typedef typename StrainFilterType::OutputImageType StrainImageType;

  typedef itk::ImageFileWriter< StrainImageType > StrainWriterType;
  typename StrainWriterType::Pointer m_StrainWriter;

  typedef itk::SplitComponentsImageFilter< StrainImageType, MetricImageType, 3 >
    StrainComponentsFilterType;
  typename StrainComponentsFilterType::Pointer m_StrainComponentsFilter;
  typedef itk::ImageFileWriter< MetricImageType > StrainComponentsWriterType;
  typename StrainComponentsWriterType::Pointer m_StrainComponentsWriter;

private:
  MultiResolutionIterationObserver( const Self& );
  void operator=( const Self & );
};

} // end namespace BlockMatching
} // end namespace itk

#include "itkBlockMatchingMultiResolutionIterationObserver.txx"

#endif
