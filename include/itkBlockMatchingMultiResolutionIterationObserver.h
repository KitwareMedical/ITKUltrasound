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
#ifndef itkBlockMatchingMultiResolutionIterationObserver_h
#define itkBlockMatchingMultiResolutionIterationObserver_h

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

/** \class MultiResolutionIterationObserver
 *
 * \brief Save status images at different resolutions.y
 *
 * \ingroup Ultrasound
 * */
template< typename TMultiResolutionMethod >
class MultiResolutionIterationObserver :
  public MultiResolutionIterationCommand< TMultiResolutionMethod >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(MultiResolutionIterationObserver);

  typedef MultiResolutionIterationObserver                        Self;
  typedef MultiResolutionIterationCommand<TMultiResolutionMethod> Superclass;
  typedef SmartPointer<Self>                                      Pointer;

  itkNewMacro( Self );

  void Execute(itk::Object *caller, const itk::EventObject & event) override
    {
    Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event) override;

  typedef TMultiResolutionMethod                                    MultiResolutionMethodType;
  typedef typename MultiResolutionMethodType::FixedImageType        FixedImageType;
  typedef typename MultiResolutionMethodType::MovingImageType       MovingImageType;
  typedef typename MultiResolutionMethodType::DisplacementImageType DisplacementImageType;

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
};

} // end namespace BlockMatching
} // end namespace itk

#include "itkBlockMatchingMultiResolutionIterationObserver.hxx"

#endif
