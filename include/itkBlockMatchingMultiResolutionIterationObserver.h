/*=========================================================================
 *
 *  Copyright NumFOCUS
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         https://www.apache.org/licenses/LICENSE-2.0.txt
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
 * \brief Save status images at different resolutions.
 *
 * \ingroup Ultrasound
 * */
template <typename TMultiResolutionMethod>
class MultiResolutionIterationObserver : public MultiResolutionIterationCommand<TMultiResolutionMethod>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(MultiResolutionIterationObserver);

  using Self = MultiResolutionIterationObserver;
  using Superclass = MultiResolutionIterationCommand<TMultiResolutionMethod>;
  using Pointer = SmartPointer<Self>;

  itkNewMacro(Self);

  using Superclass::Execute;
  void
  Execute(itk::Object * caller, const itk::EventObject & event) override;
  //{
  // Execute( (const itk::Object *)caller, event);
  //}

  // void Execute(const itk::Object * object, const itk::EventObject & event) override;

  using MultiResolutionMethodType = TMultiResolutionMethod;
  using FixedImageType = typename MultiResolutionMethodType::FixedImageType;
  using MovingImageType = typename MultiResolutionMethodType::MovingImageType;
  using DisplacementImageType = typename MultiResolutionMethodType::DisplacementImageType;

  itkSetStringMacro(OutputFilePrefix);
  itkGetConstMacro(OutputFilePrefix, std::string);

protected:
  MultiResolutionIterationObserver();
  ~MultiResolutionIterationObserver() override;

  // search region radius
  std::string m_OutputFilePrefix;

  std::ofstream m_CSVFile;

  using FixedImageWriterType = ImageFileWriter<FixedImageType>;
  typename FixedImageWriterType::Pointer m_FixedImageWriter;

  using MovingImageWriterType = ImageFileWriter<MovingImageType>;
  typename MovingImageWriterType::Pointer m_MovingImageWriter;

  using DisplacementWriterType = ImageFileWriter<DisplacementImageType>;
  typename DisplacementWriterType::Pointer m_DisplacementWriter;

  using MetricImageType = typename MultiResolutionMethodType::MetricImageType;
  using MetricPixelType = typename MetricImageType::PixelType;
  using DisplacementComponentsFilterType =
    itk::SplitComponentsImageFilter<DisplacementImageType, MetricImageType, MetricImageType::ImageDimension>;
  typename DisplacementComponentsFilterType::Pointer m_DisplacementComponentsFilter;
  using DisplacementComponentsWriterType = typename itk::ImageFileWriter<MetricImageType>;
  typename DisplacementComponentsWriterType::Pointer m_DisplacementComponentsWriter;

  using StrainFilterType = itk::StrainImageFilter<DisplacementImageType, MetricPixelType, MetricPixelType>;
  typename StrainFilterType::Pointer m_StrainFilter;
  using StrainImageType = typename StrainFilterType::OutputImageType;

  using StrainWriterType = itk::ImageFileWriter<StrainImageType>;
  typename StrainWriterType::Pointer m_StrainWriter;

  using StrainComponentsFilterType = itk::SplitComponentsImageFilter<StrainImageType, MetricImageType, 3>;
  typename StrainComponentsFilterType::Pointer m_StrainComponentsFilter;
  using StrainComponentsWriterType = itk::ImageFileWriter<MetricImageType>;
  typename StrainComponentsWriterType::Pointer m_StrainComponentsWriter;

private:
};

} // end namespace BlockMatching
} // end namespace itk

#include "itkBlockMatchingMultiResolutionIterationObserver.hxx"

#endif
