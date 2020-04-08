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
#ifndef itkBlockMatchingMultiResolutionSearchRegionWriterCommand_h
#define itkBlockMatchingMultiResolutionSearchRegionWriterCommand_h

#include "itkBlockMatchingMultiResolutionIterationCommand.h"

#include "itkImageFileWriter.h"

namespace itk
{
namespace BlockMatching
{

/** \class MultiResolutionSearchRegionWriterCommand
 *
 * \brief Write the search region images size images for each level.
 *
 * \ingroup Ultrasound
 * */
template <typename TMultiResolutionMethod>
class MultiResolutionSearchRegionWriterCommand : public MultiResolutionIterationCommand<TMultiResolutionMethod>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(MultiResolutionSearchRegionWriterCommand);

  using Self = MultiResolutionSearchRegionWriterCommand;
  using Superclass = MultiResolutionIterationCommand<TMultiResolutionMethod>;
  using Pointer = SmartPointer<Self>;

  itkNewMacro(Self);

  void
  Execute(itk::Object * caller, const itk::EventObject & event) override
  {
    Execute((const itk::Object *)caller, event);
  }

  void
  Execute(const itk::Object * object, const itk::EventObject & event) override;

  using MultiResolutionMethodType = TMultiResolutionMethod;

  using SearchRegionImageSourceType = typename MultiResolutionMethodType::SearchRegionImageSourceType;
  using SearchRegionImageType = typename SearchRegionImageSourceType::OutputImageType;
  using SearchRegionType = typename SearchRegionImageType::PixelType;

  using SearchRegionImageComponentType = Image<unsigned short, 2>;
  using SearchRegionImageComponentPointer = SearchRegionImageComponentType::Pointer;

  using SearchRegionComponentWriterType = ImageFileWriter<SearchRegionImageComponentType>;

  itkSetStringMacro(OutputFilePrefix);
  itkGetConstMacro(OutputFilePrefix, std::string);

protected:
  MultiResolutionSearchRegionWriterCommand();
  ~MultiResolutionSearchRegionWriterCommand();

  SearchRegionImageComponentPointer        m_SearchRegionImageComponent;
  SearchRegionComponentWriterType::Pointer m_SearchRegionComponentWriter;

  std::string m_OutputFilePrefix;

private:
};

} // end namespace BlockMatching
} // end namespace itk

#include "itkBlockMatchingMultiResolutionSearchRegionWriterCommand.hxx"

#endif
