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
#ifndef itkSpecialCoordinatesImageToVTKStructuredGridFilter_h
#define itkSpecialCoordinatesImageToVTKStructuredGridFilter_h

#include "itkProcessObject.h"
#include "itkSimpleDataObjectDecorator.h"
#include "vtkStructuredGrid.h"
#include "vtkSmartPointer.h"

namespace itk
{

/** \class itkSpecialCoordinatesImageToVTKStructuredGridFilter
 *
 * \brief Convert an itk::SpecialCoordinatesImage to a vtkStructuredGrid.
 *
 * \ingroup Ultrasound
 */
template <typename TInputImage>
class ITK_TEMPLATE_EXPORT SpecialCoordinatesImageToVTKStructuredGridFilter : public ProcessObject
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(SpecialCoordinatesImageToVTKStructuredGridFilter);

  /** Standard class type alias. */
  using Self = SpecialCoordinatesImageToVTKStructuredGridFilter;
  using Superclass = ProcessObject;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkOverrideGetNameOfClassMacro(SpecialCoordinatesImageToVTKStructuredGridFilter);

  using InputImageType = TInputImage;

  using StructuredGridPointerType = vtkSmartPointer<vtkStructuredGrid>;
  using DecoratedStructuredGridPointerType = SimpleDataObjectDecorator<StructuredGridPointerType>;

  itkStaticConstMacro(ImageDimension, unsigned int, InputImageType::ImageDimension);


  using Superclass::SetInput;
  virtual void
  SetInput(const InputImageType * image);
  const InputImageType *
  GetInput() const;

  vtkStructuredGrid *
  GetOutput();
  itkGetDecoratedOutputMacro(StructuredGrid, StructuredGridPointerType);

protected:
  SpecialCoordinatesImageToVTKStructuredGridFilter();
  ~SpecialCoordinatesImageToVTKStructuredGridFilter() override = default;

  /** Make a DataObject of the correct type to be used as the specified output. */
  using DataObjectPointerArraySizeType = ProcessObject::DataObjectPointerArraySizeType;
  using Superclass::MakeOutput;
  virtual DataObjectPointer MakeOutput(DataObjectPointerArraySizeType) override;

  virtual void
  GenerateData() override;

private:
  vtkSmartPointer<vtkStructuredGrid> m_StructuredGrid;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkSpecialCoordinatesImageToVTKStructuredGridFilter.hxx"
#endif

#endif
