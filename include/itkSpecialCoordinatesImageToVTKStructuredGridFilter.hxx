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
#ifndef itkSpecialCoordinatesImageToVTKStructuredGridFilter_hxx
#define itkSpecialCoordinatesImageToVTKStructuredGridFilter_hxx

#include "itkSpecialCoordinatesImageToVTKStructuredGridFilter.h"

#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkPixelTraits.h"
#include "itkProgressReporter.h"

#include "vtkPointData.h"
#include "vtkNew.h"
#include "vtkPoints.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkLongArray.h"
#include "vtkUnsignedLongArray.h"
#include "vtkIntArray.h"
#include "vtkUnsignedIntArray.h"
#include "vtkShortArray.h"
#include "vtkUnsignedShortArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkSignedCharArray.h"

#include <cstring>

namespace itk
{

template <typename TInputImage>
SpecialCoordinatesImageToVTKStructuredGridFilter<TInputImage>::SpecialCoordinatesImageToVTKStructuredGridFilter()
{
  this->SetNumberOfRequiredInputs(1);

  this->SetPrimaryOutputName("StructuredGrid");
  this->SetNumberOfRequiredOutputs(1);
  this->ReleaseDataBeforeUpdateFlagOff();

  typename DecoratedStructuredGridPointerType::Pointer decoratedStructuredGrid =
    itkDynamicCastInDebugMode<DecoratedStructuredGridPointerType *>(this->MakeOutput(0).GetPointer());
  decoratedStructuredGrid->Set(m_StructuredGrid);
  this->ProcessObject::SetNthOutput(0, decoratedStructuredGrid);
}


template <typename TInputImage>
DataObject::Pointer
SpecialCoordinatesImageToVTKStructuredGridFilter<TInputImage>::MakeOutput(DataObjectPointerArraySizeType output)
{
  switch (output)
  {
    case 0:
    {
      m_StructuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
      DecoratedStructuredGridPointerType::Pointer decoratedStructuredGrid = DecoratedStructuredGridPointerType::New();
      decoratedStructuredGrid->Set(m_StructuredGrid);
      return decoratedStructuredGrid.GetPointer();
    }
    default:
      itkExceptionMacro("MakeOutput request for an output number larger than the expected number of outputs.");
      return nullptr;
  }
}


template <typename TInputImage>
void
SpecialCoordinatesImageToVTKStructuredGridFilter<TInputImage>::SetInput(const InputImageType * input)
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput(0, const_cast<InputImageType *>(input));
}


template <typename TInputImage>
const typename SpecialCoordinatesImageToVTKStructuredGridFilter<TInputImage>::InputImageType *
SpecialCoordinatesImageToVTKStructuredGridFilter<TInputImage>::GetInput() const
{
  return itkDynamicCastInDebugMode<const TInputImage *>(this->GetPrimaryInput());
}


template <typename TInputImage>
vtkStructuredGrid *
SpecialCoordinatesImageToVTKStructuredGridFilter<TInputImage>::GetOutput()
{
  return this->m_StructuredGrid.GetPointer();
}


template <typename TInputImage>
void
SpecialCoordinatesImageToVTKStructuredGridFilter<TInputImage>::GenerateData()
{
  this->m_StructuredGrid->PrepareForNewData();
  // Const-cast because there is not (yet) a const version of GetBufferPointer()
  InputImageType * inputImage = const_cast<InputImageType *>(this->GetInput());

  const typename InputImageType::RegionType region = inputImage->GetBufferedRegion();
  const typename InputImageType::SizeType   size = region.GetSize();
  int                                       dims[3];
  dims[2] = 1;
  for (unsigned int ii = 0; ii < ImageDimension; ++ii)
  {
    dims[ii] = size[ii];
  }
  this->m_StructuredGrid->SetDimensions(dims);

  vtkNew<vtkPoints> gridPoints;
  gridPoints->SetDataTypeToDouble();
  gridPoints->SetNumberOfPoints(region.GetNumberOfPixels());

  ProgressReporter progress(this, 0, region.GetNumberOfPixels());

  ImageRegionConstIteratorWithIndex<InputImageType> imageIt(inputImage, region);
  imageIt.GoToBegin();
  for (SizeValueType ii = 0; !imageIt.IsAtEnd(); ++imageIt, ++ii)
  {
    const typename InputImageType::IndexType & index = imageIt.GetIndex();
    typename InputImageType::PointType         point;
    inputImage->TransformIndexToPhysicalPoint(index, point);
    gridPoints->SetPoint(ii, point.GetVnlVector().data_block());
    progress.CompletedPixel();
  }


  vtkDataArray * pointDataArray;
  using PixelType = typename InputImageType::PixelType;
  using ScalarType = typename PixelTraits<PixelType>::ValueType;
  const unsigned int  components = inputImage->GetNumberOfComponentsPerPixel();
  const SizeValueType tuples = region.GetNumberOfPixels();
  if (typeid(ScalarType) == typeid(double))
  {
    vtkDoubleArray * array = vtkDoubleArray::New();
    array->SetNumberOfComponents(components);
    array->SetNumberOfTuples(tuples);
    std::memcpy(array->GetVoidPointer(0),
                reinterpret_cast<void *>(inputImage->GetPixelContainer()->GetBufferPointer()),
                tuples * components * sizeof(double));
    pointDataArray = array;
  }
  else if (typeid(ScalarType) == typeid(float))
  {
    vtkFloatArray * array = vtkFloatArray::New();
    array->SetNumberOfComponents(components);
    array->SetNumberOfTuples(tuples);
    std::memcpy(array->GetVoidPointer(0),
                reinterpret_cast<void *>(inputImage->GetPixelContainer()->GetBufferPointer()),
                tuples * components * sizeof(float));
    pointDataArray = array;
  }
  else if (typeid(ScalarType) == typeid(long))
  {
    vtkLongArray * array = vtkLongArray::New();
    array->SetNumberOfComponents(components);
    array->SetNumberOfTuples(tuples);
    std::memcpy(array->GetVoidPointer(0),
                reinterpret_cast<void *>(inputImage->GetPixelContainer()->GetBufferPointer()),
                tuples * components * sizeof(long));
    pointDataArray = array;
  }
  else if (typeid(ScalarType) == typeid(unsigned long))
  {
    vtkUnsignedLongArray * array = vtkUnsignedLongArray::New();
    array->SetNumberOfComponents(components);
    array->SetNumberOfTuples(tuples);
    std::memcpy(array->GetVoidPointer(0),
                reinterpret_cast<void *>(inputImage->GetPixelContainer()->GetBufferPointer()),
                tuples * components * sizeof(unsigned long));
    pointDataArray = array;
  }
  else if (typeid(ScalarType) == typeid(int))
  {
    vtkIntArray * array = vtkIntArray::New();
    array->SetNumberOfComponents(components);
    array->SetNumberOfTuples(tuples);
    std::memcpy(array->GetVoidPointer(0),
                reinterpret_cast<void *>(inputImage->GetPixelContainer()->GetBufferPointer()),
                tuples * components * sizeof(int));
    pointDataArray = array;
  }
  else if (typeid(ScalarType) == typeid(unsigned int))
  {
    vtkUnsignedIntArray * array = vtkUnsignedIntArray::New();
    array->SetNumberOfComponents(components);
    array->SetNumberOfTuples(tuples);
    std::memcpy(array->GetVoidPointer(0),
                reinterpret_cast<void *>(inputImage->GetPixelContainer()->GetBufferPointer()),
                tuples * components * sizeof(unsigned int));
    pointDataArray = array;
  }
  else if (typeid(ScalarType) == typeid(short))
  {
    vtkShortArray * array = vtkShortArray::New();
    array->SetNumberOfComponents(components);
    array->SetNumberOfTuples(tuples);
    std::memcpy(array->GetVoidPointer(0),
                reinterpret_cast<void *>(inputImage->GetPixelContainer()->GetBufferPointer()),
                tuples * components * sizeof(short));
    pointDataArray = array;
  }
  else if (typeid(ScalarType) == typeid(unsigned short))
  {
    vtkUnsignedShortArray * array = vtkUnsignedShortArray::New();
    array->SetNumberOfComponents(components);
    array->SetNumberOfTuples(tuples);
    std::memcpy(array->GetVoidPointer(0),
                reinterpret_cast<void *>(inputImage->GetPixelContainer()->GetBufferPointer()),
                tuples * components * sizeof(unsigned short));
    pointDataArray = array;
  }
  else if (typeid(ScalarType) == typeid(unsigned char))
  {
    vtkUnsignedCharArray * array = vtkUnsignedCharArray::New();
    array->SetNumberOfComponents(components);
    array->SetNumberOfTuples(tuples);
    std::memcpy(array->GetVoidPointer(0),
                reinterpret_cast<void *>(inputImage->GetPixelContainer()->GetBufferPointer()),
                tuples * components * sizeof(unsigned char));
    pointDataArray = array;
  }
  else if (typeid(ScalarType) == typeid(signed char))
  {
    vtkSignedCharArray * array = vtkSignedCharArray::New();
    array->SetNumberOfComponents(components);
    array->SetNumberOfTuples(tuples);
    std::memcpy(array->GetVoidPointer(0),
                reinterpret_cast<void *>(inputImage->GetPixelContainer()->GetBufferPointer()),
                tuples * components * sizeof(signed char));
    pointDataArray = array;
  }
  else
  {
    itkExceptionMacro(<< "Type currently not supported");
  }

  if (components == 1)
  {
    pointDataArray->SetName("Scalars");
    this->m_StructuredGrid->GetPointData()->SetScalars(pointDataArray);
  }
  else
  {
    pointDataArray->SetName("Vectors");
    this->m_StructuredGrid->GetPointData()->SetVectors(pointDataArray);
  }
  pointDataArray->Delete();

  this->m_StructuredGrid->SetPoints(gridPoints.GetPointer());
}

} // end namespace itk

#endif
