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
/*=========================================================================
 *
 *  Portions of this file are subject to the VTK Toolkit Version 3 copyright.
 *
 *  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 *
 *  For complete copyright, license and disclaimer of warranty information
 *  please refer to the NOTICE file at the top of the ITK source tree.
 *
 *=========================================================================*/
#ifndef itkSliceSeriesSpecialCoordinatesImage_hxx
#define itkSliceSeriesSpecialCoordinatesImage_hxx

namespace itk
{

template <typename TSliceImage, typename TTransform, typename TPixel, unsigned int VDimension>
SliceSeriesSpecialCoordinatesImage<TSliceImage, TTransform, TPixel, VDimension>::SliceSeriesSpecialCoordinatesImage()
{
  this->m_SliceImage = SliceImageType::New();
  this->m_SliceTransforms = SliceTransformsType::New();
  this->m_SliceInverseTransforms = SliceTransformsType::New();
}


template <typename TSliceImage, typename TTransform, typename TPixel, unsigned int VDimension>
void
SliceSeriesSpecialCoordinatesImage<TSliceImage, TTransform, TPixel, VDimension>::SetLargestPossibleRegion(
  const RegionType & region)
{
  Superclass::SetLargestPossibleRegion(region);
  const SizeType & largestSize = this->GetLargestPossibleRegion().GetSize();
  this->m_SliceTransforms->Reserve(largestSize[ImageDimension - 1] + 1);
  this->m_SliceInverseTransforms->Reserve(largestSize[ImageDimension - 1] + 1);
}


template <typename TSliceImage, typename TTransform, typename TPixel, unsigned int VDimension>
void
SliceSeriesSpecialCoordinatesImage<TSliceImage, TTransform, TPixel, VDimension>::CopyInformation(
  const DataObject * data)
{
  Superclass::CopyInformation(data);

  if (data)
  {
    // Attempt to cast data to an ImageBase
    const SliceSeriesSpecialCoordinatesImage<TSliceImage, TTransform, TPixel, VDimension> * const sliceSeries =
      dynamic_cast<const SliceSeriesSpecialCoordinatesImage<TSliceImage, TTransform, TPixel, VDimension> *>(data);

    if (sliceSeries != nullptr)
    {
      // Copy the meta data for this data type
      typename SliceImageType::Pointer sliceImage = SliceImageType::New();
      sliceImage->CopyInformation(sliceSeries->GetSliceImage());
      this->SetSliceImage(sliceImage);

      const SizeValueType numberOfTransforms = sliceSeries->GetLargestPossibleRegion().GetSize()[ImageDimension - 1];
      for (SizeValueType transformIndex = 0; transformIndex < numberOfTransforms; ++transformIndex)
      {
        const TransformType * transform = sliceSeries->GetSliceTransform(transformIndex);
        if (transform != nullptr)
        {
          typename TransformType::Pointer transformClone = transform->Clone();
          this->SetSliceTransform(transformIndex, transformClone);
        }
      }
    }
    else
    {
      // pointer could not be cast back down
      itkExceptionMacro(<< "itk::SliceSeriesSpecialCoordinatesImage::CopyInformation() cannot cast "
                        << typeid(data).name() << " to " << typeid(const SliceSeriesSpecialCoordinatesImage *).name());
    }
  }
}


template <typename TSliceImage, typename TTransform, typename TPixel, unsigned int VDimension>
void
SliceSeriesSpecialCoordinatesImage<TSliceImage, TTransform, TPixel, VDimension>::SetSliceTransform(
  IndexValueType  sliceIndex,
  TransformType * transform)
{
  const RegionType &   largestRegion = this->GetLargestPossibleRegion();
  const IndexValueType largestIndex = largestRegion.GetIndex(ImageDimension - 1);
  const SizeValueType  transformsIndex = static_cast<SizeValueType>(sliceIndex - largestIndex);
  if (this->m_SliceTransforms->Size() > transformsIndex &&
      this->m_SliceTransforms->GetElement(transformsIndex).GetPointer() == transform)
  {
    return;
  }
  this->m_SliceTransforms->SetElement(transformsIndex, transform);
  typename TransformType::Pointer inverse = TransformType::New();
  if (transform->GetInverse(inverse))
  {
    this->m_SliceInverseTransforms->SetElement(transformsIndex, inverse);
  }
  else
  {
    itkExceptionMacro("Could not get inverse for transform: " << transform);
  }
  this->Modified();
}


template <typename TSliceImage, typename TTransform, typename TPixel, unsigned int VDimension>
const typename SliceSeriesSpecialCoordinatesImage<TSliceImage, TTransform, TPixel, VDimension>::TransformType *
SliceSeriesSpecialCoordinatesImage<TSliceImage, TTransform, TPixel, VDimension>::GetSliceTransform(
  IndexValueType sliceIndex) const
{
  const RegionType &   largestRegion = this->GetLargestPossibleRegion();
  const IndexValueType largestIndex = largestRegion.GetIndex(ImageDimension - 1);
  const SizeValueType  transformsIndex = static_cast<SizeValueType>(sliceIndex - largestIndex);
  if (this->m_SliceTransforms->Size() > transformsIndex)
  {
    return this->m_SliceTransforms->GetElement(transformsIndex).GetPointer();
  }
  return nullptr;
}


template <typename TSliceImage, typename TTransform, typename TPixel, unsigned int VDimension>
const typename SliceSeriesSpecialCoordinatesImage<TSliceImage, TTransform, TPixel, VDimension>::TransformType *
SliceSeriesSpecialCoordinatesImage<TSliceImage, TTransform, TPixel, VDimension>::GetSliceInverseTransform(
  IndexValueType sliceIndex) const
{
  const RegionType &   largestRegion = this->GetLargestPossibleRegion();
  const IndexValueType largestIndex = largestRegion.GetIndex(ImageDimension - 1);
  const SizeValueType  transformsIndex = static_cast<SizeValueType>(sliceIndex - largestIndex);
  if (this->m_SliceInverseTransforms->Size() > transformsIndex)
  {
    return this->m_SliceInverseTransforms->GetElement(transformsIndex).GetPointer();
  }
  return nullptr;
}


template <typename TSliceImage, typename TTransform, typename TPixel, unsigned int VDimension>
void
SliceSeriesSpecialCoordinatesImage<TSliceImage, TTransform, TPixel, VDimension>::PrintSelf(std::ostream & os,
                                                                                           Indent         indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "SliceImage = " << m_SliceImage << std::endl;
  os << indent << "SliceTransforms = " << m_SliceTransforms << std::endl;
}


template <typename TSliceImage, typename TTransform, typename TPixel, unsigned int VDimension>
void
SliceSeriesSpecialCoordinatesImage<TSliceImage, TTransform, TPixel, VDimension>::Graft(const DataObject * data)
{
  // call the superclass' implementation
  Superclass::Graft(data);

  if (data)
  {
    // Attempt to cast data to a SliceSeriesSpecialCoordinatesImage
    const Self * const imgData = dynamic_cast<const Self *>(data);

    if (imgData)
    {
      // Now copy anything remaining that is needed
      this->SetPixelContainer(const_cast<PixelContainer *>(imgData->GetPixelContainer()));
    }
    else
    {
      // pointer could not be cast back down
      itkExceptionMacro(<< "itk::Image::Graft() cannot cast " << typeid(data).name() << " to "
                        << typeid(const Self *).name());
    }
  }
}

} // end namespace itk

#endif
