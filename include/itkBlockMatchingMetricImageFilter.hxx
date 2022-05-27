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
#ifndef itkBlockMatchingMetricImageFilter_hxx
#define itkBlockMatchingMetricImageFilter_hxx


namespace itk
{
namespace BlockMatching
{

template <typename TFixedImage, typename TMovingImage, typename TMetricImage>
MetricImageFilter<TFixedImage, TMovingImage, TMetricImage>::MetricImageFilter()
  : m_FixedImageRegionDefined(false)
  , m_MovingImageRegionDefined(false)
{}


template <typename TFixedImage, typename TMovingImage, typename TMetricImage>
void
MetricImageFilter<TFixedImage, TMovingImage, TMetricImage>::SetFixedImage(FixedImageType * fixedImage)
{
  this->SetInput(0, fixedImage);
}


template <typename TFixedImage, typename TMovingImage, typename TMetricImage>
void
MetricImageFilter<TFixedImage, TMovingImage, TMetricImage>::SetMovingImage(MovingImageType * movingImage)
{
  this->SetInput(1, movingImage);
}


template <typename TFixedImage, typename TMovingImage, typename TMetricImage>
void
MetricImageFilter<TFixedImage, TMovingImage, TMetricImage>::SetFixedImageRegion(const FixedImageRegionType & region)
{
  FixedImageType * fixedImage = const_cast<TFixedImage *>(this->GetInput(0));
  if (!fixedImage)
  {
    itkExceptionMacro(<< "The FixedImage must be set before specifying the fixed image region.");
  }
  fixedImage->UpdateOutputInformation();
  m_FixedImageRegion = region;
  if (!m_FixedImageRegion.Crop(fixedImage->GetLargestPossibleRegion()))
  {
    itkExceptionMacro(<< "Requested block is outside of the fixed image."
                      << " block: " << region << " fixed image: " << fixedImage->GetLargestPossibleRegion());
  }
  typename FixedImageRegionType::SizeType fixedSize = m_FixedImageRegion.GetSize();
  for (unsigned int i = 0; i < ImageDimension; ++i)
  {
    // The radius may have been truncated if the fixed region was outside the
    // fixed image's LargestPossibleRegion.
    if (fixedSize[i] % 2 == 0)
      fixedSize[i]--;
    m_FixedRadius[i] = (fixedSize[i] - 1) / 2;
  }
  m_FixedImageRegion.SetSize(fixedSize);
  m_FixedImageRegionDefined = true;

  MovingImageType * moving = const_cast<TMovingImage *>(this->GetInput(1));
  if (!moving)
  {
    itkExceptionMacro(<< "The MovingImage must be set before specifying the fixed image region.");
  }
  moving->UpdateOutputInformation();
  m_MovingRadius = m_FixedRadius;
  typename FixedImageType::SpacingType  fixedSpacing = fixedImage->GetSpacing();
  typename MovingImageType::SpacingType movingSpacing = moving->GetSpacing();
  if (!(fixedSpacing == movingSpacing))
  {
    for (unsigned int i = 0; i < ImageDimension; ++i)
    {
      m_MovingRadius[i] =
        Math::Ceil<typename RadiusType::SizeValueType>(fixedSpacing[i] * m_FixedRadius[i] / movingSpacing[i]);
    }
  }
  this->Modified();
}


template <typename TFixedImage, typename TMovingImage, typename TMetricImage>
void
MetricImageFilter<TFixedImage, TMovingImage, TMetricImage>::SetMovingImageRegion(const MovingImageRegionType & region)
{
  m_MovingImageRegion = region;
  m_MovingImageRegionDefined = true;
  this->Modified();
}


template <typename TFixedImage, typename TMovingImage, typename TMetricImage>
void
MetricImageFilter<TFixedImage, TMovingImage, TMetricImage>::GenerateOutputInformation()
{
  const MovingImageType * moving = this->GetInput(1);

  MetricImageType * output = this->GetOutput();

  if (!m_MovingImageRegionDefined)
  {
    itkExceptionMacro(<< "MovingImageRegion has not been set");
  }

  MetricImageRegionType                     metricRegion;
  typename MetricImageRegionType::IndexType metricIndex;
  metricIndex.Fill(0);
  metricRegion.SetIndex(metricIndex);
  // the Default is to to use the moving image size and spacing.
  metricRegion.SetSize(m_MovingImageRegion.GetSize());
  output->SetLargestPossibleRegion(metricRegion);
  output->SetSpacing(moving->GetSpacing());

  typename MetricImageType::IndexType metricStart(m_MovingImageRegion.GetIndex());

  typename MetricImageType::PointType origin;
  moving->TransformIndexToPhysicalPoint(metricStart, origin);
  output->SetOrigin(origin);

  // The metric image direction is the same as the fixed image direction.
  output->SetDirection(moving->GetDirection());
}


template <typename TFixedImage, typename TMovingImage, typename TMetricImage>
void
MetricImageFilter<TFixedImage, TMovingImage, TMetricImage>::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();

  // const cast so we can set the requested region.
  FixedImageType * fixedImage = const_cast<TFixedImage *>(this->GetInput(0));
  if (!fixedImage)
  {
    return;
  }

  MovingImageType * moving = const_cast<TMovingImage *>(this->GetInput(1));
  if (!moving)
  {
    return;
  }

  if (!m_FixedImageRegionDefined)
  {
    itkExceptionMacro(<< "FixedImageRegion has not been set");
  }

  if (!m_MovingImageRegionDefined)
  {
    itkExceptionMacro(<< "MovingImageRegion has not been set");
  }

  fixedImage->SetRequestedRegion(m_FixedImageRegion);
  MovingImageRegionType movingImageRequestedRegion = m_MovingImageRegion;
  movingImageRequestedRegion.PadByRadius(m_MovingRadius);
  // make sure the requested region is within the largest possible.
  if (movingImageRequestedRegion.Crop(moving->GetLargestPossibleRegion()))
  {
    moving->SetRequestedRegion(movingImageRequestedRegion);
    return;
  }
  else
  {
    // store what we tried( prior to try to crop )
    moving->SetRequestedRegion(movingImageRequestedRegion);

    itkExceptionMacro(<< "Moving image requested region is at least partially outside the LargestPossibleRegion.");
  }
}


template <typename TFixedImage, typename TMovingImage, typename TMetricImage>
void
MetricImageFilter<TFixedImage, TMovingImage, TMetricImage>::EnlargeOutputRequestedRegion(DataObject * data)
{
  Superclass::EnlargeOutputRequestedRegion(data);
  data->SetRequestedRegionToLargestPossibleRegion();
}

} // end namespace BlockMatching
} // end namespace itk

#endif
