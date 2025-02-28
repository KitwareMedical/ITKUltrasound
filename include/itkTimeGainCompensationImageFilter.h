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
#ifndef itkTimeGainCompensationImageFilter_h
#define itkTimeGainCompensationImageFilter_h

#include "itkImageToImageFilter.h"

#include "itkArray2D.h"

namespace itk
{

/**
 * \class TimeGainCompensationImageFilter
 * \brief Applies a linear piecewise time gain compensation.
 *
 * This filter applies a linear piecewise gain with depth.  The depth
 * direction is assumed to be the first direction (0th direction).
 *
 * \ingroup Ultrasound
 * */
template <typename TInputImage, typename TOutputImage = TInputImage>
class ITK_TEMPLATE_EXPORT TimeGainCompensationImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(TimeGainCompensationImageFilter);

  /** Standard class type alias. */
  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;

  using Self = TimeGainCompensationImageFilter;
  using Superclass = ImageToImageFilter<InputImageType, OutputImageType>;

  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  itkOverrideGetNameOfClassMacro(TimeGainCompensationImageFilter);
  itkNewMacro(Self);

  using GainType = Array2D<double>;

  /** Set/Get the gain.  The first column specifies the depth. The second
   * column specifies the gain. */
  itkSetMacro(Gain, GainType);
  itkGetConstReferenceMacro(Gain, GainType);

protected:
  using OutputImageRegionType = typename OutputImageType::RegionType;

  TimeGainCompensationImageFilter();
  ~TimeGainCompensationImageFilter() override = default;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  void
  BeforeThreadedGenerateData() override;
  void
  DynamicThreadedGenerateData(const OutputImageRegionType & outputRegionForThread) override;

private:
  GainType m_Gain;
};

} // namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkTimeGainCompensationImageFilter.hxx"
#endif

#endif // itkTimeGainCompensationImageFilter_h
