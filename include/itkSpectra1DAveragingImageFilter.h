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
#ifndef itkSpectra1DAveragingImageFilter_h
#define itkSpectra1DAveragingImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkDefaultConvertPixelTraits.h"
#include "itkVectorImage.h"

#include <type_traits>

namespace itk
{

/** \class Spectra1DAveragingImageFilter
 * \brief Average multiple spectra in input images into a single response line.
 *
 *
 * \ingroup Ultrasound
 */
template <typename TInputImage, typename TOutputImage = TInputImage>
class ITK_TEMPLATE_EXPORT Spectra1DAveragingImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(Spectra1DAveragingImageFilter);

  static constexpr unsigned int InputImageDimension = TInputImage::ImageDimension;
  static constexpr unsigned int OutputImageDimension = TOutputImage::ImageDimension;

  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;

  using InputPixelType = typename InputImageType::PixelType;
  using OutputPixelType = typename OutputImageType::PixelType;

  using InputScalarType = typename DefaultConvertPixelTraits<typename InputImageType::PixelType>::ComponentType;
  using OutputScalarType = typename DefaultConvertPixelTraits<typename OutputImageType::PixelType>::ComponentType;

  /** Standard class type alias. */
  using Self = Spectra1DAveragingImageFilter;
  using Superclass = ImageToImageFilter<InputImageType, OutputImageType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  itkTypeMacro(Spectra1DAveragingImageFilter, ImageToImageFilter);
  itkNewMacro(Self);

#ifdef ITK_USE_CONCEPT_CHECKING
  // Begin concept checking
  itkConceptMacro(FloatConvertibleToInputScalarTypeCheck, (Concept::Convertible<float, InputScalarType>));
  itkConceptMacro(FloatConvertibleToOutputScalarTypeCheck, (Concept::Convertible<float, OutputScalarType>));
  itkConceptMacro(InputScalarTypeConvertibleToFloatCheck, (Concept::Convertible<InputScalarType, float>));
  itkConceptMacro(OutputScalarTypeConvertibleToFloatCheck, (Concept::Convertible<OutputScalarType, float>));
  // End concept checking
#endif

protected:
  Spectra1DAveragingImageFilter() = default;
  ~Spectra1DAveragingImageFilter() override = default;

  void
  GenerateOutputInformation() override;

  void
  GenerateData() override;

  void
  VerifyInputInformation() const override
  {}


  using InputVectorImage = VectorImage<InputScalarType, InputImageDimension>;
  using OutputVectorImage = VectorImage<OutputScalarType, OutputImageDimension>;

  // sets the vector length if output is a VectorImage
  template <typename TIn, typename TOut>
  std::enable_if_t<!std::is_same<OutputImageType, VectorImage<OutputScalarType, TOut::ImageDimension>>::value>
  PrepareOutput(TIn *, TOut *)
  {
    // nothing to do if output is not a VectorImage
  }
  template <typename TIn, typename TOut>
  std::enable_if_t<!std::is_same<InputImageType, VectorImage<InputScalarType, TIn::ImageDimension>>::value &&
                   std::is_same<OutputImageType, VectorImage<OutputScalarType, TOut::ImageDimension>>::value>
  PrepareOutput(TIn * in, TOut * out)
  {
    this->GetOutput()->SetVectorLength(TIn::PixelType::Dimension);
  }
  template <typename TIn, typename TOut>
  std::enable_if_t<std::is_same<InputImageType, VectorImage<InputScalarType, TIn::ImageDimension>>::value &&
                   std::is_same<OutputImageType, VectorImage<OutputScalarType, TOut::ImageDimension>>::value>
  PrepareOutput(TIn * in, TOut * out)
  {
    this->GetOutput()->SetVectorLength(in->GetVectorLength());
  }

  using IndexType = typename InputImageType::IndexType;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkSpectra1DAveragingImageFilter.hxx"
#endif

#endif // itkSpectra1DAveragingImageFilter_h
