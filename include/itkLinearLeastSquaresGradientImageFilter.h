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
#ifndef itkLinearLeastSquaresGradientImageFilter_h
#define itkLinearLeastSquaresGradientImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkCovariantVector.h"
#include "itkNeighborhoodAlgorithm.h"

namespace itk
{

/** \class LinearLeastSquaresGradientImageFilter
 *
 * \brief Calculate the gradient using a linear least squares fit to values in
 * the neighborhood slice.
 *
 * The size of the neighborhood is determined by the Radius.
 *
 * \warning this algorithm is adapted such that it is sensitive to changes in
 * sign of the underlying function.  If there are a majority of consecutive
 * values with the same sign, then only those values are used in the linear fit.
 *
 * \sa GradientImageFilter
 *
 * \ingroup GradientFilters
 * \ingroup Ultrasound
 *
 */
template <typename TInputImage, typename TOperatorValueType = double, typename TOutputValueType = double>
class ITK_TEMPLATE_EXPORT LinearLeastSquaresGradientImageFilter
  : public ImageToImageFilter<
      TInputImage,
      Image<CovariantVector<TOutputValueType, TInputImage::ImageDimension>, TInputImage::ImageDimension>>
{
public:
  /** Extract dimension from input image. */
  itkStaticConstMacro(ImageDimension, unsigned int, TInputImage::ImageDimension);

  /** Standard class type alias. */
  using Self = LinearLeastSquaresGradientImageFilter;

  /** Convenient type alias for simplifying declarations. */
  using InputImageType = TInputImage;
  using InputImagePointer = typename InputImageType::Pointer;
  using OutputImageType = Image<CovariantVector<TOutputValueType, itkGetStaticConstMacro(OutputImageDimension)>,
                                itkGetStaticConstMacro(OutputImageDimension)>;
  using OutputImagePointer = typename OutputImageType::Pointer;

  /** Standard class type alias. */
  using Superclass = ImageToImageFilter<InputImageType, OutputImageType>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(LinearLeastSquaresGradientImageFilter, ImageToImageFilter);

  /** Image type alias support. */
  using InputPixelType = typename InputImageType::PixelType;
  using OperatorValueType = TOperatorValueType;
  using OutputValueType = TOutputValueType;
  using OutputPixelType = CovariantVector<OutputValueType, itkGetStaticConstMacro(OutputImageDimension)>;
  using OutputImageRegionType = typename OutputImageType::RegionType;
  using RadiusType = typename InputImageType::SizeType;

  /** Set/Get whether or not the filter will use the spacing of the input
      image in its calculations */
  itkSetMacro(UseImageSpacing, bool);
  itkGetConstMacro(UseImageSpacing, bool);
  itkBooleanMacro(UseImageSpacing);

  /** Set/Get the radius of the linear least squares fit. */
  virtual void
  SetRadius(const RadiusType & radius)
  {
    m_Radius = radius;
    this->Modified();
  }
  /** Set the radius to the given value in all directions. */
  virtual void
  SetRadius(const unsigned int rad)
  {
    RadiusType radius;
    radius.Fill(rad);
    this->SetRadius(radius);
  }
  itkGetConstReferenceMacro(Radius, RadiusType);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputConvertibleToOutputCheck, (Concept::Convertible<InputPixelType, OutputValueType>));
  itkConceptMacro(OutputHasNumericTraitsCheck, (Concept::HasNumericTraits<OutputValueType>));
  /** End concept checking */
#endif

protected:
  LinearLeastSquaresGradientImageFilter();
  ~LinearLeastSquaresGradientImageFilter() override = default;

  using NeighborhoodIteratorType = ConstNeighborhoodIterator<InputImageType>;

  inline TOutputValueType
  GetDerivative(const std::slice & s, const NeighborhoodIteratorType & nit);

  /** GradientImageFilter needs a larger input requested region than
   * the output requested region.  As such, GradientImageFilter needs
   * to provide an implementation for GenerateInputRequestedRegion()
   * in order to inform the pipeline execution model.
   *
   * \sa ImageToImageFilter::GenerateInputRequestedRegion() */
  void
  GenerateInputRequestedRegion() override;


  /** GradientImageFilter can be implemented as a multithreaded filter.
   * Therefore, this implementation provides a ThreadedGenerateData()
   * routine which is called for each processing thread. The output
   * image data is allocated automatically by the superclass prior to
   * calling ThreadedGenerateData().  ThreadedGenerateData can only
   * write to the portion of the output image specified by the
   * parameter "outputRegionForThread"
   *
   * \sa ImageToImageFilter::ThreadedGenerateData(),
   *     ImageToImageFilter::GenerateData() */
  void
  DynamicThreadedGenerateData(const OutputImageRegionType & outputRegionForThread) override;

  /** The UseImageDirection flag determines whether image derivatives are
   * computed with respect to the image grid or with respect to the physical
   * space. When this flag is ON the derivatives are computed with respect to
   * the coodinate system of physical space. The difference is whether we take
   * into account the image Direction or not. The flag ON will take into
   * account the image direction and will result in an extra matrix
   * multiplication compared to the amount of computation performed when the
   * flag is OFF.
   * The default value of this flag is On.
   */
  itkSetMacro(UseImageDirection, bool);
  itkGetConstMacro(UseImageDirection, bool);
  itkBooleanMacro(UseImageDirection);

private:
  LinearLeastSquaresGradientImageFilter(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented

  bool       m_UseImageSpacing;
  RadiusType m_Radius;
  // flag to take or not the image direction into account
  // when computing the derivatives.
  bool m_UseImageDirection;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkLinearLeastSquaresGradientImageFilter.hxx"
#endif

#endif
