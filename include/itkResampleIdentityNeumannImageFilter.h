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
#ifndef itkResampleIdentityNeumannImageFilter_h
#define itkResampleIdentityNeumannImageFilter_h

#include "itkFixedArray.h"
#include "itkTransform.h"
#include "itkImageFunction.h"
#include "itkImageRegionIterator.h"
#include "itkImageToImageFilter.h"
#include "itkInterpolateImageFunction.h"
#include "itkSize.h"

namespace itk
{

/** \class ResampleIdentityNeumannImageFilter
 * \brief Resample an image via an identity transform and use Neumann boundary
 * conditions.
 *
 * Instead of using a default pixel value when a requested pixel value is
 * outside the image boundary, like
 * ResampleImageFilter, the filter uses the value at the nearest pixel.
 *
 * Note that the choice of interpolator function can be important.
 * This function is set via SetInterpolator().  The default is
 * itk::LinearInterpolateImageFunction<InputImageType, TInterpolatorPrecisionType>, which
 * is reasonable for ordinary medical images.  However, some synthetic
 * images have pixels drawn from a finite prescribed set.  An example
 * would be a mask indicating the segmentation of a brain into a small
 * number of tissue types.  For such an image, one does not want to
 * interpolate between different pixel values, and so
 * itk::NearestNeighborInterpolateImageFunction< InputImageType,
 * TCoordRep > would be a better choice.
 *
 * Output information (spacing, size and direction) for the output
 * image should be set. This information has the normal defaults of
 * unit spacing, zero origin and identity direction. Optionally, the
 * output information can be obtained from a reference image. If the
 * reference image is provided and UseReferenceImage is On, then the
 * spacing, origin and direction of the reference image will be used.
 *
 * Since this filter produces an image which is a different size than
 * its input, it needs to override several of the methods defined
 * in ProcessObject in order to properly manage the pipeline execution model.
 * In particular, this filter overrides
 * ProcessObject::GenerateInputRequestedRegion() and
 * ProcessObject::GenerateOutputInformation().
 *
 * This filter is implemented as a multithreaded filter.  It provides a
 * ThreadedGenerateData() method for its implementation.
 *
 * \sa ResampleImageFilter
 *
 * \ingroup Ultrasound
 */
template <class TInputImage, class TOutputImage,
          class TInterpolatorPrecisionType = double>
class ITK_TEMPLATE_EXPORT ResampleIdentityNeumannImageFilter :
  public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class type alias. */
  using Self = ResampleIdentityNeumannImageFilter;
  using Superclass = ImageToImageFilter<TInputImage, TOutputImage>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;
  using InputImagePointer = typename InputImageType::Pointer;
  using InputImageConstPointer = typename InputImageType::ConstPointer;
  using OutputImagePointer = typename OutputImageType::Pointer;
  using InputImageRegionType = typename InputImageType::RegionType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ResampleIdentityNeumannImageFilter, ImageToImageFilter);

  /** Number of dimensions. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension);
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Interpolator type alias. */
  using InterpolatorType = InterpolateImageFunction<InputImageType, TInterpolatorPrecisionType>;
  using InterpolatorPointerType = typename InterpolatorType::Pointer;

  /** Image size type alias. */
  using SizeType = Size<itkGetStaticConstMacro(ImageDimension)>;

  /** Image index type alias. */
  using IndexType = typename TOutputImage::IndexType;

  /** Image point type alias. */
  using PointType = typename InterpolatorType::PointType;
  //using PointType = typename TOutputImage::PointType;

  /** Image pixel value type alias. */
  using PixelType = typename TOutputImage::PixelType;
  using InputPixelType = typename TInputImage::PixelType;

  /** Typedef to describe the output image region type. */
  using OutputImageRegionType = typename TOutputImage::RegionType;

  /** Image spacing,origin and direction type alias */
  using SpacingType = typename TOutputImage::SpacingType;
  using OriginPointType = typename TOutputImage::PointType;
  using DirectionType = typename TOutputImage::DirectionType;

  /** base type for images of the current ImageDimension */
  using ImageBaseType = ImageBase<itkGetStaticConstMacro(ImageDimension)>;

  /** Set the interpolator function.  The default is
   * itk::LinearInterpolateImageFunction<InputImageType, TInterpolatorPrecisionType>. Some
   * other options are itk::NearestNeighborInterpolateImageFunction
   * (useful for binary masks and other images with a small number of
   * possible pixel values), and itk::BSplineInterpolateImageFunction
   * (which provides a higher order of interpolation).  */
  itkSetObjectMacro( Interpolator, InterpolatorType );

  /** Get a pointer to the interpolator function. */
  itkGetConstObjectMacro( Interpolator, InterpolatorType );

  /** Set the size of the output image. */
  itkSetMacro( Size, SizeType );

  /** Get the size of the output image. */
  itkGetConstReferenceMacro( Size, SizeType );

  /** Set the output image spacing. */
  itkSetMacro( OutputSpacing, SpacingType );
  virtual void SetOutputSpacing( const double* values );

  /** Get the output image spacing. */
  itkGetConstReferenceMacro( OutputSpacing, SpacingType );

  /** Set the output image origin. */
  itkSetMacro( OutputOrigin, OriginPointType );
  virtual void SetOutputOrigin( const double* values);

  /** Get the output image origin. */
  itkGetConstReferenceMacro( OutputOrigin, OriginPointType );

  /** Set the output direciton cosine matrix. */
  itkSetMacro( OutputDirection, DirectionType );
  itkGetConstReferenceMacro( OutputDirection, DirectionType );

  /** Helper method to set the output parameters based on this image */
  void SetOutputParametersFromImage ( const ImageBaseType * image );

  /** Set the start index of the output largest possible region.
   * The default is an index of all zeros. */
  itkSetMacro( OutputStartIndex, IndexType );

  /** Get the start index of the output largest possible region. */
  itkGetConstReferenceMacro( OutputStartIndex, IndexType );

  /** Copy the output information from another Image.  By default,
   *  the information is specified with the SetOutputSpacing, Origin,
   *  and Direction methods. UseReferenceImage must be On and a
   *  Reference image must be present to override the defaul behavior.
   *  NOTE: This function seems redundant with the
   *  SetOutputParametersFromImage( image ) function */
  void SetReferenceImage ( const TOutputImage *image );
  const TOutputImage * GetReferenceImage( void ) const;

  itkSetMacro( UseReferenceImage, bool );
  itkBooleanMacro( UseReferenceImage );
  itkGetConstMacro( UseReferenceImage, bool );

protected:
  ResampleIdentityNeumannImageFilter( void );
  ~ResampleIdentityNeumannImageFilter( void ) {};

  /** ResampleIdentityNeumannImageFilter produces an image which is a different size
   * than its input.  As such, it needs to provide an implementation
   * for GenerateOutputInformation() in order to inform the pipeline
   * execution model.  The original documentation of this method is
   * below. \sa ProcessObject::GenerateOutputInformaton() */
  void GenerateOutputInformation() override;

  /** ResampleIdentityNeumannImageFilter needs a different input requested region than
   * the output requested region.  As such, ResampleIdentityNeumannImageFilter needs
   * to provide an implementation for GenerateInputRequestedRegion()
   * in order to inform the pipeline execution model.
   * \sa ProcessObject::GenerateInputRequestedRegion() */
  void GenerateInputRequestedRegion() override;

  /** This method is used to set the state of the filter before
   * multi-threading. */
  void BeforeThreadedGenerateData() override;

  /** This method is used to set the state of the filter after
   * multi-threading. */
  void AfterThreadedGenerateData() override;

  /** ResampleIdentityNeumannImageFilter can be implemented as a multithreaded filter.  Therefore,
   * this implementation provides a ThreadedGenerateData() routine which
   * is called for each processing thread. The output image data is allocated
   * automatically by the superclass prior to calling ThreadedGenerateData().
   * ThreadedGenerateData can only write to the portion of the output image
   * specified by the parameter "outputRegionForThread"
   * \sa ImageToImageFilter::ThreadedGenerateData(),
   *     ImageToImageFilter::GenerateData() */
  void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread,
                             ThreadIdType threadId ) override;

  /** Method Compute the Modified Time based on changed to the components. */
  unsigned long GetMTime( void ) const override;

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(OutputHasNumericTraitsCheck,
                  (Concept::HasNumericTraits<PixelType>));
  /** End concept checking */
#endif


  void PrintSelf( std::ostream& os, Indent indent ) const override;

private:
  ResampleIdentityNeumannImageFilter( const Self& ); //purposely not implemented
  void operator=( const Self& ); //purposely not implemented

  SizeType                m_Size;              // Size of the output image
  InterpolatorPointerType m_Interpolator;      // Image function for
                                               // interpolation
  SpacingType             m_OutputSpacing;     // output image spacing
  OriginPointType         m_OutputOrigin;      // output image origin
  DirectionType           m_OutputDirection;   // output image direction cosines
  IndexType               m_OutputStartIndex;  // output image start index
  bool                    m_UseReferenceImage;
};


} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkResampleIdentityNeumannImageFilter.hxx"
#endif

#endif
