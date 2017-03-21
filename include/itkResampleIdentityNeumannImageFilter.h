/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkResampleIdentityNeumannImageFilter.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkResampleIdentityNeumannImageFilter_h
#define __itkResampleIdentityNeumannImageFilter_h

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
 */
template <class TInputImage, class TOutputImage,
          class TInterpolatorPrecisionType = double>
class ITK_EXPORT ResampleIdentityNeumannImageFilter :
  public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef ResampleIdentityNeumannImageFilter            Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  typedef TInputImage                           InputImageType;
  typedef TOutputImage                          OutputImageType;
  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;
  typedef typename OutputImageType::Pointer     OutputImagePointer;
  typedef typename InputImageType::RegionType   InputImageRegionType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ResampleIdentityNeumannImageFilter, ImageToImageFilter);

  /** Number of dimensions. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension);
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Interpolator typedef. */
  typedef InterpolateImageFunction<InputImageType, TInterpolatorPrecisionType> InterpolatorType;
  typedef typename InterpolatorType::Pointer  InterpolatorPointerType;

  /** Image size typedef. */
  typedef Size<itkGetStaticConstMacro(ImageDimension)> SizeType;

  /** Image index typedef. */
  typedef typename TOutputImage::IndexType IndexType;

  /** Image point typedef. */
  typedef typename InterpolatorType::PointType    PointType;
  //typedef typename TOutputImage::PointType    PointType;

  /** Image pixel value typedef. */
  typedef typename TOutputImage::PixelType   PixelType;
  typedef typename TInputImage::PixelType    InputPixelType;
  
  /** Typedef to describe the output image region type. */
  typedef typename TOutputImage::RegionType OutputImageRegionType;

  /** Image spacing,origin and direction typedef */
  typedef typename TOutputImage::SpacingType   SpacingType;
  typedef typename TOutputImage::PointType     OriginPointType;
  typedef typename TOutputImage::DirectionType DirectionType;

  /** base type for images of the current ImageDimension */
  typedef ImageBase<itkGetStaticConstMacro(ImageDimension)> ImageBaseType;

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

  /** ResampleIdentityNeumannImageFilter produces an image which is a different size
   * than its input.  As such, it needs to provide an implementation
   * for GenerateOutputInformation() in order to inform the pipeline
   * execution model.  The original documentation of this method is
   * below. \sa ProcessObject::GenerateOutputInformaton() */
  virtual void GenerateOutputInformation( void );

  /** ResampleIdentityNeumannImageFilter needs a different input requested region than
   * the output requested region.  As such, ResampleIdentityNeumannImageFilter needs
   * to provide an implementation for GenerateInputRequestedRegion()
   * in order to inform the pipeline execution model.
   * \sa ProcessObject::GenerateInputRequestedRegion() */
  virtual void GenerateInputRequestedRegion( void );

  /** This method is used to set the state of the filter before 
   * multi-threading. */
  virtual void BeforeThreadedGenerateData( void );

  /** This method is used to set the state of the filter after 
   * multi-threading. */
  virtual void AfterThreadedGenerateData( void );

  /** Method Compute the Modified Time based on changed to the components. */
  unsigned long GetMTime( void ) const;

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(OutputHasNumericTraitsCheck,
                  (Concept::HasNumericTraits<PixelType>));
  /** End concept checking */
#endif

protected:
  ResampleIdentityNeumannImageFilter( void );
  ~ResampleIdentityNeumannImageFilter( void ) {};

  void PrintSelf( std::ostream& os, Indent indent ) const;

  /** ResampleIdentityNeumannImageFilter can be implemented as a multithreaded filter.  Therefore,
   * this implementation provides a ThreadedGenerateData() routine which
   * is called for each processing thread. The output image data is allocated
   * automatically by the superclass prior to calling ThreadedGenerateData().
   * ThreadedGenerateData can only write to the portion of the output image
   * specified by the parameter "outputRegionForThread"
   * \sa ImageToImageFilter::ThreadedGenerateData(),
   *     ImageToImageFilter::GenerateData() */
  void ThreadedGenerateData( const OutputImageRegionType& outputRegionForThread,
                             int threadId );

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
#include "itkResampleIdentityNeumannImageFilter.txx"
#endif

#endif