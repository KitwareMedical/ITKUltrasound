/*=========================================================================
 *
 *  Copyright Insight Software Consortium
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
#ifndef itkSliceSeriesSpecialCoordinatesImage_h
#define itkSliceSeriesSpecialCoordinatesImage_h

#include "itkSpecialCoordinatesImage.h"
#include "itkPoint.h"
#include "itkNeighborhoodAccessorFunctor.h"
#include "itkVectorContainer.h"
#include "itkMath.h"

namespace itk
{
/** \class SliceSeriesSpecialCoordinatesImage
 *
 *  \brief An image composed of a series of adjacent slices that are not
 *  necessarily uniformly spaced.
 *
 * This is an itk::SpecialCoordinatesImage comprised of a series of N-1 dimension image
 * slices. It is assumed that these slices are adjacent to each other without
 * overlap or intersection.
 *
 * The physical location of a pixel Index is determined by first calling
 * TransformIndexToPhysicalPoint for dimensions 0 to N-1, then by running the
 * result through a Transform. The metadata for the image that defines its
 * spatial domain are the SliceImage and the SliceTransform's. It is possible
 * for the Slice image to be a specialized type like the
 * CurvilinearArraySpecialCoordinatesImage.
 *
 * \sa SpecialCoordinatesImage
 * \sa CurvilinearArraySpecialCoordinatesImage
 *
 * \ingroup Ultrasound
 *
 * \ingroup ImageObjects
 */
template< typename TSliceImage, typename TTransform, typename TPixel=typename TSliceImage::PixelType, unsigned int VDimension=TSliceImage::ImageDimension + 1 >
class ITK_TEMPLATE_EXPORT SliceSeriesSpecialCoordinatesImage:
  public SpecialCoordinatesImage< TPixel, VDimension >
{
public:
  /** Standard class typedefs */
  typedef SliceSeriesSpecialCoordinatesImage            Self;
  typedef SpecialCoordinatesImage< TPixel, VDimension > Superclass;
  typedef SmartPointer< Self >                          Pointer;
  typedef SmartPointer< const Self >                    ConstPointer;
  typedef WeakPointer< const Self >                     ConstWeakPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SliceSeriesSpecialCoordinatesImage, SpecialCoordinatesImage);

  /** Pixel typedef support. Used to declare pixel type in filters
   * or other operations. */
  typedef TPixel PixelType;

  /** Typedef alias for PixelType */
  typedef TPixel ValueType;

  /** Typedef of the image slice. This is usually one less dimension of the
   * image.*/
  typedef TSliceImage SliceImageType;

  /** Typedef of the transform used to transform each slice. This could also
   * be a base class transform type. */
  typedef TTransform TransformType;

  /** Typedef of the array of transform used to store the per-slice
   * transforms. */
  typedef VectorContainer< IdentifierType, typename TransformType::Pointer > SliceTransformsType;

  /** Internal Pixel representation. Used to maintain a uniform API
   * with Image Adaptors and allow to keep a particular internal
   * representation of data while showing a different external
   * representation. */
  typedef TPixel InternalPixelType;

  typedef typename Superclass::IOPixelType IOPixelType;

  /** Accessor type that convert data between internal and external
   *  representations.  */
  typedef DefaultPixelAccessor< PixelType > AccessorType;

  /** Accessor functor to choose between accessors: DefaultPixelAccessor for
   * the Image, and DefaultVectorPixelAccessor for the vector image. The
   * functor provides a generic API between the two accessors. */
  typedef DefaultPixelAccessorFunctor< Self > AccessorFunctorType;

  /** Typedef for the functor used to access a neighborhood of pixel
   * pointers. */
  typedef NeighborhoodAccessorFunctor< Self > NeighborhoodAccessorFunctorType;

  /** Dimension of the image.  This constant is used by functions that are
   * templated over image type (as opposed to being templated over pixel type
   * and dimension) when they need compile time access to the dimension of
   * the image. */
  itkStaticConstMacro(ImageDimension, unsigned int, VDimension);

  /** Index typedef support. An index is used to access pixel values. */
  typedef typename Superclass::IndexType      IndexType;
  typedef typename Superclass::IndexValueType IndexValueType;

  /** Offset typedef support. An offset is used to access pixel values. */
  typedef typename Superclass::OffsetType    OffsetType;

  /** Size typedef support. A size is used to define region bounds. */
  typedef typename Superclass::SizeType      SizeType;
  typedef typename Superclass::SizeValueType SizeValueType;

  /** Container used to store pixels in the image. */
  typedef ImportImageContainer< SizeValueType, PixelType > PixelContainer;

  /** Region typedef support. A region is used to specify a subset of
   *  an image.
   */
  typedef typename Superclass::RegionType RegionType;

  /** Spacing typedef support.  Spacing holds the "fake" size of a
   *  pixel, making each pixel look like a 1 unit hyper-cube to filters
   *  that were designed for normal images and that therefore use
   *  m_Spacing.  The spacing is the geometric distance between image
   *  samples.
   */
  typedef typename Superclass::SpacingType SpacingType;

  /** Origin typedef support.  The origin is the "fake" geometric
   *  coordinates of the index (0,0).  Also for use w/ filters designed
   *  for normal images.
   */
  typedef typename Superclass::PointType PointType;

  /** A pointer to the pixel container. */
  typedef typename PixelContainer::Pointer      PixelContainerPointer;
  typedef typename PixelContainer::ConstPointer PixelContainerConstPointer;

  /** Set/Get the Slice Image. This image is used to find the relative pixel
   * location locations within an image. Its size must be the same size as the
   * the size in the corresponding dimenions of this image. The pixel values
   * in this image are not used. */
  itkSetObjectMacro(SliceImage, SliceImageType);
  itkGetConstObjectMacro(SliceImage, SliceImageType);

  void SetSliceTransform( IndexValueType sliceIndex, TransformType * transform );
  const TransformType * GetSliceTransform( IndexValueType sliceIndex ) const;

  /** Graft the data and information from one image to another. This
   * is a convenience method to setup a second image with all the meta
   * information of another image and use the same pixel
   * container. Note that this method is different than just using two
   * SmartPointers to the same image since separate DataObjects are
   * still maintained. This method is similar to
   * ImageSource::GraftOutput(). The implementation in ImageBase
   * simply calls CopyInformation() and copies the region ivars.
   * The implementation here refers to the superclass' implementation
   * and then copies over the pixel container. */
  virtual void Graft(const DataObject *data) override;

  /** \brief Get the continuous index from a physical point
   *
   * Returns true if the resulting index is within the image, false otherwise.
   * \sa Transform */
  template< typename TCoordRep, typename TIndexRep >
  bool TransformPhysicalPointToContinuousIndex(
    const Point< TCoordRep, VDimension > & point,
    ContinuousIndex< TIndexRep, VDimension > & index) const
  {
    const RegionType & region = this->GetLargestPossibleRegion();
    const unsigned int sliceDimensionIndex = ImageDimension - 1;
    IndexValueType lowerIndex = region.GetIndex( sliceDimensionIndex );
    IndexValueType upperIndex = lowerIndex + region.GetSize( sliceDimensionIndex ) - 1;
    IndexValueType nextIndex = lowerIndex;
    PointType lowerPoint;
    const TransformType * transform = this->GetSliceInverseTransform( lowerIndex );
    if( transform == ITK_NULLPTR )
      {
      itkExceptionMacro("Inverse slice transform not available for index: " << lowerIndex);
      }
    lowerPoint = transform->TransformPoint( point );
    int lowerSign = Math::sgn( lowerPoint[sliceDimensionIndex] );
    PointType upperPoint;
    transform = this->GetSliceInverseTransform( upperIndex );
    if( transform == ITK_NULLPTR )
      {
      itkExceptionMacro("Inverse slice transform not available for index: " << upperIndex);
      }
    upperPoint = transform->TransformPoint( point );
    int upperSign = Math::sgn( upperPoint[sliceDimensionIndex] );
    PointType nextPoint = lowerPoint;
    int nextSign = 0;
    if( lowerSign == 0 )
      {
      nextSign = 0;
      nextPoint = lowerPoint;
      }
    else if( upperSign == 0 )
      {
      nextSign = 0;
      nextPoint = upperPoint;
      }
    else
      {
      if( lowerSign == upperSign )
        {
        // outside the image
        return false;
        }

      // Binary search for the transforms that bounds the slice
      while( upperIndex - lowerIndex > 1 )
        {
        nextIndex = lowerIndex + ( upperIndex - lowerIndex ) / 2;
        transform = this->GetSliceInverseTransform( nextIndex );
        nextPoint = transform->TransformPoint( point );
        nextSign = Math::sgn( nextPoint[sliceDimensionIndex] );
        if( nextSign == 0 )
          {
          break;
          }
        else if( nextSign == lowerSign )
          {
          lowerIndex = nextIndex;
          lowerPoint = nextPoint;
          }
        else
          {
          upperIndex = nextIndex;
          upperPoint = nextPoint;
          }
        }
      }

    if( nextSign != 0 )
      {
      const TCoordRep fraction = - lowerPoint[sliceDimensionIndex] / ( upperPoint[sliceDimensionIndex] - lowerPoint[sliceDimensionIndex] );
      nextPoint[sliceDimensionIndex] = lowerIndex + fraction * (upperIndex - lowerIndex);
      for( unsigned int ii = 0; ii < SliceImageType::ImageDimension; ++ii )
        {
        nextPoint[ii] = lowerPoint[ii] + fraction * (upperPoint[ii] - lowerPoint[ii]);
        }
      }
    else
      {
      nextPoint[sliceDimensionIndex] = nextIndex;
      }

    typename SliceImageType::PointType slicePoint;
    for( unsigned int ii = 0; ii < SliceImageType::ImageDimension; ++ii )
      {
      slicePoint[ii] = nextPoint[ii];
      }
    ContinuousIndex< TIndexRep, SliceImageType::ImageDimension > sliceIndex;
    this->m_SliceImage->TransformPhysicalPointToContinuousIndex( slicePoint, sliceIndex );
    for( unsigned int ii = 0; ii < SliceImageType::ImageDimension; ++ii )
      {
      index[ii] = sliceIndex[ii];
      }
    index[sliceDimensionIndex] = nextPoint[sliceDimensionIndex];
    // Now, check to see if the index is within allowed bounds
    const bool isInside = region.IsInside( index );

    return isInside;
  }

  /** Get the index (discrete) from a physical point.
   * Floating point index results are truncated to integers.
   * Returns true if the resulting index is within the image, false otherwise
   * \sa Transform */
  template< typename TCoordRep >
  bool TransformPhysicalPointToIndex(
    const Point< TCoordRep, VDimension > & point,
    IndexType & index) const
  {
    const RegionType & region = this->GetLargestPossibleRegion();
    const unsigned int sliceDimensionIndex = ImageDimension - 1;
    IndexValueType lowerIndex = region.GetIndex( sliceDimensionIndex );
    IndexValueType upperIndex = lowerIndex + region.GetSize( sliceDimensionIndex ) - 1;
    IndexValueType nextIndex = lowerIndex;
    PointType lowerPoint;
    const TransformType * transform = this->GetSliceInverseTransform( lowerIndex );
    if( transform == ITK_NULLPTR )
      {
      itkExceptionMacro("Inverse slice transform not available for index: " << lowerIndex);
      }
    lowerPoint = transform->TransformPoint( point );
    int lowerSign = Math::sgn( lowerPoint[sliceDimensionIndex] );
    PointType upperPoint;
    transform = this->GetSliceInverseTransform( upperIndex );
    if( transform == ITK_NULLPTR )
      {
      itkExceptionMacro("Inverse slice transform not available for index: " << upperIndex);
      }
    upperPoint = transform->TransformPoint( point );
    int upperSign = Math::sgn( upperPoint[sliceDimensionIndex] );
    PointType nextPoint = lowerPoint;
    int nextSign = 0;
    if( lowerSign == 0 )
      {
      nextSign = 0;
      nextPoint = lowerPoint;
      }
    else if( upperSign == 0 )
      {
      nextSign = 0;
      nextPoint = upperPoint;
      }
    else
      {
      if( lowerSign == upperSign )
        {
        // outside the image
        return false;
        }

      // Binary search for the transforms that bounds the slice
      while( upperIndex - lowerIndex > 1 )
        {
        nextIndex = lowerIndex + ( upperIndex - lowerIndex ) / 2;
        transform = this->GetSliceInverseTransform( nextIndex );
        nextPoint = transform->TransformPoint( point );
        nextSign = Math::sgn( nextPoint[sliceDimensionIndex] );
        if( nextSign == 0 )
          {
          break;
          }
        else if( nextSign == lowerSign )
          {
          lowerIndex = nextIndex;
          lowerPoint = nextPoint;
          }
        else
          {
          upperIndex = nextIndex;
          upperPoint = nextPoint;
          }
        }
      }

    if( nextSign != 0 )
      {
      const TCoordRep fraction = - lowerPoint[sliceDimensionIndex] / ( upperPoint[sliceDimensionIndex] - lowerPoint[sliceDimensionIndex] );
      nextPoint[sliceDimensionIndex] = lowerIndex + fraction * (upperIndex - lowerIndex);
      for( unsigned int ii = 0; ii < SliceImageType::ImageDimension; ++ii )
        {
        nextPoint[ii] = lowerPoint[ii] + fraction * (upperPoint[ii] - lowerPoint[ii]);
        }
      }
    else
      {
      nextPoint[sliceDimensionIndex] = nextIndex;
      }

    typename SliceImageType::PointType slicePoint;
    for( unsigned int ii = 0; ii < SliceImageType::ImageDimension; ++ii )
      {
      slicePoint[ii] = nextPoint[ii];
      }
    typename SliceImageType::IndexType sliceIndex;
    this->m_SliceImage->TransformPhysicalPointToIndex( slicePoint, sliceIndex );
    for( unsigned int ii = 0; ii < SliceImageType::ImageDimension; ++ii )
      {
      index[ii] = sliceIndex[ii];
      }
    index[sliceDimensionIndex] = Math::RoundHalfIntegerUp< IndexValueType >( nextPoint[sliceDimensionIndex] );
    // Now, check to see if the index is within allowed bounds
    const bool isInside = region.IsInside( index );

    return isInside;
  }

  /** Get a physical point (in the space which
   * the origin and spacing information comes from)
   * from a continuous index (in the index space)
   * \sa Transform */
  template< typename TCoordRep, typename TIndexRep >
  void TransformContinuousIndexToPhysicalPoint(
    const ContinuousIndex< TIndexRep, VDimension > & index,
    Point< TCoordRep, VDimension > & point) const
  {
    point.Fill( 0.0 );
    ContinuousIndex< TIndexRep, SliceImageType::ImageDimension > sliceIndex;
    for( unsigned int ii = 0; ii < SliceImageType::ImageDimension; ++ii )
      {
      sliceIndex[ii] = index[ii];
      }
    Point< TCoordRep, SliceImageType::ImageDimension > slicePoint;
    this->m_SliceImage->TransformContinuousIndexToPhysicalPoint( sliceIndex, slicePoint );
    for( unsigned int ii = 0; ii < SliceImageType::ImageDimension; ++ii )
      {
      point[ii] += slicePoint[ii];
      }
    typedef Point< TCoordRep, VDimension > PointType;
    PointType lowerPoint;
    const IndexValueType floor = Math::Floor< IndexValueType, TIndexRep >( index[ImageDimension - 1] );
    const IndexValueType ceil = Math::Ceil< IndexValueType, TIndexRep >( index[ImageDimension - 1] );
    const TransformType * transform = this->GetSliceTransform( floor );
    if( transform != ITK_NULLPTR )
      {
      lowerPoint = transform->TransformPoint( point );
      }
    else
      {
      const RegionType & largestRegion = this->GetLargestPossibleRegion();
      const IndexType & largestIndex = largestRegion.GetIndex();
      if( index[ImageDimension - 1] < largestIndex[ImageDimension - 1] )
        {
        point[ImageDimension - 1] = index[ImageDimension - 1] - largestIndex[ImageDimension - 1];
        point = transform->TransformPoint( point );
        return;
        }

      const SizeType & largestSize = largestRegion.GetSize();
      if( index[ImageDimension - 1] > largestIndex[ImageDimension - 1] + largestSize[ImageDimension - 1] - 1 )
        {
        point[ImageDimension - 1] = index[ImageDimension - 1] - largestIndex[ImageDimension - 1] + largestSize[ImageDimension - 1] - 1;
        point = transform->TransformPoint( point );
        }
      return;
      }

    transform = this->GetSliceTransform( ceil );
    PointType upperPoint;
    if( transform != ITK_NULLPTR )
      {
      upperPoint = transform->TransformPoint( point );
      }
    else
      {
      const RegionType & largestRegion = this->GetLargestPossibleRegion();
      const IndexType & largestIndex = largestRegion.GetIndex();
      if( index[ImageDimension - 1] < largestIndex[ImageDimension - 1] )
        {
        point[ImageDimension - 1] = index[ImageDimension - 1] - largestIndex[ImageDimension - 1];
        point = transform->TransformPoint( point );
        return;
        }

      const SizeType & largestSize = largestRegion.GetSize();
      if( index[ImageDimension - 1] > largestIndex[ImageDimension - 1] + largestSize[ImageDimension - 1] - 1 )
        {
        point[ImageDimension - 1] = index[ImageDimension - 1] - largestIndex[ImageDimension - 1] + largestSize[ImageDimension - 1] - 1;
        point = transform->TransformPoint( point );
        }
      return;
      }

    const TIndexRep fraction = index[ImageDimension - 1] - floor;
    for( unsigned int ii = 0; ii < ImageDimension; ++ii )
      {
      point[ii] = lowerPoint[ii] + fraction * (upperPoint[ii] - lowerPoint[ii]);
      }
  }

  /** Get a physical point (in the space which
   * the origin and spacing information comes from)
   * from a discrete index (in the index space)
   *
   * \sa Transform */
  template< typename TCoordRep >
  void TransformIndexToPhysicalPoint(
    const IndexType & index,
    Point< TCoordRep, VDimension > & point) const
  {
    point.Fill( 0.0 );
    typename SliceImageType::IndexType sliceIndex;
    for( unsigned int ii = 0; ii < SliceImageType::ImageDimension; ++ii )
      {
      sliceIndex[ii] = index[ii];
      }
    Point< TCoordRep, SliceImageType::ImageDimension > slicePoint;
    this->m_SliceImage->TransformIndexToPhysicalPoint( sliceIndex, slicePoint );
    for( unsigned int ii = 0; ii < SliceImageType::ImageDimension; ++ii )
      {
      point[ii] += slicePoint[ii];
      }
    const TransformType * transform = this->GetSliceTransform( index[ImageDimension - 1] );
    if( transform != ITK_NULLPTR )
      {
      point = transform->TransformPoint( point );
      return;
      }
    const RegionType & largestRegion = this->GetLargestPossibleRegion();
    const IndexType & largestIndex = largestRegion.GetIndex();
    if( index[ImageDimension - 1] < largestIndex[ImageDimension - 1] )
      {
      point[ImageDimension - 1] = index[ImageDimension - 1] - largestIndex[ImageDimension - 1];
      transform = this->GetSliceTransform( largestIndex[ImageDimension - 1] );
      point = transform->TransformPoint( point );
      return;
      }

    const SizeType & largestSize = largestRegion.GetSize();
    if( index[ImageDimension - 1] > static_cast< IndexValueType >( largestIndex[ImageDimension - 1] + largestSize[ImageDimension - 1] - 1 ) )
      {
      point[ImageDimension - 1] = index[ImageDimension - 1] - largestIndex[ImageDimension - 1] + largestSize[ImageDimension - 1] - 1;
      transform = this->GetSliceTransform( largestIndex[ImageDimension - 1] + largestSize[ImageDimension - 1] - 1 );
      point = transform->TransformPoint( point );
      }
  }

  template< typename TCoordRep >
  void TransformLocalVectorToPhysicalVector(
    FixedArray< TCoordRep, VDimension > & ) const
    {}

  template< typename TCoordRep >
  void TransformPhysicalVectorToLocalVector(
    const FixedArray< TCoordRep, VDimension > & ,
    FixedArray< TCoordRep, VDimension > & ) const
    {}

  /** Return the Pixel Accessor object */
  AccessorType GetPixelAccessor(void)
  { return AccessorType(); }

  /** Return the Pixel Accesor object */
  const AccessorType GetPixelAccessor(void) const
  { return AccessorType(); }

  /** Return the NeighborhoodAccessor functor */
  NeighborhoodAccessorFunctorType GetNeighborhoodAccessor()
  { return NeighborhoodAccessorFunctorType(); }

  /** Return the NeighborhoodAccessor functor */
  const NeighborhoodAccessorFunctorType GetNeighborhoodAccessor() const
  { return NeighborhoodAccessorFunctorType(); }

  virtual void SetLargestPossibleRegion(const RegionType & region) override;

  virtual void CopyInformation( const DataObject * data ) override;

protected:
  SliceSeriesSpecialCoordinatesImage();

  virtual ~SliceSeriesSpecialCoordinatesImage() {}
  virtual void PrintSelf(std::ostream & os, Indent indent) const override;

  const TransformType * GetSliceInverseTransform( IndexValueType sliceIndex ) const;

private:
  SliceSeriesSpecialCoordinatesImage(const Self &); // purposely not implemented
  void operator=(const Self &);                          // purposely not implemented

  typename SliceImageType::Pointer      m_SliceImage;
  typename SliceTransformsType::Pointer m_SliceTransforms;
  typename SliceTransformsType::Pointer m_SliceInverseTransforms;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSliceSeriesSpecialCoordinatesImage.hxx"
#endif

#endif
