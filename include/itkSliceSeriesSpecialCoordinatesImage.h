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
template <typename TSliceImage,
          typename TTransform,
          typename TPixel = typename TSliceImage::PixelType,
          unsigned int VDimension = TSliceImage::ImageDimension + 1>
class ITK_TEMPLATE_EXPORT SliceSeriesSpecialCoordinatesImage : public SpecialCoordinatesImage<TPixel, VDimension>
{
public:
  /** Standard class type alias */
  using Self = SliceSeriesSpecialCoordinatesImage;
  using Superclass = SpecialCoordinatesImage<TPixel, VDimension>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;
  using ConstWeakPointer = WeakPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SliceSeriesSpecialCoordinatesImage, SpecialCoordinatesImage);

  /** Pixel type alias support. Used to declare pixel type in filters
   * or other operations. */
  using PixelType = TPixel;

  /** Typedef alias for PixelType */
  using ValueType = TPixel;

  /** Typedef of the image slice. This is usually one less dimension of the
   * image.*/
  using SliceImageType = TSliceImage;

  /** Typedef of the transform used to transform each slice. This could also
   * be a base class transform type. */
  using TransformType = TTransform;

  /** Typedef of the array of transform used to store the per-slice
   * transforms. */
  using SliceTransformsType = VectorContainer<IdentifierType, typename TransformType::Pointer>;

  /** Internal Pixel representation. Used to maintain a uniform API
   * with Image Adaptors and allow to keep a particular internal
   * representation of data while showing a different external
   * representation. */
  using InternalPixelType = TPixel;

  using IOPixelType = typename Superclass::IOPixelType;

  /** Accessor type that convert data between internal and external
   *  representations.  */
  using AccessorType = DefaultPixelAccessor<PixelType>;

  /** Accessor functor to choose between accessors: DefaultPixelAccessor for
   * the Image, and DefaultVectorPixelAccessor for the vector image. The
   * functor provides a generic API between the two accessors. */
  using AccessorFunctorType = DefaultPixelAccessorFunctor<Self>;

  /** Typedef for the functor used to access a neighborhood of pixel
   * pointers. */
  using NeighborhoodAccessorFunctorType = NeighborhoodAccessorFunctor<Self>;

  /** Dimension of the image.  This constant is used by functions that are
   * templated over image type (as opposed to being templated over pixel type
   * and dimension) when they need compile time access to the dimension of
   * the image. */
  itkStaticConstMacro(ImageDimension, unsigned int, VDimension);

  /** Index type alias support. An index is used to access pixel values. */
  using IndexType = typename Superclass::IndexType;
  using IndexValueType = typename Superclass::IndexValueType;

  /** Offset type alias support. An offset is used to access pixel values. */
  using OffsetType = typename Superclass::OffsetType;

  /** Size type alias support. A size is used to define region bounds. */
  using SizeType = typename Superclass::SizeType;
  using SizeValueType = typename Superclass::SizeValueType;

  /** Container used to store pixels in the image. */
  using PixelContainer = ImportImageContainer<SizeValueType, PixelType>;

  /** Region type alias support. A region is used to specify a subset of
   *  an image.
   */
  using RegionType = typename Superclass::RegionType;

  /** Spacing type alias support.  Spacing holds the "fake" size of a
   *  pixel, making each pixel look like a 1 unit hyper-cube to filters
   *  that were designed for normal images and that therefore use
   *  m_Spacing.  The spacing is the geometric distance between image
   *  samples.
   */
  using SpacingType = typename Superclass::SpacingType;

  /** Origin type alias support.  The origin is the "fake" geometric
   *  coordinates of the index (0,0).  Also for use w/ filters designed
   *  for normal images.
   */
  using PointType = typename Superclass::PointType;

  /** A pointer to the pixel container. */
  using PixelContainerPointer = typename PixelContainer::Pointer;
  using PixelContainerConstPointer = typename PixelContainer::ConstPointer;

  /** Set/Get the Slice Image. This image is used to find the relative pixel
   * location locations within an image. Its size must be the same size as the
   * the size in the corresponding dimenions of this image. The pixel values
   * in this image are not used. */
  itkSetObjectMacro(SliceImage, SliceImageType);
  itkGetConstObjectMacro(SliceImage, SliceImageType);

  void
  SetSliceTransform(IndexValueType sliceIndex, TransformType * transform);
  const TransformType *
  GetSliceTransform(IndexValueType sliceIndex) const;

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
  virtual void
  Graft(const DataObject * data) override;

  /** \brief Get the continuous index from a physical point
   *
   * Returns true if the resulting index is within the image, false otherwise.
   * \sa Transform */
  template <typename TCoordRep, typename TIndexRep>
  bool
  TransformPhysicalPointToContinuousIndex(const Point<TCoordRep, VDimension> &     point,
                                          ContinuousIndex<TIndexRep, VDimension> & index) const
  {
    const RegionType &    region = this->GetLargestPossibleRegion();
    const unsigned int    sliceDimensionIndex = ImageDimension - 1;
    IndexValueType        lowerIndex = region.GetIndex(sliceDimensionIndex);
    IndexValueType        upperIndex = lowerIndex + region.GetSize(sliceDimensionIndex) - 1;
    IndexValueType        nextIndex = lowerIndex;
    PointType             lowerPoint;
    const TransformType * transform = this->GetSliceInverseTransform(lowerIndex);
    if (transform == nullptr)
    {
      itkExceptionMacro("Inverse slice transform not available for index: " << lowerIndex);
    }
    lowerPoint = transform->TransformPoint(point);
    int       lowerSign = Math::sgn(lowerPoint[sliceDimensionIndex]);
    PointType upperPoint;
    transform = this->GetSliceInverseTransform(upperIndex);
    if (transform == nullptr)
    {
      itkExceptionMacro("Inverse slice transform not available for index: " << upperIndex);
    }
    upperPoint = transform->TransformPoint(point);
    int       upperSign = Math::sgn(upperPoint[sliceDimensionIndex]);
    PointType nextPoint = lowerPoint;
    int       nextSign = 0;
    if (lowerSign == 0)
    {
      nextSign = 0;
      nextPoint = lowerPoint;
    }
    else if (upperSign == 0)
    {
      nextSign = 0;
      nextPoint = upperPoint;
    }
    else
    {
      if (lowerSign == upperSign)
      {
        // outside the image
        return false;
      }

      // Binary search for the transforms that bounds the slice
      while (upperIndex - lowerIndex > 1)
      {
        nextIndex = lowerIndex + (upperIndex - lowerIndex) / 2;
        transform = this->GetSliceInverseTransform(nextIndex);
        nextPoint = transform->TransformPoint(point);
        nextSign = Math::sgn(nextPoint[sliceDimensionIndex]);
        if (nextSign == 0)
        {
          break;
        }
        else if (nextSign == lowerSign)
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

    if (nextSign != 0)
    {
      const TCoordRep fraction =
        -lowerPoint[sliceDimensionIndex] / (upperPoint[sliceDimensionIndex] - lowerPoint[sliceDimensionIndex]);
      nextPoint[sliceDimensionIndex] = lowerIndex + fraction * (upperIndex - lowerIndex);
      for (unsigned int ii = 0; ii < SliceImageType::ImageDimension; ++ii)
      {
        nextPoint[ii] = lowerPoint[ii] + fraction * (upperPoint[ii] - lowerPoint[ii]);
      }
    }
    else
    {
      nextPoint[sliceDimensionIndex] = nextIndex;
    }

    typename SliceImageType::PointType slicePoint;
    for (unsigned int ii = 0; ii < SliceImageType::ImageDimension; ++ii)
    {
      slicePoint[ii] = nextPoint[ii];
    }
    ContinuousIndex<TIndexRep, SliceImageType::ImageDimension> sliceIndex;
    this->m_SliceImage->TransformPhysicalPointToContinuousIndex(slicePoint, sliceIndex);
    for (unsigned int ii = 0; ii < SliceImageType::ImageDimension; ++ii)
    {
      index[ii] = sliceIndex[ii];
    }
    index[sliceDimensionIndex] = nextPoint[sliceDimensionIndex];
    // Now, check to see if the index is within allowed bounds
    const bool isInside = region.IsInside(index);

    return isInside;
  }

  /** \brief Returns the continuous index from a physical point
   * \note This specific overload does not figure out whether or not
   *  the returned index is inside the image. Of course, the user can
   * still test this afterwards by calling ImageRegion::IsInside(index):
     \code
     auto index = image->TransformPhysicalPointToContinuousIndex<double>(point);
     if (image->GetLargestPossibleRegion().IsInside(index)) // Et cetera...
     \endcode
   * Which is equivalent to the following code, which calls the other overload:
     \code
     itk::ContinuousIndex<double, ImageDimension> index;
     if (image->TransformPhysicalPointToContinuousIndex(point, index)) // Et cetera...
     \endcode
   * \sa Transform */
  template <typename TIndexRep, typename TCoordRep>
  ContinuousIndex<TIndexRep, VDimension>
  TransformPhysicalPointToContinuousIndex(const Point<TCoordRep, VDimension> & point) const
  {
    ContinuousIndex<TIndexRep, VDimension> index;
    TransformPhysicalPointToContinuousIndex<TIndexRep>(point, index);
    return index;
  }

  /** Get the index (discrete) from a physical point.
   * Floating point index results are truncated to integers.
   * Returns true if the resulting index is within the image, false otherwise
   * \sa Transform */
  template <typename TCoordRep>
  bool
  TransformPhysicalPointToIndex(const Point<TCoordRep, VDimension> & point, IndexType & index) const
  {
    const RegionType &    region = this->GetLargestPossibleRegion();
    const unsigned int    sliceDimensionIndex = ImageDimension - 1;
    IndexValueType        lowerIndex = region.GetIndex(sliceDimensionIndex);
    IndexValueType        upperIndex = lowerIndex + region.GetSize(sliceDimensionIndex) - 1;
    IndexValueType        nextIndex = lowerIndex;
    PointType             lowerPoint;
    const TransformType * transform = this->GetSliceInverseTransform(lowerIndex);
    if (transform == nullptr)
    {
      itkExceptionMacro("Inverse slice transform not available for index: " << lowerIndex);
    }
    lowerPoint = transform->TransformPoint(point);
    int       lowerSign = Math::sgn(lowerPoint[sliceDimensionIndex]);
    PointType upperPoint;
    transform = this->GetSliceInverseTransform(upperIndex);
    if (transform == nullptr)
    {
      itkExceptionMacro("Inverse slice transform not available for index: " << upperIndex);
    }
    upperPoint = transform->TransformPoint(point);
    int       upperSign = Math::sgn(upperPoint[sliceDimensionIndex]);
    PointType nextPoint = lowerPoint;
    int       nextSign = 0;
    if (lowerSign == 0)
    {
      nextSign = 0;
      nextPoint = lowerPoint;
    }
    else if (upperSign == 0)
    {
      nextSign = 0;
      nextPoint = upperPoint;
    }
    else
    {
      if (lowerSign == upperSign)
      {
        // outside the image
        return false;
      }

      // Binary search for the transforms that bounds the slice
      while (upperIndex - lowerIndex > 1)
      {
        nextIndex = lowerIndex + (upperIndex - lowerIndex) / 2;
        transform = this->GetSliceInverseTransform(nextIndex);
        nextPoint = transform->TransformPoint(point);
        nextSign = Math::sgn(nextPoint[sliceDimensionIndex]);
        if (nextSign == 0)
        {
          break;
        }
        else if (nextSign == lowerSign)
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

    if (nextSign != 0)
    {
      const TCoordRep fraction =
        -lowerPoint[sliceDimensionIndex] / (upperPoint[sliceDimensionIndex] - lowerPoint[sliceDimensionIndex]);
      nextPoint[sliceDimensionIndex] = lowerIndex + fraction * (upperIndex - lowerIndex);
      for (unsigned int ii = 0; ii < SliceImageType::ImageDimension; ++ii)
      {
        nextPoint[ii] = lowerPoint[ii] + fraction * (upperPoint[ii] - lowerPoint[ii]);
      }
    }
    else
    {
      nextPoint[sliceDimensionIndex] = nextIndex;
    }

    typename SliceImageType::PointType slicePoint;
    for (unsigned int ii = 0; ii < SliceImageType::ImageDimension; ++ii)
    {
      slicePoint[ii] = nextPoint[ii];
    }
    typename SliceImageType::IndexType sliceIndex;
    this->m_SliceImage->TransformPhysicalPointToIndex(slicePoint, sliceIndex);
    for (unsigned int ii = 0; ii < SliceImageType::ImageDimension; ++ii)
    {
      index[ii] = sliceIndex[ii];
    }
    index[sliceDimensionIndex] = Math::RoundHalfIntegerUp<IndexValueType>(nextPoint[sliceDimensionIndex]);
    // Now, check to see if the index is within allowed bounds
    const bool isInside = region.IsInside(index);

    return isInside;
  }

  /** Returns the index (discrete) of a voxel from a physical point.
   * Floating point index results are rounded to integers
   * \note This specific overload does not figure out whether or not
   * the returned index is inside the image. Of course, the user can
   * still test this afterwards by calling ImageRegion::IsInside(index):
     \code
      auto index = image->TransformPhysicalPointToIndex(point);
      if (image->GetLargestPossibleRegion().IsInside(index)) // Et cetera...
     \endcode
   * Which is equivalent to the following code, which calls the other overload:
     \code
      IndexType index;
      if (image->TransformPhysicalPointToIndex(point, index)) // Et cetera...
     \endcode
   * \sa Transform */
  template <typename TCoordRep>
  IndexType
  TransformPhysicalPointToIndex(const Point<TCoordRep, VDimension> & point) const
  {
    IndexType index;
    TransformPhysicalPointToIndex(point, index);
    return index;
  }

  /** Get a physical point (in the space which
   * the origin and spacing information comes from)
   * from a continuous index (in the index space)
   * \sa Transform */
  template <typename TCoordRep, typename TIndexRep>
  void
  TransformContinuousIndexToPhysicalPoint(const ContinuousIndex<TIndexRep, VDimension> & index,
                                          Point<TCoordRep, VDimension> &                 point) const
  {
    point.Fill(0.0);
    ContinuousIndex<TIndexRep, SliceImageType::ImageDimension> sliceIndex;
    for (unsigned int ii = 0; ii < SliceImageType::ImageDimension; ++ii)
    {
      sliceIndex[ii] = index[ii];
    }
    Point<TCoordRep, SliceImageType::ImageDimension> slicePoint;
    this->m_SliceImage->TransformContinuousIndexToPhysicalPoint(sliceIndex, slicePoint);
    for (unsigned int ii = 0; ii < SliceImageType::ImageDimension; ++ii)
    {
      point[ii] += slicePoint[ii];
    }
    using PointType = Point<TCoordRep, VDimension>;
    PointType             lowerPoint;
    const IndexValueType  floor = Math::Floor<IndexValueType, TIndexRep>(index[ImageDimension - 1]);
    const IndexValueType  ceil = Math::Ceil<IndexValueType, TIndexRep>(index[ImageDimension - 1]);
    const TransformType * transform = this->GetSliceTransform(floor);
    if (transform != nullptr)
    {
      lowerPoint = transform->TransformPoint(point);
    }
    else
    {
      const RegionType & largestRegion = this->GetLargestPossibleRegion();
      const IndexType &  largestIndex = largestRegion.GetIndex();
      if (index[ImageDimension - 1] < largestIndex[ImageDimension - 1])
      {
        point[ImageDimension - 1] = index[ImageDimension - 1] - largestIndex[ImageDimension - 1];
        point = transform->TransformPoint(point);
        return;
      }

      const SizeType & largestSize = largestRegion.GetSize();
      if (index[ImageDimension - 1] > largestIndex[ImageDimension - 1] + largestSize[ImageDimension - 1] - 1)
      {
        point[ImageDimension - 1] =
          index[ImageDimension - 1] - largestIndex[ImageDimension - 1] + largestSize[ImageDimension - 1] - 1;
        point = transform->TransformPoint(point);
      }
      return;
    }

    transform = this->GetSliceTransform(ceil);
    PointType upperPoint;
    if (transform != nullptr)
    {
      upperPoint = transform->TransformPoint(point);
    }
    else
    {
      const RegionType & largestRegion = this->GetLargestPossibleRegion();
      const IndexType &  largestIndex = largestRegion.GetIndex();
      if (index[ImageDimension - 1] < largestIndex[ImageDimension - 1])
      {
        point[ImageDimension - 1] = index[ImageDimension - 1] - largestIndex[ImageDimension - 1];
        point = transform->TransformPoint(point);
        return;
      }

      const SizeType & largestSize = largestRegion.GetSize();
      if (index[ImageDimension - 1] > largestIndex[ImageDimension - 1] + largestSize[ImageDimension - 1] - 1)
      {
        point[ImageDimension - 1] =
          index[ImageDimension - 1] - largestIndex[ImageDimension - 1] + largestSize[ImageDimension - 1] - 1;
        point = transform->TransformPoint(point);
      }
      return;
    }

    const TIndexRep fraction = index[ImageDimension - 1] - floor;
    for (unsigned int ii = 0; ii < ImageDimension; ++ii)
    {
      point[ii] = lowerPoint[ii] + fraction * (upperPoint[ii] - lowerPoint[ii]);
    }
  }

  /** Returns a physical point (in the space which
   * the origin and spacing information comes from)
   * from a continuous index (in the index space)
   * \sa Transform */
  template <typename TCoordRep, typename TIndexRep>
  Point<TCoordRep, VDimension>
  TransformContinuousIndexToPhysicalPoint(const ContinuousIndex<TIndexRep, VDimension> & index) const
  {
    Point<TCoordRep, VDimension> point;
    TransformContinuousIndexToPhysicalPoint(index, point);
    return point;
  }

  /** Get a physical point (in the space which
   * the origin and spacing information comes from)
   * from a discrete index (in the index space)
   *
   * \sa Transform */
  template <typename TCoordRep>
  void
  TransformIndexToPhysicalPoint(const IndexType & index, Point<TCoordRep, VDimension> & point) const
  {
    point.Fill(0.0);
    typename SliceImageType::IndexType sliceIndex;
    for (unsigned int ii = 0; ii < SliceImageType::ImageDimension; ++ii)
    {
      sliceIndex[ii] = index[ii];
    }
    Point<TCoordRep, SliceImageType::ImageDimension> slicePoint;
    this->m_SliceImage->TransformIndexToPhysicalPoint(sliceIndex, slicePoint);
    for (unsigned int ii = 0; ii < SliceImageType::ImageDimension; ++ii)
    {
      point[ii] += slicePoint[ii];
    }
    const TransformType * transform = this->GetSliceTransform(index[ImageDimension - 1]);
    if (transform != nullptr)
    {
      point = transform->TransformPoint(point);
      return;
    }
    const RegionType & largestRegion = this->GetLargestPossibleRegion();
    const IndexType &  largestIndex = largestRegion.GetIndex();
    if (index[ImageDimension - 1] < largestIndex[ImageDimension - 1])
    {
      point[ImageDimension - 1] = index[ImageDimension - 1] - largestIndex[ImageDimension - 1];
      transform = this->GetSliceTransform(largestIndex[ImageDimension - 1]);
      point = transform->TransformPoint(point);
      return;
    }

    const SizeType & largestSize = largestRegion.GetSize();
    if (index[ImageDimension - 1] >
        static_cast<IndexValueType>(largestIndex[ImageDimension - 1] + largestSize[ImageDimension - 1] - 1))
    {
      point[ImageDimension - 1] =
        index[ImageDimension - 1] - largestIndex[ImageDimension - 1] + largestSize[ImageDimension - 1] - 1;
      transform = this->GetSliceTransform(largestIndex[ImageDimension - 1] + largestSize[ImageDimension - 1] - 1);
      point = transform->TransformPoint(point);
    }
  }

  /** Returns a physical point (in the space which
   * the origin and spacing information comes from)
   * from a discrete index (in the index space)
   *
   * \sa Transform */
  template <typename TCoordRep>
  Point<TCoordRep, VDimension>
  TransformIndexToPhysicalPoint(const IndexType & index) const
  {
    Point<TCoordRep, VDimension> point;
    TransformIndexToPhysicalPoint(index, point);
    return point;
  }

  template <typename TCoordRep>
  void
  TransformLocalVectorToPhysicalVector(FixedArray<TCoordRep, VDimension> &) const
  {}

  template <typename TCoordRep>
  void
  TransformPhysicalVectorToLocalVector(const FixedArray<TCoordRep, VDimension> &,
                                       FixedArray<TCoordRep, VDimension> &) const
  {}

  /** Return the Pixel Accessor object */
  AccessorType
  GetPixelAccessor(void)
  {
    return AccessorType();
  }

  /** Return the Pixel Accesor object */
  const AccessorType
  GetPixelAccessor(void) const
  {
    return AccessorType();
  }

  /** Return the NeighborhoodAccessor functor */
  NeighborhoodAccessorFunctorType
  GetNeighborhoodAccessor()
  {
    return NeighborhoodAccessorFunctorType();
  }

  /** Return the NeighborhoodAccessor functor */
  const NeighborhoodAccessorFunctorType
  GetNeighborhoodAccessor() const
  {
    return NeighborhoodAccessorFunctorType();
  }

  virtual void
  SetLargestPossibleRegion(const RegionType & region) override;

  virtual void
  CopyInformation(const DataObject * data) override;

protected:
  SliceSeriesSpecialCoordinatesImage();
  ~SliceSeriesSpecialCoordinatesImage() override = default;

  virtual void
  PrintSelf(std::ostream & os, Indent indent) const override;

  const TransformType *
  GetSliceInverseTransform(IndexValueType sliceIndex) const;

private:
  SliceSeriesSpecialCoordinatesImage(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented

  typename SliceImageType::Pointer      m_SliceImage;
  typename SliceTransformsType::Pointer m_SliceTransforms;
  typename SliceTransformsType::Pointer m_SliceInverseTransforms;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkSliceSeriesSpecialCoordinatesImage.hxx"
#endif

#endif
