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
#ifndef itkCurvilinearArraySpecialCoordinatesImage_h
#define itkCurvilinearArraySpecialCoordinatesImage_h

#include "itkSpecialCoordinatesImage.h"
#include "itkPoint.h"
#include "vnl/vnl_math.h"
#include "itkNeighborhoodAccessorFunctor.h"

namespace itk
{
/** \class CurvilinearArraySpecialCoordinatesImage
 *
 *  \brief Templated 2D nonrectilinear-coordinate image class for
 *  curvilinear/phased-array "range" images.
 *
 * \verbatim
 *
 *                             +---------------------> x-axis
 *                             |\
 *                          /  | \
 *                             |-~\
 *                       /     |   \
 *                             |    \
 *                    /        |     \
 *                             | lateral
 *                             |
 *                             v y-axis
 *
 * \endverbatim
 *
 * The equations form performing the conversion from Cartesian coordinates to
 * curvilinear/phased array coordinates are as follows:
 *
 * lateral = arctan(x/y)
 * radius = std::sqrt(x^2 + y^2)
 *
 * The reversed transforms are:
 *
 * x = radius * std::sin(lateral)
 * y = radius * std::cos(lateral)
 *
 * CurvilinearArraySpecialCoordinatesImages are templated over a pixel
 * type and follow the SpecialCoordinatesImage interface.  The data in
 * an image is  arranged in a 1D array as
 * [lateral-index][radius-index] with radius-index
 * varying most rapidly.  The Index type reverses the order so that
 * Index[0] = radius-index, Index[1] = lateral-index.
 *
 * Lateral is discretized into m_LateralAngularSeparation intervals
 * per angular voxel, the most negative lateral interval containing
 * data is then mapped to lateral-index=0, and the largest lateral
 * interval containing data is then mapped to lateral-index=( number
 * of samples along lateral axis - 1 ). Elevation is discretized in
 * the same manner.  This way, the mapping to Cartesian space is
 * symmetric about the x axis such that the line defined by
 * lateral/2 = x-axis.  Radius is discretized into
 * m_RadiusSampleSize units per angular voxel.  The smallest range
 * interval containing data is then mapped to radius-index=0, such
 * that radius = m_FirstSampleDistance + (radius-index *
 * m_RadiusSampleSize).
 *
 * \sa SpecialCoordinatesImage
 * \sa PhasedArray3DSpecialCoordinatesImage
 * \sa SliceSeriesSpecialCoordinatesImage
 *
 * \ingroup Ultrasound
 *
 * \ingroup ImageObjects
 */
template <typename TPixel, unsigned int VDimension>
class ITK_TEMPLATE_EXPORT CurvilinearArraySpecialCoordinatesImage : public SpecialCoordinatesImage<TPixel, VDimension>
{
public:
  /** Standard class type alias */
  using Self = CurvilinearArraySpecialCoordinatesImage;
  using Superclass = SpecialCoordinatesImage<TPixel, VDimension>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;
  using ConstWeakPointer = WeakPointer<const Self>;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(CurvilinearArraySpecialCoordinatesImage, SpecialCoordinatesImage);

  /** Pixel type alias support. Used to declare pixel type in filters
   * or other operations. */
  using PixelType = TPixel;

  /** Typedef alias for PixelType */
  using ValueType = TPixel;

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
    const RegionType & region = this->GetLargestPossibleRegion();
    const double       maxLateral = region.GetSize(1) - 1;

    // Convert Cartesian coordinates into angular coordinates
    TCoordRep lateral = Math::pi_over_2;
    if (point[1] != 0.0)
    {
      lateral = std::atan(point[0] / point[1]);
    }
    const TCoordRep radius = std::sqrt(point[0] * point[0] + point[1] * point[1]);

    // Convert the "proper" angular coordinates into index format
    index[0] = static_cast<TCoordRep>(((radius - m_FirstSampleDistance) / m_RadiusSampleSize));
    index[1] = static_cast<TCoordRep>((lateral / m_LateralAngularSeparation) + (maxLateral / 2.0));
    Vector<SpacePrecisionType, VDimension> cvector;
    for (unsigned int kk = 0; kk < VDimension; ++kk)
    {
      cvector[kk] = point[kk] - this->m_Origin[kk];
    }
    cvector = this->m_PhysicalPointToIndex * cvector;
    for (unsigned int ii = 2; ii < VDimension; ++ii)
    {
      index[ii] = static_cast<TIndexRep>(cvector[ii]);
    }

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
    const RegionType & region = this->GetLargestPossibleRegion();
    const double       maxLateral = region.GetSize(1) - 1;

    // Convert Cartesian coordinates into angular coordinates
    TCoordRep lateral = Math::pi_over_2;
    if (point[1] != 0.0)
    {
      lateral = std::atan(point[0] / point[1]);
    }
    const TCoordRep radius = std::sqrt(point[0] * point[0] + point[1] * point[1]);

    // Convert the "proper" angular coordinates into index format
    index[0] = static_cast<IndexValueType>(((radius - m_FirstSampleDistance) / m_RadiusSampleSize));
    index[1] = static_cast<IndexValueType>((lateral / m_LateralAngularSeparation) + (maxLateral / 2.0));
    for (unsigned int ii = 2; ii < VDimension; ++ii)
    {
      TCoordRep sum = NumericTraits<TCoordRep>::ZeroValue();
      for (unsigned int jj = 0; jj < VDimension; ++jj)
      {
        sum += this->m_PhysicalPointToIndex[ii][jj] * (point[jj] - this->m_Origin[jj]);
      }
      index[ii] = Math::RoundHalfIntegerUp<IndexValueType>(sum);
    }

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
    const RegionType & region = this->GetLargestPossibleRegion();
    const double       maxLateral = region.GetSize(1) - 1;

    // Convert the index into proper angular coordinates
    const TCoordRep radius = (index[0] * m_RadiusSampleSize) + m_FirstSampleDistance;
    const TCoordRep lateral = (index[1] - (maxLateral / 2.0)) * m_LateralAngularSeparation;

    // Convert the angular coordinates into Cartesian coordinates
    point[0] = static_cast<TCoordRep>(radius * std::sin(lateral));
    point[1] = static_cast<TCoordRep>(radius * std::cos(lateral));
    for (unsigned int rr = 2; rr < VDimension; ++rr)
    {
      TCoordRep sum = NumericTraits<TCoordRep>::ZeroValue();
      for (unsigned int cc = 0; cc < VDimension; ++cc)
      {
        sum += this->m_IndexToPhysicalPoint(rr, cc) * index[cc];
      }
      point[rr] = sum + this->m_Origin[rr];
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
    const RegionType & region = this->GetLargestPossibleRegion();
    const double       maxLateral = region.GetSize(1) - 1;

    // Convert the index into proper angular coordinates
    const TCoordRep radius = (index[0] * m_RadiusSampleSize) + m_FirstSampleDistance;
    const TCoordRep lateral = (index[1] - (maxLateral / 2.0)) * m_LateralAngularSeparation;

    // Convert the angular coordinates into Cartesian coordinates
    point[0] = static_cast<TCoordRep>(radius * std::sin(lateral));
    point[1] = static_cast<TCoordRep>(radius * std::cos(lateral));
    for (unsigned int ii = 2; ii < VDimension; ++ii)
    {
      point[ii] = this->m_Origin[ii];
      for (unsigned int jj = 0; jj < VDimension; ++jj)
      {
        point[ii] += this->m_IndexToPhysicalPoint[ii][jj] * index[jj];
      }
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

  /** Set/Get the number of radians between each lateral unit.   */
  itkSetMacro(LateralAngularSeparation, double);
  itkGetConstMacro(LateralAngularSeparation, double);

  /** Set/Get the number of cartesian units between each unit along the R .  */
  itkSetMacro(RadiusSampleSize, double);
  itkGetConstMacro(RadiusSampleSize, double);

  /** Set the distance to add to the radius. */
  itkSetMacro(FirstSampleDistance, double);
  itkGetConstMacro(FirstSampleDistance, double);

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
  CopyInformation(const DataObject * data) override;

protected:
  CurvilinearArraySpecialCoordinatesImage()
  {
    m_RadiusSampleSize = 1;
    m_LateralAngularSeparation = 1 * (2.0 * vnl_math::pi / 360.0); // 1 degree
    m_FirstSampleDistance = 0;
  }
  ~CurvilinearArraySpecialCoordinatesImage() override = default;
  
  virtual void
  PrintSelf(std::ostream & os, Indent indent) const override;

private:
  CurvilinearArraySpecialCoordinatesImage(const Self &); // purposely not implemented
  void
  operator=(const Self &); // purposely not implemented

  double m_LateralAngularSeparation; // in radians
  double m_RadiusSampleSize;
  double m_FirstSampleDistance;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itkCurvilinearArraySpecialCoordinatesImage.hxx"
#endif

#endif
