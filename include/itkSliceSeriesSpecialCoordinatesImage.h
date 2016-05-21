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
class SliceSeriesSpecialCoordinatesImage:
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

  void SetSliceTransform( SizeValueType sliceIndex, TransformType * transform );
  const TransformType * GetSliceTransform( SizeValueType sliceIndex ) const;

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
  virtual void Graft(const DataObject *data) ITK_OVERRIDE;

  /** \brief Get the continuous index from a physical point
   *
   * Returns true if the resulting index is within the image, false otherwise.
   * \sa Transform */
  template< typename TCoordRep, typename TIndexRep >
  bool TransformPhysicalPointToContinuousIndex(
    const Point< TCoordRep, VDimension > & point,
    ContinuousIndex< TIndexRep, VDimension > & index) const
  {
    //const RegionType & region = this->GetLargestPossibleRegion();
    //const double maxLateral = region.GetSize(1) - 1;

    //// Convert Cartesian coordinates into angular coordinates
    //TCoordRep lateral = Math::pi_over_2;
    //if( point[1] != 0.0 )
      //{
      //lateral = std::atan(point[0] / point[1]);
      //}
    //const TCoordRep radius  = std::sqrt(point[0] * point[0] + point[1] * point[1] );

    //// Convert the "proper" angular coordinates into index format
    //index[0] = static_cast< TCoordRep >( ( ( radius - m_FirstSampleDistance )
                                           /// m_RadiusSampleSize ) );
    //index[1] = static_cast< TCoordRep >( ( lateral / m_LateralAngularSeparation )
                                         //+ ( maxLateral / 2.0 ) );
    //Vector< SpacePrecisionType, VDimension > cvector;
    //for ( unsigned int kk = 0; kk < VDimension; ++kk )
      //{
      //cvector[kk] = point[kk] - this->m_Origin[kk];
      //}
    //cvector = this->m_PhysicalPointToIndex * cvector;
    //for ( unsigned int ii = 2; ii < VDimension; ++ii )
      //{
      //index[ii] = static_cast< TIndexRep >( cvector[ii] );
      //}

    //// Now, check to see if the index is within allowed bounds
    //const bool isInside = region.IsInside(index);

    //return isInside;
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
    //const RegionType & region = this->GetLargestPossibleRegion();
    //const double maxLateral = region.GetSize(1) - 1;

    //// Convert Cartesian coordinates into angular coordinates
    //TCoordRep lateral = Math::pi_over_2;
    //if( point[1] != 0.0 )
      //{
      //lateral = std::atan(point[0] / point[1]);
      //}
    //const TCoordRep radius  = std::sqrt(point[0] * point[0] + point[1] * point[1] );

    //// Convert the "proper" angular coordinates into index format
    //index[0] = static_cast< IndexValueType >( ( ( radius - m_FirstSampleDistance )
                                           /// m_RadiusSampleSize ) );
    //index[1] = static_cast< IndexValueType >( ( lateral / m_LateralAngularSeparation )
                                         //+ ( maxLateral / 2.0 ) );
    //for ( unsigned int ii = 2; ii < VDimension; ++ii )
      //{
      //TCoordRep sum = NumericTraits< TCoordRep >::ZeroValue();
      //for ( unsigned int jj = 0; jj < VDimension; ++jj )
        //{
        //sum += this->m_PhysicalPointToIndex[ii][jj] * ( point[jj] - this->m_Origin[jj] );
        //}
      //index[ii] = Math::RoundHalfIntegerUp< IndexValueType >(sum);
      //}

    //// Now, check to see if the index is within allowed bounds
    //const bool isInside = region.IsInside(index);

    //return isInside;
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
    typename SliceImageType::ContinuousIndexType sliceIndex;
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
    const TransformType * transform = this->GetSliceTransform( Math::RoundHalfIntegerUp< IndexValueType, TIndexRep >( index[ImageDimension - 1] ) );
    if( transform != ITK_NULLPTR )
      {
      point = transform->TransformPoint( point );
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

  virtual void SetLargestPossibleRegion(const RegionType & region) ITK_OVERRIDE;

protected:
  SliceSeriesSpecialCoordinatesImage();

  virtual ~SliceSeriesSpecialCoordinatesImage() {}
  virtual void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

  const TransformType * GetSliceInverseTransform( SizeValueType sliceIndex ) const;

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
