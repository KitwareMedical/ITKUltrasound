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
#ifndef itkBlockMatchingMetricImageToDisplacementCalculator_h
#define itkBlockMatchingMetricImageToDisplacementCalculator_h

#include "itkImage.h"
#include "itkImageDuplicator.h"
#include "itkMultiThreader.h"

namespace itk
{
namespace BlockMatching
{

/** \class MetricImageToDisplacementCalculator
 *
 * \brief Calculates the displacement of a block from the MetricImage(s).
 *
 * This class tightly integrates with BlockMatching::ImageRegistrationMethod.
 * It is the responsibility of this class to calculate the displacement given
 * MetricImage(s) in a block matching registration method.
 *
 * This is an abstract class with a required, protected  method:
 * void Compute()
 *
 * This method should be overridden and the calculation performed here if the
 * local displacment depends on more than the local MetricImage.  Alternatively,
 * if the local displacement depends only on the metric image,
 * ComputePixel( index, MetricImage* ) should be overridden and Compute()
 * set to do nothing.
 *
 * Caching of the MetricImage can be enabled by SetCacheMetricImageOn();
 *
 * The behavior of the associated BlockMatching::ImageRegistrationMethod
 * GenerateInputRequestedRegion() and EnlargeOutputRequestedRegion() with
 * ModifyGenerateInputRequestedRegion() and
 * ModifyEnlargeOutputRequestedRegion().
 *
 * \ingroup Ultrasound
 * */
template< typename TMetricImage, typename TDisplacementImage >
class ITK_TEMPLATE_EXPORT MetricImageToDisplacementCalculator :
  public Object
{
public:
  /** Standard class typedefs. */
  typedef MetricImageToDisplacementCalculator Self;
  typedef Object                              Superclass;
  typedef SmartPointer<Self>                  Pointer;
  typedef SmartPointer<const Self>            ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( MetricImageToDisplacementCalculator, Object );

  /** Type of the metric image (input pixels). */
  typedef TMetricImage                        MetricImageType;
  typedef typename MetricImageType::Pointer   MetricImagePointerType;
  typedef typename MetricImageType::IndexType IndexType;

  /** Type of the displacement image (output). */
  typedef TDisplacementImage                         DisplacementImageType;
  typedef typename DisplacementImageType::Pointer    DisplacementImagePointerType;
  typedef typename DisplacementImageType::RegionType RegionType;
  typedef typename DisplacementImageType::PointType  PointType;

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TDisplacementImage::ImageDimension);

  /** Type of an image of center points of the fixed image blocks. */
  typedef itk::Image<PointType, ImageDimension>   CenterPointsImageType;
  typedef typename CenterPointsImageType::Pointer CenterPointsImagePointerType;

  /** Type of the metric image image. */
  typedef itk::Image<MetricImagePointerType, ImageDimension>
                                                 MetricImageImageType;
  typedef typename MetricImageImageType::Pointer MetricImageImagePointerType;

  /** Ensure all the metric images are stored in the class's MetricImageImage.
   * */
  itkSetMacro( CacheMetricImage, bool );
  itkGetConstMacro( CacheMetricImage, bool );
  itkBooleanMacro( CacheMetricImage );

  /** Set a metric image pixel.  This is the way to supply input to this class.
   * The point is the center point of the corresponding fixed block.
   * The index is the index of the corresponding displacement image pixel.  The
   * metric image may or may not be stored in the class's MetricImageImage
   * depending on the needs of the implementation and whether
   * SetCacheMetricImageOn() has been set. */
  virtual void SetMetricImagePixel( const PointType & point, const IndexType& index, MetricImageType * image );

  /** Set/Get the displacement image.  Get should only be called after calling
   * Compute() and all the metric image pixels have been set. */
  virtual void SetDisplacementImage( DisplacementImageType * image );

  const DisplacementImageType * GetDisplacementImage() const
  {
    return this->m_DisplacementImage.GetPointer();
  }

  /** Get the MetricImage image (The image of metric images.)  This should only
   * be called after all the the metric image pixels have been set with
   * SetMetricImagePixel(). */
  itkSetObjectMacro( MetricImageImage, MetricImageImageType );
  itkGetConstObjectMacro( MetricImageImage, MetricImageImageType );

  /** Get the center points of the kernels.  The points in the metric image that
   * corresponds to the kernels movement subtracted from these points is the
   * displacement. */
  itkSetObjectMacro( CenterPointsImage, CenterPointsImageType );
  itkGetConstObjectMacro( CenterPointsImage, CenterPointsImageType );

  /** Subclasses must implement this method.  If the displacement calculation takes place in
   * SetMetricImagePixel(), then this can do nothing. */
  virtual void Compute() = 0;

  /** Modify the associated BlockMatching::ImageRegistrationMethod's
   * GenerateInputRequestedRegion().  */
  virtual void ModifyGenerateInputRequestedRegion( RegionType& region ) {}

  /** Modify the associated BlockMatching::ImageRegistrationMethod's
   * EnlargeOutputRequestedRegion().  */
  virtual void ModifyEnlargeOutputRequestedRegion( DataObject* data ) {}

  /** Get/Set the number of threads to create when executing. */
  itkSetClampMacro( NumberOfThreads, ThreadIdType, 1, ITK_MAX_THREADS );
  itkGetConstReferenceMacro( NumberOfThreads, ThreadIdType );

  /** Return the multithreader used by this class. */
  MultiThreader * GetMultiThreader()
  {
    return m_Threader;
  }
protected:
  typedef ImageDuplicator<MetricImageType> MetricImageDuplicatorType;
  MetricImageToDisplacementCalculator();

  /** An inherited class provides the operation in a thread. ThreaderCallback()
   * calls this functor on the thread region.  Watch out when inheriting this --
   * the first argument needs type Superclass and you want to dynamic_cast it
   * down to the inherited class to access private or protected members.
   *
   * \ingroup Ultrasound
   * */
  class ThreadFunctor
  {
  public:
    virtual ITK_THREAD_RETURN_TYPE operator() ( Self *self,
      RegionType& region, ThreadIdType threadId ) = 0;
    virtual ~ThreadFunctor() {};
  };

  /** This is the base thread structure. Populate this with your ThreadFunctor
   * of choice.  */
  struct ThreadStruct
    {
    ThreadFunctor * functor;
    Self * self;
    };

  /** The thread functor of choice is specified here. The functor is applied to
   * the RequestedRegion of the displacement image.  */
  void ApplyThreadFunctor( ThreadFunctor& functor );

  /** Static function used as a "callback" by the MultiThreader.  The threading
   * library will call this routine for each thread, which will delegate the
   * control to the ThreadFunctor. */
  static ITK_THREAD_RETURN_TYPE ThreaderCallback( void *arg );

  /** This is taken from ImageSource and modified.
   *
   * Split the output's RequestedRegion into "num" pieces, returning region "i"
   * as "splitRegion". This method is called "num" times. The regions must not
   * overlap. The method returns the number of pieces that the routine is
   * capable of splitting the output RequestedRegion, i.e. return value is less
   * than or equal to "num". */
  virtual unsigned int SplitRequestedRegion( unsigned int i, unsigned int num, RegionType & splitRegion);

  MultiThreader::Pointer m_Threader;
  ThreadIdType           m_NumberOfThreads;

  CenterPointsImagePointerType m_CenterPointsImage;
  MetricImageImagePointerType  m_MetricImageImage;
  DisplacementImagePointerType m_DisplacementImage;

  bool m_CacheMetricImage;
  bool m_RegionsDefined;

  typename MetricImageDuplicatorType::Pointer m_MetricImageDuplicator;

private:
  ITK_DISALLOW_COPY_AND_ASSIGN(MetricImageToDisplacementCalculator);

};

} // end namespace BlockMatching
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlockMatchingMetricImageToDisplacementCalculator.hxx"
#endif

#endif
