#ifndef __itkBlockMatchingMultiResolutionMinMaxSearchRegionImageSource_h
#define __itkBlockMatchingMultiResolutionMinMaxSearchRegionImageSource_h

#include "itkBlockMatchingMultiResolutionSearchRegionImageSource.h"

namespace itk
{
namespace BlockMatching
{

/** \class MultiResolutionMinMaxSearchRegionImageSource
 *
 * \brief The search region is a factor of the matching block size.
 *
 * The search region size is a factor of the matching block size ( should be 1.0
 * or higher ).  This factor is set at the bottom level at top level and
 * linearly interpolated at in between levels.
 *
 * \todo rename this from 'MinMax' to 'TopBottomFactor'.
 */
template < class TFixedImage, class TMovingImage, class TDisplacementImage >
class ITK_EXPORT MultiResolutionMinMaxSearchRegionImageSource :
  public MultiResolutionSearchRegionImageSource< TFixedImage, TMovingImage, TDisplacementImage >
{
public:
  /** Standard class typedefs. */
  typedef MultiResolutionMinMaxSearchRegionImageSource   Self;
  typedef MultiResolutionSearchRegionImageSource< TFixedImage,
          TMovingImage, TDisplacementImage >             Superclass;
  typedef SmartPointer< Self >                           Pointer;
  typedef SmartPointer< const Self >                     ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( MultiResolutionMinMaxSearchRegionImageSource, MultiResolutionSearchRegionImageSource );

  /** New macro for creation of through a Smart Pointer. */
  itkNewMacro( Self );

  /** ImageDimension enumeration. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TMovingImage::ImageDimension);

  /** Type of the fixed image. */
  typedef typename Superclass::FixedImageType  FixedImageType;
  typedef typename Superclass::FixedRegionType FixedRegionType;

  /** Type of the radius used to characterized the fixed image block. */
  typedef typename Superclass::RadiusType RadiusType;

  /** Type of the moving image. */
  typedef typename Superclass::MovingImageType  MovingImageType;
  typedef typename Superclass::MovingRegionType MovingRegionType;

  /** Type of the search region image. */
  typedef typename Superclass::OutputImageType  OutputImageType;
  typedef typename Superclass::OutputRegionType OutputRegionType;

  /** Type of the displacement image. */
  typedef typename Superclass::DisplacementImageType DisplacementImageType;

  /** ScheduleType typedef support. */
  typedef typename Superclass::PyramidScheduleType PyramidScheduleType;

  /** Type of the search region to block radius ratio. */
  typedef FixedArray< double, ImageDimension > FactorType;

  /** Set the ration of the search region radius to the matching block radius at
   * the bottom level. */
  void SetMinFactor( const FactorType& factor )
    {
    m_MinFactor = factor;
    this->Modified();
    }
  void SetMinFactor( const double& factor )
    {
    FactorType f;
    f.Fill( factor );
    this->SetMinFactor( f );
    }
  itkGetConstReferenceMacro( MinFactor, FactorType );

  /** Set the ration of the search region radius to the matching block radius at
   * the top level. */
  void SetMaxFactor( const FactorType& factor )
    {
    m_MaxFactor = factor;
    this->Modified();
    }
  void SetMaxFactor( const double& factor )
    {
    FactorType f;
    f.Fill( factor );
    this->SetMaxFactor( f );
    }
  itkGetConstReferenceMacro( MaxFactor, FactorType );

protected:
  MultiResolutionMinMaxSearchRegionImageSource():
    m_MinFactor( 1.1 ),
    m_MaxFactor( 3.0 )
    {
    }

  virtual void ThreadedGenerateData( const OutputRegionType& outputRegion,
                                     int threadID );

  FactorType m_MinFactor;
  FactorType m_MaxFactor;

private:
  MultiResolutionMinMaxSearchRegionImageSource( const Self & );
  void operator=( const Self & );

};

} // end namespace BlockMatching
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBlockMatchingMultiResolutionMinMaxSearchRegionImageSource.txx"
#endif

#endif
