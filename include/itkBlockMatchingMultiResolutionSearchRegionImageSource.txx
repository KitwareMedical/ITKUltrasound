#ifndef __itkBlockMatchingMultiResolutionSearchRegionImageSource_txx
#define __itkBlockMatchingMultiResolutionSearchRegionImageSource_txx

#include "itkBlockMatchingMultiResolutionSearchRegionImageSource.h"

namespace itk
{
namespace BlockMatching
{

template < class TFixedImage, class TMovingImage, class TDisplacementImage >
MultiResolutionSearchRegionImageSource< TFixedImage, TMovingImage,
  TDisplacementImage >
::MultiResolutionSearchRegionImageSource():
  m_CurrentLevel( 0 )
{
  m_PreviousDisplacements  = DisplacementImageType::New();
  m_DisplacementDuplicator = DisplacementDuplicatorType::New(); 
  m_DisplacementResampler  = DisplacementResamplerType::New();
}

template < class TFixedImage, class TMovingImage, class TDisplacementImage >
void
MultiResolutionSearchRegionImageSource< TFixedImage, TMovingImage,
                                        TDisplacementImage >
::SetPreviousDisplacements( const DisplacementImageType* displacements )
{
  // We make a copy because the ImageRegistrationMethods will mess with the
  // output information on GenerateOutputInformation before the
  // SearchRegionImageSource has a chance to use it in GenerateData.
  m_DisplacementDuplicator->SetInputImage( displacements );
  m_DisplacementDuplicator->Update();
  m_PreviousDisplacements = m_DisplacementDuplicator->GetOutput();
}

template < class TFixedImage, class TMovingImage, class TDisplacementImage >
void
MultiResolutionSearchRegionImageSource< TFixedImage, TMovingImage,
                                        TDisplacementImage >
::SetOverlapSchedule( const double& schedule )
{
  // Check to make sure the PyramidSchedule has been set.
  if( m_PyramidSchedule.size() == 0 )
    {
    itkExceptionMacro(<<"The PyramidSchedule must be set before calling this method.");
    }

  m_OverlapSchedule.SetSize( m_PyramidSchedule.rows(), m_PyramidSchedule.cols() );
  m_OverlapSchedule.Fill( schedule );
}

template < class TFixedImage, class TMovingImage, class TDisplacement >
void
MultiResolutionSearchRegionImageSource< TFixedImage, TMovingImage,
                                        TDisplacement >
::GenerateOutputInformation()
{
  typename OutputImageType::Pointer outputPtr = this->GetOutput();
  if( !outputPtr )
    return;

  if( m_FixedImage.GetPointer() == NULL )
    {
    itkExceptionMacro( << "Fixed Image is not present." );
    }
  m_FixedImage->UpdateOutputInformation();
  outputPtr->SetDirection( m_FixedImage->GetDirection() );

  // Set origin.  The first block is in the corner of the fixed image.
  typename FixedImageType::PointType  fixedOrigin = m_FixedImage->GetOrigin();
  typename OutputImageType::PointType origin;

  typename FixedImageType::SpacingType fixedSpacing = m_FixedImage->GetSpacing();

  RadiusType nullRadius;
  nullRadius.Fill( 0 );
  if( m_FixedBlockRadius == nullRadius )
    {
    itkExceptionMacro( << "The FixedBlockRadius has not been set." );
    }

  if( m_OverlapSchedule.size() == 0 )
    {
    itkExceptionMacro( << "OverlapSchedule is not present." );
    }

  unsigned int i;
  typename FixedImageType::IndexType fixedIndex = m_FixedImage->GetLargestPossibleRegion().GetIndex();
  for( i = 0; i < ImageDimension; i++ )
    {
    origin[i] = fixedOrigin[i] + 
      + ( fixedIndex[i] + m_FixedBlockRadius[i] ) * fixedSpacing[i];
    }
  outputPtr->SetOrigin( origin );

  typename OutputImageType::SpacingType spacing;
  for( i = 0; i < ImageDimension; i++ )
    {
    spacing[i] = 2 * m_FixedBlockRadius[i] * fixedSpacing[i] * m_OverlapSchedule( m_CurrentLevel, i );
    }
  outputPtr->SetSpacing( spacing );

  OutputRegionType                     region;
  typename OutputRegionType::IndexType index;
  index.Fill( 0 );
  region.SetIndex( index );
  typename OutputRegionType::SizeType size;
  typename FixedImageType::SizeType fixedSize = m_FixedImage->GetLargestPossibleRegion().GetSize();
  for( i = 0; i < ImageDimension; i++ )
    {
    size[i] = static_cast< typename OutputRegionType::SizeType::SizeValueType >( vcl_floor((
    fixedSize[i] - m_FixedBlockRadius[i] - 1 ) * fixedSpacing[i] / spacing[i] ));
    }
  region.SetSize( size );
  outputPtr->SetLargestPossibleRegion( region );
}

template < class TFixedImage, class TMovingImage, class TDisplacement >
void
MultiResolutionSearchRegionImageSource< TFixedImage, TMovingImage,
                                        TDisplacement >
::BeforeThreadedGenerateData()
{
  if( this->m_CurrentLevel != 0 )
    {
    // ! @todo these resampler should be replaced by resamplers for each
    // component that can specify a neumann boundary condition
    // ditto with FixedSearchRegionImageSource
    m_DisplacementResampler->SetInput( this->m_PreviousDisplacements );
    typename OutputImageType::Pointer outputPtr = this->GetOutput();
    if( !outputPtr )
      {
      return;
      }
    m_DisplacementResampler->SetSize(             outputPtr->GetRequestedRegion().GetSize() );
    m_DisplacementResampler->SetOutputStartIndex( outputPtr->GetRequestedRegion().GetIndex() );
    m_DisplacementResampler->SetOutputSpacing(    outputPtr->GetSpacing() );
    m_DisplacementResampler->SetOutputOrigin(     outputPtr->GetOrigin() );
    m_DisplacementResampler->SetOutputDirection(  outputPtr->GetDirection() );
    m_DisplacementResampler->UpdateLargestPossibleRegion();
    }
}

} // end namespace itk
} // end namespace BlockMatching

#endif
