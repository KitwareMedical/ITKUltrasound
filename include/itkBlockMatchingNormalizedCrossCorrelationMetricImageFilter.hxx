#ifndef itkBlockMatchingNormalizedCrossCorrelationMetricImageFilter_hxx
#define itkBlockMatchingNormalizedCrossCorrelationMetricImageFilter_hxx

#include "itkBlockMatchingImageToImageMetricMetricImageFilter.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkImageKernelOperator.h"
#include "itkNeighborhoodAlgorithm.h"
#include "itkNeighborhoodInnerProduct.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

namespace itk
{
namespace BlockMatching
{

template <class TFixedImage, class TMovingImage,
          class TMetricImage >
NormalizedCrossCorrelationMetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
::NormalizedCrossCorrelationMetricImageFilter()
{
  // Only #1 is needed by the user in the output.  The others are for subclasses
  // or to simply use the pipeline memory allocation system.
  // 1.  Metric image.
  // 2.  Denominator of the normalized cross correlation coefficient.
  // 3.  Fixed Kernel - Fixed Mean
  // 4.  Moving Search Region - Moving Kernel Mean
  // 5.  Moving image of ones.
  // 6.  FixedPseudoSigmaImage.
  // 7.  FixedMinusMeanSquared.
  this->SetNumberOfOutputs( 7 );
  // ImageSource only does this for the first output.
  for( unsigned int i = 1; i < 7; i++ )
    this->SetNthOutput( i, this->MakeOutput( i ) );

  m_BoxMeanFilter   = BoxMeanFilterType::New();
  m_BoxPseudoSigmaFilter = BoxPseudoSigmaFilterType::New();

  m_BoundaryCondition.SetConstant( NumericTraits< MetricImagePixelType >::Zero );
}

template <class TFixedImage, class TMovingImage,
          class TMetricImage >
void
NormalizedCrossCorrelationMetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
::GenerateOutputInformation()
{
  // We generate the information for the first output ( the metric image );
  Superclass::GenerateOutputInformation();

  FixedImageConstPointerType fixedPtr = this->GetInput( 0 );
  if( !fixedPtr )
    return;

  MovingImageConstPointerType movingPtr = this->GetInput( 1 );
  if( !movingPtr )
    return;

  // Then we copy the information to all the other outputs.
  MetricImagePointerType metricPtr =  this->GetOutput();
  if( !metricPtr )
    {
    return;
    }

  if( !this->m_MovingImageRegionDefined )
    {
    itkExceptionMacro( << "Moving image Region has not been set." );
    }

  // Denominator of normalized cross correlation coefficient.
  MetricImagePointerType outputPtr;
  outputPtr = this->GetOutput( 1 );
  outputPtr->CopyInformation( movingPtr );
  outputPtr->SetRegions( this->m_MovingImageRegion );

  // Fixed kernel - fixed mean.
  outputPtr = this->GetOutput( 2 );
  outputPtr->CopyInformation( fixedPtr );
  outputPtr->SetRegions( this->m_FixedImageRegion );

  // Moving search region - moving kernel mean.
  outputPtr = this->GetOutput( 3 );
  outputPtr->CopyInformation( movingPtr );
  MovingImageRegionType movingRegion = this->m_MovingImageRegion;
  movingRegion.PadByRadius( this->m_MovingRadius );
  // Make sure the requested region is within the largest possible.
  if ( movingRegion.Crop( movingPtr->GetLargestPossibleRegion() ) )
    {
    outputPtr->SetRegions( movingRegion );
    }
  else
    {
    outputPtr->SetRegions( movingRegion );

    itkExceptionMacro(<< "Moving image requested region is at least partially outside the LargestPossibleRegion.");
    }

  if( !this->m_FixedImageRegionDefined )
    {
    itkExceptionMacro( << "Fixed image Region has not been set." );
    }

  // Moving image of ones.
  outputPtr = this->GetOutput( 4 );
  outputPtr->CopyInformation( movingPtr );

  // FixedPseudoSigmaImage.
  outputPtr = this->GetOutput( 5 );
  outputPtr->CopyInformation( movingPtr );
  outputPtr->SetRegions( this->m_MovingImageRegion );

  // FixedMinusMeanSquared.
  outputPtr = this->GetOutput( 6 );
  outputPtr->CopyInformation( movingPtr );
  outputPtr->SetRegions( this->m_FixedImageRegion );
}

template <class TFixedImage, class TMovingImage,
          class TMetricImage >
void
NormalizedCrossCorrelationMetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
::EnlargeOutputRequestedRegion(DataObject *data)
{
  MetricImagePointerType outputPtr = this->GetOutput( 0 );

  outputPtr->SetRequestedRegionToLargestPossibleRegion();
}

template <class TFixedImage, class TMovingImage,
          class TMetricImage >
void
NormalizedCrossCorrelationMetricImageFilter< TFixedImage, TMovingImage, TMetricImage >
::GenerateHelperImages()
{
  FixedImageConstPointerType  fixedPtr  = this->GetInput( 0 );
  MovingImageConstPointerType movingPtr = this->GetInput( 1 );
  if( !fixedPtr || !movingPtr )
    return;

  // It will screw up all the requested region assumptions.
  if( fixedPtr.GetPointer() == movingPtr.GetPointer() )
    {
    itkExceptionMacro( << "The fixed image and the moving image must be different." );
    }

  if( ! (fixedPtr->GetSpacing() == movingPtr->GetSpacing() ) )
    {
    itkExceptionMacro( << "This metric image filter assumes the moving and fixed image have the same spacing." );
    }

  // This is the m_MovingImageRegion dilated by the radius and cropped by the
  // LargestPossibleRegion.
  MovingImageRegionType movingRequestedRegion = movingPtr->GetRequestedRegion();

  bool movingImageRegionIsSmall = false;
  for( unsigned int i = 0; i < ImageDimension; ++i )
    {
    if( 2 * this->m_MovingRadius[i] + 1 >= this->m_MovingImageRegion.GetSize()[i] )
      {
      movingImageRegionIsSmall = true;
      }
    }
  if( movingImageRegionIsSmall )
    {
    MetricImagePixelType movingMean = NumericTraits< MetricImagePixelType >::Zero;
    typedef ImageRegionConstIterator< MovingImageType > movingImageIteratorType;
    movingImageIteratorType movingIt( movingPtr, movingRequestedRegion );
    for( movingIt.GoToBegin(); !movingIt.IsAtEnd(); ++movingIt )
      {
      movingMean += static_cast< MetricImagePixelType >( movingIt.Get() );
      }
    movingMean /= static_cast< MetricImagePixelType >( movingRequestedRegion.GetNumberOfPixels() );
    MetricImagePointerType movingMeanImg = m_BoxMeanFilter->GetOutput();
    movingMeanImg->SetBufferedRegion( movingRequestedRegion );
    movingMeanImg->Allocate();
    movingMeanImg->FillBuffer( movingMean );
    }
  else
    {
    // Calculate the means.
    m_BoxMeanFilter->SetRadius( this->m_MovingRadius );

    m_BoxMeanFilter->SetNumberOfThreads( this->GetNumberOfThreads() );
    m_BoxMeanFilter->SetInput( movingPtr );
    m_BoxMeanFilter->GetOutput()->SetRequestedRegion( movingRequestedRegion );
    m_BoxMeanFilter->Update();
    }

  MetricImagePixelType fixedMean = NumericTraits< MetricImagePixelType >::Zero;
  typedef ImageRegionConstIterator< FixedImageType > FixedImageIteratorType;
  FixedImageIteratorType fixedIt( fixedPtr, this->m_FixedImageRegion );
  for( fixedIt.GoToBegin(); !fixedIt.IsAtEnd(); ++fixedIt )
    {
    fixedMean += static_cast< MetricImagePixelType >( fixedIt.Get() );
    }
  fixedMean /= static_cast< MetricImagePixelType >( this->m_FixedImageRegion.GetNumberOfPixels() );

  // Calculate the fixed image less the fixed image mean.
  MetricImagePointerType fixedMinusMean = this->GetOutput( 2 );

  ImageRegionIterator< MetricImageType > fixedMinusMeanIt( fixedMinusMean, this->m_FixedImageRegion );
  for( fixedIt.GoToBegin(), fixedMinusMeanIt.GoToBegin();
       !fixedIt.IsAtEnd();
       ++fixedIt, ++fixedMinusMeanIt )
    {
    fixedMinusMeanIt.Set(
      static_cast< MetricImagePixelType >( fixedIt.Get() ) - fixedMean );
    }

  // Calculate the psuedo standard deviations.
  // We have to be careful about the border.

  // Fill an image with ones to signify we are within the image.
  MetricImagePointerType fixedOnes = this->GetOutput( 4 );
  fixedOnes->FillBuffer( NumericTraits< MetricImagePixelType >::One );

  // The value of the fixed pseudo sigma when we are away from the image
  // boundary.
  ImageRegionConstIterator< MetricImageType > fixedMinusMeanConstIt( fixedMinusMean, this->m_FixedImageRegion );
  MetricImagePixelType fixedPseudoSigma = NumericTraits< MetricImagePixelType >::Zero;
  MetricImagePixelType temp;
  for( fixedMinusMeanConstIt.GoToBegin();
       !fixedMinusMeanConstIt.IsAtEnd();
       ++fixedMinusMeanConstIt )
    {
    temp = static_cast< MetricImagePixelType >( fixedMinusMeanConstIt.Get() );
    fixedPseudoSigma += temp * temp;
    }
  fixedPseudoSigma = vcl_sqrt( fixedPseudoSigma );

  // Has the same value everywhere but the border.
  MetricImagePointerType fixedPseudoSigmaImage = this->GetOutput( 5 );
  // Set the value in the interior.
  fixedPseudoSigmaImage->FillBuffer( fixedPseudoSigma );
  // Modify the value at the borders.
  MetricImagePointerType fixedMinusMeanSquared = this->GetOutput( 6 );
  ImageRegionIterator< MetricImageType > fixedMinusMeanSquaredIt( fixedMinusMeanSquared, this->m_FixedImageRegion );
  for( fixedMinusMeanConstIt.GoToBegin(), fixedMinusMeanSquaredIt.GoToBegin();
    !fixedMinusMeanConstIt.IsAtEnd();
    ++fixedMinusMeanConstIt, ++fixedMinusMeanSquaredIt )
    {
    fixedMinusMeanSquaredIt.Set( fixedMinusMeanConstIt.Get() * fixedMinusMeanConstIt.Get() );
    }
  ImageKernelOperator< typename MetricImageType::PixelType, ImageDimension > fixedKernelOperator;
  fixedKernelOperator.SetImageKernel( fixedMinusMeanSquared );
  fixedKernelOperator.CreateToRadius( this->m_FixedRadius );
  NeighborhoodInnerProduct< MetricImageType > innerProduct;
  typedef ConstNeighborhoodIterator< MetricImageType > NeighborhoodIteratorType;
  typedef typename NeighborhoodAlgorithm::
    ImageBoundaryFacesCalculator< MovingImageType > FaceCalculatorType;
  FaceCalculatorType faceCalculator;
  typename FaceCalculatorType::FaceListType faceList = faceCalculator(
    movingPtr, this->m_MovingImageRegion, this->m_MovingRadius );
  typename FaceCalculatorType::FaceListType::iterator fit;
  for( fit = faceList.begin(), ++fit; fit != faceList.end(); ++fit )
    {
    if( (*fit).IsInside( this->m_MovingImageRegion ) )
      {
      ImageRegionIterator< MetricImageType > fixedPseudoSigmaFaceIt( fixedPseudoSigmaImage, *fit );
      NeighborhoodIteratorType nIt( this->m_FixedRadius, fixedOnes, *fit );
      nIt.OverrideBoundaryCondition( &m_BoundaryCondition );
      for( fixedPseudoSigmaFaceIt.GoToBegin(), nIt.GoToBegin();
        !fixedPseudoSigmaFaceIt.IsAtEnd();
        ++fixedPseudoSigmaFaceIt, ++nIt )
        {
        fixedPseudoSigmaFaceIt.Set( vcl_sqrt( innerProduct( nIt, fixedKernelOperator ) ) );
        }
      }
    }

  //Calculate the moving search region less the moving kernel means.
  MetricImagePointerType movingMinusMean = this->GetOutput( 3 );
  ImageRegionIterator< MetricImageType > movingMinusMeanIt( movingMinusMean, movingRequestedRegion );
  ImageRegionConstIterator< MetricImageType > meanIt( m_BoxMeanFilter->GetOutput(), movingRequestedRegion );
  ImageRegionConstIterator< MovingImageType > movingIt( movingPtr, movingRequestedRegion );
  for( movingMinusMeanIt.GoToBegin(), meanIt.GoToBegin(), movingIt.GoToBegin();
       !movingMinusMeanIt.IsAtEnd();
       ++movingMinusMeanIt, ++meanIt, ++movingIt )
    {
    movingMinusMeanIt.Set( movingIt.Get() - meanIt.Get() );
    }

  if( movingImageRegionIsSmall )
    {
    // The cropped region is too small.

    MetricImagePixelType   movingPseudoSigmaVal = NumericTraits< MetricImagePixelType >::Zero;
    MetricImagePointerType movingPseudoSigma = m_BoxPseudoSigmaFilter->GetOutput();
    for( movingMinusMeanIt.GoToBegin(); !movingMinusMeanIt.IsAtEnd(); ++movingMinusMeanIt )
      {
      movingPseudoSigmaVal += movingMinusMeanIt.Get() * movingMinusMeanIt.Get();
      }
    movingPseudoSigmaVal = vcl_sqrt( movingPseudoSigmaVal );
    movingPseudoSigma->SetBufferedRegion( this->m_MovingImageRegion );
    movingPseudoSigma->Allocate();
    movingPseudoSigma->FillBuffer( movingPseudoSigmaVal );
    }
  else
    {
    // Calculate the pseudo sigma in the moving image.
    m_BoxPseudoSigmaFilter->SetNumberOfThreads( this->GetNumberOfThreads() );
    m_BoxPseudoSigmaFilter->SetRadius( this->m_MovingRadius );
    m_BoxPseudoSigmaFilter->SetInput( movingPtr );
    m_BoxPseudoSigmaFilter->GetOutput()->SetRequestedRegion( this->m_MovingImageRegion );
    m_BoxPseudoSigmaFilter->Update();
    }

  MetricImagePointerType denom = this->GetOutput( 1 );
  ImageRegionConstIterator< MetricImageType > fixedPseudoSigmaConstIt( fixedPseudoSigmaImage, this->m_MovingImageRegion );
  ImageRegionConstIterator< MetricImageType > movingPseudoSigmaConstIt( m_BoxPseudoSigmaFilter->GetOutput(), this->m_MovingImageRegion );
  ImageRegionIterator< MetricImageType > denomIt( denom, this->m_MovingImageRegion );
  for( fixedPseudoSigmaConstIt.GoToBegin(), denomIt.GoToBegin(), movingPseudoSigmaConstIt.GoToBegin();
       !fixedPseudoSigmaConstIt.IsAtEnd();
       ++fixedPseudoSigmaConstIt, ++denomIt, ++movingPseudoSigmaConstIt )
    {
    denomIt.Set( movingPseudoSigmaConstIt.Get() * fixedPseudoSigmaConstIt.Get() );
    }
}


} // end namespace BlockMatching
} // end namespace itk

#endif
