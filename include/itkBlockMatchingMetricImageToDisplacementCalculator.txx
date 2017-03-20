#ifndef __itkBlockMatchingMetricImageToDisplacementCalculator_txx
#define __itkBlockMatchingMetricImageToDisplacementCalculator_txx

#include "itkBlockMatchingMetricImageToDisplacementCalculator.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

namespace itk
{
namespace BlockMatching
{

template < class TMetricImage, class TDisplacementImage >
MetricImageToDisplacementCalculator< TMetricImage, TDisplacementImage >
::MetricImageToDisplacementCalculator():
  m_CacheMetricImage( false ),
  m_RegionsDefined( false )
{
  m_MetricImageImage  = NULL;
  m_DisplacementImage = NULL;
  m_MetricImageDuplicator = MetricImageDuplicatorType::New();

  m_Threader = MultiThreader::New();
  m_NumberOfThreads = m_Threader->GetNumberOfThreads();
}


template < class TMetricImage, class TDisplacementImage >
void
MetricImageToDisplacementCalculator< TMetricImage, TDisplacementImage >
::SetMetricImagePixel( const PointType& point, const IndexType& index,
                       MetricImageType * metricImage )
{
  if( m_CacheMetricImage )
    {
    // Copy the image contents to a new image.
    m_MetricImageDuplicator->SetInputImage( metricImage );
    m_MetricImageDuplicator->Update();
    m_MetricImageImage->SetPixel ( index, m_MetricImageDuplicator->GetOutput() );
    m_CenterPointsImage->SetPixel( index, point );
    }
}


template < class TMetricImage, class TDisplacementImage >
void
MetricImageToDisplacementCalculator< TMetricImage, TDisplacementImage >
::SetDisplacementImage( DisplacementImageType * image )
{
  if( this->m_DisplacementImage.GetPointer() != image )
    {
    this->m_DisplacementImage = image;
    this->Modified();
    }

  if( m_CacheMetricImage && ( !m_MetricImageImage.GetPointer() ||
      m_MetricImageImage->GetLargestPossibleRegion() != image->GetLargestPossibleRegion() ))
    {
    // We create a new image and reallocate here every time because if we try to
    // keep the same image reallocation appears to do something funny with the
    // SmartPointer pixels.
    m_MetricImageImage = MetricImageImageType::New();
    m_CenterPointsImage = CenterPointsImageType::New();
    m_MetricImageImage->CopyInformation( image );
    m_MetricImageImage->SetRegions( image->GetLargestPossibleRegion() );
    m_MetricImageImage->Allocate();
    m_CenterPointsImage->CopyInformation( image );
    m_CenterPointsImage->SetRegions( image->GetLargestPossibleRegion() );
    m_CenterPointsImage->Allocate();
    this->Modified();
    }
}


template < class TMetricImage, class TDisplacementImage >
void
MetricImageToDisplacementCalculator< TMetricImage, TDisplacementImage >
::ApplyThreadFunctor( ThreadFunctor& functor )
{
  ThreadStruct str;
  str.self = this;
  str.functor = &functor;
  this->m_Threader->SetSingleMethod(
    this->ThreaderCallback, &str );
  this->m_Threader->SingleMethodExecute();
}


template < class TMetricImage, class TDisplacementImage >
ITK_THREAD_RETURN_TYPE
MetricImageToDisplacementCalculator< TMetricImage, TDisplacementImage >
::ThreaderCallback( void *arg )
{
  ThreadStruct *str;
  int total, threadId, threadCount;

  threadId = ((MultiThreader::ThreadInfoStruct *)(arg))->ThreadID;
  threadCount = ((MultiThreader::ThreadInfoStruct *)(arg))->NumberOfThreads;

  str = (ThreadStruct *)(((MultiThreader::ThreadInfoStruct *)(arg))->UserData);

  // execute the actual method with appropriate output region
  // first find out how many pieces extent can be split into.
  RegionType splitRegion;
  total = str->self->SplitRequestedRegion( threadId,
    threadCount, splitRegion);

  if (threadId < total)
    {
    return (*(str->functor))( str->self, splitRegion, threadId );
    }
  // else
  //   {
  //   otherwise don't use this thread. Sometimes the threads dont
  //   break up very well and it is just as efficient to leave a 
  //   few threads idle.
  //   }
  
  return ITK_THREAD_RETURN_VALUE;
}


template < class TMetricImage, class TDisplacementImage >
int
MetricImageToDisplacementCalculator< TMetricImage, TDisplacementImage >
::SplitRequestedRegion( int i, int num,
                        RegionType& splitRegion)
{
  RegionType                                   regionToSplit = this->m_DisplacementImage->GetRequestedRegion();
  const typename TDisplacementImage::SizeType& requestedRegionSize
    = regionToSplit.GetSize();

  int                                    splitAxis;
  typename TDisplacementImage::IndexType splitIndex;
  typename TDisplacementImage::SizeType  splitSize;

  // Initialize the splitRegion to the output requested region
  splitRegion = regionToSplit;
  splitIndex = splitRegion.GetIndex();
  splitSize = splitRegion.GetSize();

  // split on the outermost dimension available
  splitAxis = regionToSplit.GetImageDimension() - 1;
  while (requestedRegionSize[splitAxis] == 1)
    {
    --splitAxis;
    if (splitAxis < 0)
      {   // cannot split
      itkDebugMacro("  Cannot Split");
      return 1;
      }
    }

  // determine the actual number of pieces that will be generated
  typename TDisplacementImage::SizeType::SizeValueType range = requestedRegionSize[splitAxis];
  int                                                  valuesPerThread = Math::Ceil<int>(range/(double)num);
  int                                                  maxThreadIdUsed =
    Math::Ceil<int>(range/(double)valuesPerThread) - 1;

  // Split the region
  if (i < maxThreadIdUsed)
    {
    splitIndex[splitAxis] += i*valuesPerThread;
    splitSize[splitAxis] = valuesPerThread;
    }
  if (i == maxThreadIdUsed)
    {
    splitIndex[splitAxis] += i*valuesPerThread;
    // last thread needs to process the "rest" dimension being split
    splitSize[splitAxis] = splitSize[splitAxis] - i*valuesPerThread;
    }

  // set the split region ivars
  splitRegion.SetIndex( splitIndex );
  splitRegion.SetSize( splitSize );

  itkDebugMacro("  Split Piece: " << splitRegion );

  return maxThreadIdUsed + 1;
}
} // end namespace BlockMatching
} // end namespace itk

#endif
