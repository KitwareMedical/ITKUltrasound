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
#ifndef itkBlockMatchingMetricImageToDisplacementCalculator_hxx
#define itkBlockMatchingMetricImageToDisplacementCalculator_hxx

#include "itkBlockMatchingMetricImageToDisplacementCalculator.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

namespace itk
{
namespace BlockMatching
{

template < typename TMetricImage, typename TDisplacementImage >
MetricImageToDisplacementCalculator< TMetricImage, TDisplacementImage >
::MetricImageToDisplacementCalculator():
  m_CacheMetricImage( false ),
  m_RegionsDefined( false )
{
  m_MetricImageImage  = ITK_NULLPTR;
  m_DisplacementImage = ITK_NULLPTR;
  m_MetricImageDuplicator = MetricImageDuplicatorType::New();

  m_Threader = MultiThreader::New();
  m_NumberOfThreads = m_Threader->GetNumberOfThreads();
}


template < typename TMetricImage, typename TDisplacementImage >
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


template < typename TMetricImage, typename TDisplacementImage >
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


template < typename TMetricImage, typename TDisplacementImage >
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


template < typename TMetricImage, typename TDisplacementImage >
ITK_THREAD_RETURN_TYPE
MetricImageToDisplacementCalculator< TMetricImage, TDisplacementImage >
::ThreaderCallback( void *arg )
{
  ThreadStruct *str;

  const ThreadIdType threadId = ((MultiThreader::ThreadInfoStruct *)(arg))->ThreadID;
  const ThreadIdType threadCount = ((MultiThreader::ThreadInfoStruct *)(arg))->NumberOfThreads;

  str = (ThreadStruct *)(((MultiThreader::ThreadInfoStruct *)(arg))->UserData);

  // execute the actual method with appropriate output region
  // first find out how many pieces extent can be split into.
  RegionType splitRegion;
  const ThreadIdType total = str->self->SplitRequestedRegion( threadId,
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


template < typename TMetricImage, typename TDisplacementImage >
unsigned int
MetricImageToDisplacementCalculator< TMetricImage, TDisplacementImage >
::SplitRequestedRegion( unsigned int i, unsigned int num,
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
  const ThreadIdType valuesPerThread = Math::Ceil< ThreadIdType >( range/(double)num );
  const ThreadIdType maxThreadIdUsed = Math::Ceil< ThreadIdType >( range/(double)valuesPerThread ) - 1;

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
