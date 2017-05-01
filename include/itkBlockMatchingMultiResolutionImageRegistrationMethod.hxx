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
#ifndef itkBlockMatchingMultiResolutionImageRegistrationMethod_hxx
#define itkBlockMatchingMultiResolutionImageRegistrationMethod_hxx

#include "itkBlockMatchingMultiResolutionImageRegistrationMethod.h"

#include "itkRecursiveMultiResolutionPyramidImageFilter.h"

namespace itk
{

namespace BlockMatching
{

template <class TFixedImage, class TMovingImage,
          class TMetricImage, class TDisplacementImage, class TCoordRep>
MultiResolutionImageRegistrationMethod<TFixedImage, TMovingImage,
                                       TMetricImage, TDisplacementImage, TCoordRep>
::MultiResolutionImageRegistrationMethod() :
  m_FixedImage( ITK_NULLPTR ),
  m_MovingImage( ITK_NULLPTR ),
  m_NumberOfLevels( 1 ),
  m_CurrentLevel( 0 ),
  m_Stop( false ),
  m_ScheduleSpecified( false ),
  m_NumberOfLevelsSpecified( false ),
  m_ImageRegistrationMethod( ITK_NULLPTR ),
  m_BlockRadiusCalculator( ITK_NULLPTR ),
  m_SearchRegionImageSource( ITK_NULLPTR )
{
  m_FixedImagePyramid  = RecursiveMultiResolutionPyramidImageFilter<
      FixedImageType, FixedImageType>::New();
  m_MovingImagePyramid = RecursiveMultiResolutionPyramidImageFilter<
      MovingImageType, MovingImageType>::New();
}


template <class TFixedImage, class TMovingImage,
          class TMetricImage, class TDisplacementImage, class TCoordRep>
unsigned long
MultiResolutionImageRegistrationMethod<TFixedImage, TMovingImage,
                                       TMetricImage, TDisplacementImage, TCoordRep>
::GetMTime() const
{
  unsigned long mtime = Superclass::GetMTime();
  unsigned long m;

  if( m_FixedImage )
    {
    m = m_FixedImage->GetMTime();
    mtime = (m > mtime ? m : mtime);
    }

  if( m_MovingImage )
    {
    m = m_MovingImage->GetMTime();
    mtime = (m > mtime ? m : mtime);
    }

  return mtime;
}


template <class TFixedImage, class TMovingImage,
          class TMetricImage, class TDisplacementImage, class TCoordRep>
void
MultiResolutionImageRegistrationMethod<TFixedImage, TMovingImage,
                                       TMetricImage, TDisplacementImage, TCoordRep>
::StopRegistration( void )
{
  m_Stop = true;
}


template <class TFixedImage, class TMovingImage,
          class TMetricImage, class TDisplacementImage, class TCoordRep>
void
MultiResolutionImageRegistrationMethod<TFixedImage, TMovingImage,
                                       TMetricImage, TDisplacementImage, TCoordRep>
::Initialize()
{
  // Sanity checks
  if( !m_ImageRegistrationMethod )
    {
    itkExceptionMacro(<< "ImageRegistrationMethod is not present");
    }

  this->PrepareBlockRadiusCalculator();
  this->PrepareSearchRegionImageSource();
}


/**
 * Set the schedules for the fixed and moving image pyramid.  Taken form
 * itkMultiResolutionImageRegistrationMethod.
 */
template <class TFixedImage, class TMovingImage,
          class TMetricImage, class TDisplacementImage, class TCoordRep>
void
MultiResolutionImageRegistrationMethod<TFixedImage, TMovingImage,
                                       TMetricImage, TDisplacementImage, TCoordRep>
::SetSchedules( const ScheduleType & fixedImagePyramidSchedule, const ScheduleType & movingImagePyramidSchedule )
{
  if( m_NumberOfLevelsSpecified )
    {
    itkExceptionMacro( "SetSchedules should not be used "
                       << "if numberOfLevelves are specified using SetNumberOfLevels" );
    }
  m_FixedImagePyramidSchedule = fixedImagePyramidSchedule;
  m_MovingImagePyramidSchedule = movingImagePyramidSchedule;
  m_ScheduleSpecified = true;

  // Set the number of levels based on the pyramid schedule specified
  if( m_FixedImagePyramidSchedule.rows() !=
      m_MovingImagePyramidSchedule.rows() )
    {
    itkExceptionMacro("The specified schedules contain unequal number of levels");
    }
  else
    {
    m_NumberOfLevels = m_FixedImagePyramidSchedule.rows();
    }

  m_BlockRadiusCalculator->SetPyramidSchedule( fixedImagePyramidSchedule );
  this->Modified();
}


template <class TFixedImage, class TMovingImage,
          class TMetricImage, class TDisplacementImage, class TCoordRep>
void
MultiResolutionImageRegistrationMethod<TFixedImage, TMovingImage,
                                       TMetricImage, TDisplacementImage, TCoordRep>
::SetNumberOfLevels( unsigned long numberOfLevels )
{
  if( m_ScheduleSpecified )
    {
    itkExceptionMacro( "SetNumberOfLevels should not be used "
                       << "if schedules have been specified using SetSchedules method " );
    }

  m_NumberOfLevels = numberOfLevels;
  m_NumberOfLevelsSpecified = true;
  this->Modified();
}


template <class TFixedImage, class TMovingImage,
          class TMetricImage, class TDisplacementImage, class TCoordRep>
void
MultiResolutionImageRegistrationMethod<TFixedImage, TMovingImage,
                                       TMetricImage, TDisplacementImage, TCoordRep>
::GenerateOutputInformation()
{
  if( !m_FixedImage )
    {
    itkExceptionMacro(<< "FixedImage is not present");
    }
  // So that SearchRegionImageSource has the right information to work with,
  // possibly others too.
  m_FixedImage->UpdateOutputInformation();

  if( !m_MovingImage )
    {
    itkExceptionMacro(<< "MovingImage is not present");
    }
  m_MovingImage->UpdateOutputInformation();

  m_SearchRegionImageSource->SetFixedImage( m_FixedImage );
  m_SearchRegionImageSource->SetMovingImage( m_MovingImage );
  m_SearchRegionImageSource->SetCurrentLevel( m_NumberOfLevels - 1 );
  m_SearchRegionImageSource->SetFixedBlockRadius( m_BlockRadiusCalculator->Compute( m_NumberOfLevels - 1 ) );
  m_SearchRegionImageSource->UpdateOutputInformation();

  TDisplacementImage* output = this->GetOutput( 0 );
  if( !output )
    {
    return;
    }

  output->CopyInformation( m_SearchRegionImageSource->GetOutput( 0 ) );
}


template <class TFixedImage, class TMovingImage,
          class TMetricImage, class TDisplacementImage, class TCoordRep>
void
MultiResolutionImageRegistrationMethod<TFixedImage, TMovingImage,
                                       TMetricImage, TDisplacementImage, TCoordRep>
::PreparePyramids()
{
  // Sanity checks
  if( !m_FixedImagePyramid )
    {
    itkExceptionMacro(<< "Fixed image pyramid is not present");
    }

  if( !m_MovingImagePyramid )
    {
    itkExceptionMacro(<< "Moving image pyramid is not present");
    }

  // Setup the fixed and moving image pyramid
  if( m_NumberOfLevelsSpecified )
    {
    m_FixedImagePyramid->SetNumberOfLevels( m_NumberOfLevels );
    m_MovingImagePyramid->SetNumberOfLevels( m_NumberOfLevels );
    }

  if( m_ScheduleSpecified )
    {
    m_FixedImagePyramid->SetNumberOfLevels( m_FixedImagePyramidSchedule.rows() );
    m_FixedImagePyramid->SetSchedule( m_FixedImagePyramidSchedule );

    m_MovingImagePyramid->SetNumberOfLevels( m_MovingImagePyramidSchedule.rows() );
    m_MovingImagePyramid->SetSchedule( m_MovingImagePyramidSchedule );
    }

  m_FixedImagePyramid->SetInput( m_FixedImage );
  m_FixedImagePyramid->UpdateLargestPossibleRegion();

  m_MovingImagePyramid->SetInput( m_MovingImage );
  m_MovingImagePyramid->UpdateLargestPossibleRegion();
}


template <class TFixedImage, class TMovingImage,
          class TMetricImage, class TDisplacementImage, class TCoordRep>
void
MultiResolutionImageRegistrationMethod<TFixedImage, TMovingImage,
                                       TMetricImage, TDisplacementImage, TCoordRep>
::PrepareBlockRadiusCalculator()
{
  // Sanity checks
  if( !m_BlockRadiusCalculator )
    {
    itkExceptionMacro(<< "BlockRadiusCalculator is not present");
    }

  m_BlockRadiusCalculator->SetFixedImage( m_FixedImage );
  m_BlockRadiusCalculator->SetPyramidSchedule( m_FixedImagePyramid->GetSchedule() );
}


template <class TFixedImage, class TMovingImage,
          class TMetricImage, class TDisplacementImage, class TCoordRep>
void
MultiResolutionImageRegistrationMethod<TFixedImage, TMovingImage,
                                       TMetricImage, TDisplacementImage, TCoordRep>
::PrepareSearchRegionImageSource()
{
  // Sanity checks
  if( !m_SearchRegionImageSource )
    {
    itkExceptionMacro(<< "SearchRegionImageSource is not present");
    }

  m_SearchRegionImageSource->SetPyramidSchedule( m_FixedImagePyramid->GetSchedule() );
  m_ImageRegistrationMethod->SetInput( m_SearchRegionImageSource->GetOutput() );
}


template <class TFixedImage, class TMovingImage,
          class TMetricImage, class TDisplacementImage, class TCoordRep>
void
MultiResolutionImageRegistrationMethod<TFixedImage, TMovingImage,
                                       TMetricImage, TDisplacementImage, TCoordRep>
::GenerateData()
{
  m_Stop = false;

  this->PreparePyramids();
  for( m_CurrentLevel = 0; m_CurrentLevel < m_NumberOfLevels; ++m_CurrentLevel )
    {
    this->Initialize();

    // Check if there has been a stop request
    if( m_Stop )
      {
      break;
      }

    m_SearchRegionImageSource->SetFixedImage( m_FixedImagePyramid->GetOutput( m_CurrentLevel ) );
    m_SearchRegionImageSource->SetMovingImage( m_MovingImagePyramid->GetOutput( m_CurrentLevel ) );
    m_SearchRegionImageSource->SetCurrentLevel( m_CurrentLevel );
    m_SearchRegionImageSource->SetFixedBlockRadius( m_BlockRadiusCalculator->Compute( m_CurrentLevel ) );
    m_ImageRegistrationMethod->SetRadius( m_BlockRadiusCalculator->Compute( m_CurrentLevel ) );
    m_ImageRegistrationMethod->SetFixedImage(  m_FixedImagePyramid->GetOutput( m_CurrentLevel ) );
    m_ImageRegistrationMethod->SetMovingImage( m_MovingImagePyramid->GetOutput( m_CurrentLevel ) );

    // Invoke an iteration event.
    // This allows a UI to reset any of the components between
    // resolution level.
    this->InvokeEvent( IterationEvent() );

    m_ImageRegistrationMethod->UpdateLargestPossibleRegion();
    m_SearchRegionImageSource->SetPreviousDisplacements( m_ImageRegistrationMethod->GetOutput() );
    }

  this->GraftOutput( m_ImageRegistrationMethod->GetOutput() );
}

} // end namespace itk
} // end namespace BlockMatching

#endif
