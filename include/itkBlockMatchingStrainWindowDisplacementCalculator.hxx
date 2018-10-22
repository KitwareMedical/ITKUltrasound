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
#ifndef itkBlockMatchingStrainWindowDisplacementCalculator_hxx
#define itkBlockMatchingStrainWindowDisplacementCalculator_hxx

#include "itkBlockMatchingStrainWindowDisplacementCalculator.h"
#include "itkBlockMatchingMaximumPixelDisplacementCalculator.h"

#include "itkStrainImageFilter.h"

#include "itkFixedArray.h"
#include "itkImageLinearConstIteratorWithIndex.h"
#include "itkImageLinearIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

#include <set>

namespace itk
{
namespace BlockMatching
{

template<class TMetricImage, class TDisplacementImage, class TStrainValueType>
StrainWindowDisplacementCalculator<TMetricImage, TDisplacementImage, TStrainValueType>
::StrainWindowDisplacementCalculator():
  m_CurrentIteration( 0 ),
  m_MaximumIterations( 1 )
{
  this->m_CacheMetricImage = true;

  typedef itk::StrainImageFilter<TDisplacementImage, TStrainValueType, TStrainValueType>
  itkStrainImageFilterType;

  // Sensible default.
  m_StrainImageFilter = itkStrainImageFilterType::New().GetPointer();

  // Sensible default.
  m_DisplacementCalculator =
    MaximumPixelDisplacementCalculator<TMetricImage,
                                       TDisplacementImage>::New();

  m_MaximumAbsStrain.Fill( NumericTraits<TStrainValueType>::max() );

  m_Mask = MaskType::New();
  m_NthElementAdaptor = NthElementAdaptorType::New();
  m_AbsFilter = AbsFilterType::New();
  m_AbsFilter->SetInput( m_NthElementAdaptor );
  m_BoxMeanFilter = BoxMeanImageFilterType::New();
  m_BoxMeanFilter->SetRadius( 1 );
  m_BoxMeanFilter->SetInput( m_AbsFilter->GetOutput() );
}


template<class TMetricImage, class TDisplacementImage, class TStrainValueType>
void
StrainWindowDisplacementCalculator<TMetricImage, TDisplacementImage, TStrainValueType>
::SetMetricImagePixel( const PointType & point,
                       const IndexType& index,
                       MetricImageType * image )
{
  Superclass::SetMetricImagePixel( point, index, image );
  this->m_DisplacementCalculator->SetMetricImagePixel( point, index, image );
}


template<class TMetricImage, class TDisplacementImage, class TStrainValueType>
unsigned long long
StrainWindowDisplacementCalculator<TMetricImage, TDisplacementImage, TStrainValueType>
::GenerateMask( const RegionType & region )
{
  // Create the strain image.
  m_StrainImageFilter->SetInput( this->m_DisplacementImage );
  typename StrainImageType::Pointer strain = m_StrainImageFilter->GetOutput();

  // Determine which displacements we need to displace.
  if( m_Mask->GetBufferedRegion() != region )
    {
    m_Mask->SetRegions( region );
    m_Mask->Allocate();
    }
  m_Mask->FillBuffer( false );
  m_NthElementAdaptor->SetImage( strain );
  m_BoxMeanFilter->GetOutput()->SetRegions( region );
  unsigned long long count = 0;
  ImageRegionIterator< MaskType >   maskIt( m_Mask, region );

  //! @todo multithread me
  constexpr unsigned int i = 0;
  //for( unsigned int i = 0; i < m_MaximumAbsStrain.GetNumberOfComponents(); ++i )
    //{
    m_NthElementAdaptor->SelectNthElement( i );
    // this will be correct after this is multi-threaded and a fixer-uper for
    // now
    m_BoxMeanFilter->SetNumberOfWorkUnits( 1 );
    m_BoxMeanFilter->Update();
    ImageRegionIterator< MetricImageType > strainCompIt( m_BoxMeanFilter->GetOutput(), region );

    for( maskIt.GoToBegin(), strainCompIt.GoToBegin();
         !maskIt.IsAtEnd();
         ++maskIt, ++strainCompIt )
      {
      if( !maskIt.Value() )
        {
        if( strainCompIt.Value() > m_MaximumAbsStrain[i] )
          {
          maskIt.Set( true );
          ++count;
          }
        }
      }
    //}

  return count;
}


template<class TMetricImage, class TDisplacementImage, class TStrainValueType>
void
StrainWindowDisplacementCalculator<TMetricImage, TDisplacementImage, TStrainValueType>
::ReplaceDisplacements( const RegionType & region )
{
  // The lines in every direction that need extrapolation.
  typedef std::set<long> NeedsExtrapolationSetType;
  NeedsExtrapolationSetType::const_iterator setIt;
  typedef itk::FixedArray<NeedsExtrapolationSetType, ImageDimension> NeedsExtrapolationType;
  NeedsExtrapolationType needsExtrapolation;
  IndexType              needsExtrapolationIndex;

  typedef itk::ImageLinearIteratorWithIndex< DisplacementImageType > DisplacementIteratorType;
  DisplacementIteratorType dispIt( this->m_DisplacementImage,
                                   region );

  typedef itk::ImageLinearConstIteratorWithIndex< DisplacementImageType >
    DisplacementConstIteratorType;
  DisplacementConstIteratorType extrapolationDispIt( this->m_DisplacementImage,
    region );
  typedef itk::ImageLinearConstIteratorWithIndex< MaskType > MaskIteratorType;
  MaskIteratorType maskIt( this->m_Mask, region );
  TStrainValueType normalStrain;
  typedef typename DisplacementImageType::PixelType DisplacementVectorType;
  DisplacementVectorType            displacementStart;
  DisplacementVectorType            displacementEnd;
  double                            spacing;
  double                            slope;
  long                              offset;
  IndexType                         offsetIndex;
  typename StrainImageType::Pointer strain = m_StrainImageFilter->GetOutput();
  for( unsigned int dim = 0; dim < ImageDimension; ++dim )
    {
    dispIt.SetDirection( dim );
    maskIt.SetDirection( dim );
    for( dispIt.GoToBegin(), maskIt.GoToBegin();
         !dispIt.IsAtEnd();
         dispIt.NextLine(), maskIt.NextLine() )
      {
      dispIt.GoToBeginOfLine();
      maskIt.GoToBeginOfLine();
      if( maskIt.Get() )
        {
        while( !maskIt.IsAtEndOfLine() && maskIt.Get() )
          {
          ++maskIt;
          }

        if( maskIt.IsAtEndOfLine() )
          {
          needsExtrapolation[dim].insert( maskIt.GetIndex()[dim] );
          continue;
          }
        else
          {
          // The normal strain in this direction.
          normalStrain = strain->GetPixel( maskIt.GetIndex() ) ( dim, dim );
          dispIt.SetIndex( maskIt.GetIndex() );
          displacementEnd = dispIt.Get();
          spacing = this->m_DisplacementImage->GetSpacing()[dim];
          slope = normalStrain * spacing;
          --dispIt;
          while( !dispIt.IsAtReverseEndOfLine() )
            {
            displacementEnd[dim] -= slope;
            dispIt.Set( displacementEnd );
            --dispIt;
            }
          }
        dispIt.SetIndex( maskIt.GetIndex() );
        }
      while( !dispIt.IsAtEndOfLine() )
        {
        if( maskIt.Get() )
          {
          --dispIt;
          displacementStart = dispIt.Get();
          while( !maskIt.IsAtEndOfLine() && maskIt.Get() )
            {
            ++maskIt;
            }

          if( maskIt.IsAtEndOfLine() )
            {
            // The normal strain in this direction.
            normalStrain = strain->GetPixel( dispIt.GetIndex() ) ( dim, dim );
            spacing = this->m_DisplacementImage->GetSpacing()[dim];
            slope = normalStrain * spacing;
            ++dispIt;
            while( !dispIt.IsAtEndOfLine() )
              {
              displacementStart[dim] += slope;
              dispIt.Set( displacementStart );
              ++dispIt;
              }

            continue;
            }
          else
            {
            displacementEnd = this->m_DisplacementImage->GetPixel( maskIt.GetIndex() );
            slope = ( displacementEnd[dim] - displacementStart[dim] )
              / ( maskIt.GetIndex()[dim] - dispIt.GetIndex()[dim] );
            ++dispIt;
            maskIt.SetIndex( dispIt.GetIndex() );
            while( maskIt.Get() )
              {
              displacementStart[dim] += slope;
              dispIt.Set( displacementStart );
              ++maskIt;
              ++dispIt;
              }
            }
          }
        ++dispIt;
        ++maskIt;
        }
      }
    for( setIt = needsExtrapolation[dim].begin(); setIt != needsExtrapolation[dim].end(); ++setIt )
      {
      needsExtrapolationIndex = region.GetIndex();
      offsetIndex             = region.GetIndex();
      needsExtrapolationIndex[dim] = *setIt;
      dispIt.SetIndex( needsExtrapolationIndex );
      offset = 1;
      while( offset < static_cast<long>( region.GetSize()[dim] ) )
        {
        offsetIndex[dim] = needsExtrapolationIndex[dim] + offset;
        if( region.IsInside( offsetIndex ) &&
            needsExtrapolation[dim].find( needsExtrapolationIndex[dim] + offset ) == needsExtrapolation[dim].end() )
          {
          extrapolationDispIt.SetIndex( offsetIndex );
          while( !dispIt.IsAtEndOfLine() )
            {
            dispIt.Set( extrapolationDispIt.Get() );
            ++dispIt;
            ++extrapolationDispIt;
            }

          break;
          }
        else
          {
          offsetIndex[dim] = needsExtrapolationIndex[dim] - offset;
          if( region.IsInside( offsetIndex ) &&
              needsExtrapolation[dim].find( needsExtrapolationIndex[dim] - offset ) == needsExtrapolation[dim].end() )
            {
            extrapolationDispIt.SetIndex( offsetIndex );
            while( !dispIt.IsAtEndOfLine() )
              {
              dispIt.Set( extrapolationDispIt.Get() );
              ++dispIt;
              ++extrapolationDispIt;
              }

            break;
            }
          }
        ++offset;
        }

      if( offset == static_cast<long>( region.GetSize()[dim] ) )
        {
        itkExceptionMacro(<< "No values inside the strain window were found.");
        }
      }
    } // end for dimension
  this->m_DisplacementImage->Modified();
}


template<class TMetricImage, class TDisplacementImage, class TStrainValueType>
void
StrainWindowDisplacementCalculator<TMetricImage, TDisplacementImage, TStrainValueType>
::Compute()
{
  this->m_DisplacementCalculator->Compute();

  this->InvokeEvent( StartEvent() );

  typename DisplacementImageType::RegionType region = this->m_DisplacementImage->GetBufferedRegion();

  unsigned long long valuesToReplace = this->GenerateMask( region );
  this->ReplaceDisplacements( region );

  this->m_CurrentIteration = 1;
  unsigned long long nextValuesToReplace = itk::NumericTraits<unsigned long long>::max();
  while( this->m_CurrentIteration < this->m_MaximumIterations &&
         nextValuesToReplace != valuesToReplace &&
         nextValuesToReplace != 0 )
    {
    valuesToReplace = nextValuesToReplace;
    nextValuesToReplace = this->GenerateMask( region );
    this->ReplaceDisplacements( region );
    ++this->m_CurrentIteration;
    }

  this->InvokeEvent( EndEvent() );
}

} // end namespace BlockMatching
} // end namespace itk

#endif
