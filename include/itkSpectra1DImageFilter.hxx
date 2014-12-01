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
#ifndef __itkSpectra1DImageFilter_hxx
#define __itkSpectra1DImageFilter_hxx

#include "itkSpectra1DImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkMetaDataObject.h"

namespace itk
{

template< typename TInputImage, typename TSupportWindowImage, typename TOutputImage >
Spectra1DImageFilter< TInputImage, TSupportWindowImage, TOutputImage >
::Spectra1DImageFilter()
{
  this->AddRequiredInputName( "InputImage" );
  this->SetPrimaryInputName( "InputImage" );
  this->AddRequiredInputName( "SupportWindowImage" );
}


template< typename TInputImage, typename TSupportWindowImage, typename TOutputImage >
void
Spectra1DImageFilter< TInputImage, TSupportWindowImage, TOutputImage >
::ThreadedGenerateData( const OutputImageRegionType & outputRegionForThread, ThreadIdType itkNotUsed( threadId ) )
{
  //OutputImageType * output = this->GetOutput();
  //const InputImageType * input = this->GetInput();

  //const OutputImageRegionType outputLargestRegion = output->GetLargestPossibleRegion();
  //typedef typename OutputImageType::IndexType IndexType;
  //const IndexType largestIndexStart = outputLargestRegion.GetIndex();
  //const IndexType largestIndexStop = largestIndexStart + outputLargestRegion.GetSize();

  //typedef ImageRegionConstIteratorWithIndex< InputImageType > InputIteratorType;
  //InputIteratorType inputIt( input, outputRegionForThread );
  //typedef ImageRegionIterator< OutputImageType > OutputIteratorType;
  //OutputIteratorType outputIt( output, outputRegionForThread );
}

} // end namespace itk

#endif
