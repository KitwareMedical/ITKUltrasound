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
#ifndef itkUltrasoundImageFileReader_hxx
#define itkUltrasoundImageFileReader_hxx

#include "itkUltrasoundImageFileReader.h"

#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"
#include "itkEuler3DTransform.h"

#include "itkCurvilinearArraySpecialCoordinatesImage.h"
#include "itkSliceSeriesSpecialCoordinatesImage.h"

#ifdef ITK_HAS_GCC_PRAGMA_DIAG_PUSHPOP
  ITK_GCC_PRAGMA_DIAG_PUSH()
#endif
ITK_GCC_PRAGMA_DIAG(ignored "-Wattributes")
namespace
{

double
GetMetaDataAsDouble(const itk::MetaDataDictionary & dictionary, const std::string & key)
{
  std::string valueAsStr;
  itk::ExposeMetaData< std::string >( dictionary, key, valueAsStr );
  std::istringstream istrm( valueAsStr );
  double value;
  istrm >> value;
  return value;
}


template< typename TImage >
class UltrasoundImageFileReaderDispatch
{
public:
  typedef TImage ImageType;
  static void ExtractMetaData( TImage * image ) {};
};


template< typename TPixel, unsigned int VDimension >
class UltrasoundImageFileReaderDispatch< itk::CurvilinearArraySpecialCoordinatesImage< TPixel, VDimension > >
{
public:
  typedef itk::CurvilinearArraySpecialCoordinatesImage< TPixel, VDimension > ImageType;
  static void ExtractMetaData( ImageType * image )
  {
    const itk::MetaDataDictionary & dictionary = image->GetMetaDataDictionary();
    if( dictionary.HasKey( "LateralAngularSeparation" ) )
      {
      image->SetLateralAngularSeparation( GetMetaDataAsDouble( dictionary, "LateralAngularSeparation" ) );
      }
    if( dictionary.HasKey( "FirstSampleDistance" ) )
      {
      image->SetFirstSampleDistance( GetMetaDataAsDouble( dictionary, "FirstSampleDistance" ) );
      }
    if( dictionary.HasKey( "RadiusSampleSize" ) )
      {
      image->SetRadiusSampleSize( GetMetaDataAsDouble( dictionary, "RadiusSampleSize" ) );
      }
  }
};


template< typename TPixel, unsigned int VDimension, typename TParametersValue >
class UltrasoundImageFileReaderDispatch< itk::SliceSeriesSpecialCoordinatesImage< itk::Image< TPixel, VDimension - 1 >, itk::Euler3DTransform< TParametersValue >, TPixel, VDimension > >
{
public:
  static const unsigned int Dimension = VDimension;
  static const unsigned int SliceDimension = VDimension -1;

  typedef TPixel                                                                   PixelType;
  typedef TParametersValue                                                         ParametersValueType;
  typedef itk::Image< PixelType, SliceDimension >                                  SliceImageType;
  typedef itk::Euler3DTransform< ParametersValueType >                             TransformType;
  typedef itk::SliceSeriesSpecialCoordinatesImage< SliceImageType, TransformType,
    PixelType, Dimension > ImageType;
  static void ExtractMetaData( ImageType * image )
  {
    const itk::MetaDataDictionary & dictionary = image->GetMetaDataDictionary();
    // Currently only an Image SliceType is supported.
    if( !dictionary.HasKey( "SliceType" ) )
      {
      return;
      }
    else
      {
      std::string sliceType;
      itk::ExposeMetaData< std::string >( dictionary, "SliceType", sliceType );
      if( sliceType != "Image" )
        {
        return;
        }
      }

    typename SliceImageType::Pointer sliceImage = SliceImageType::New();

    if( dictionary.HasKey( "SliceSpacing" ) )
      {
      typedef itk::Array< double > SliceSpacingArrayType;
      SliceSpacingArrayType sliceSpacingArray( SliceDimension );
      itk::ExposeMetaData< SliceSpacingArrayType >( dictionary, "SliceSpacing", sliceSpacingArray );
      typename SliceImageType::SpacingType sliceSpacing;
      for( unsigned int ii = 0; ii < SliceDimension; ++ii )
        {
        sliceSpacing[ii] = sliceSpacingArray[ii];
        }
      sliceImage->SetSpacing( sliceSpacing );
      }

    if( dictionary.HasKey( "SliceOrigin" ) )
      {
      typedef itk::Array< double > SliceOriginArrayType;
      SliceOriginArrayType sliceOriginArray( SliceDimension );
      itk::ExposeMetaData< SliceOriginArrayType >( dictionary, "SliceOrigin", sliceOriginArray );
      typename SliceImageType::PointType sliceOrigin;
      for( unsigned int ii = 0; ii < SliceDimension; ++ii )
        {
        sliceOrigin[ii] = sliceOriginArray[ii];
        }
      sliceImage->SetOrigin( sliceOrigin );
      }

    image->SetSliceImage( sliceImage );

    if( dictionary.HasKey( "ElevationalSliceAngles" ) )
      {
      typedef itk::Array< double > ElevationalSliceAnglesType;
      const typename ImageType::RegionType & largestRegion = image->GetLargestPossibleRegion();
      const typename ImageType::SizeType & largestSize = largestRegion.GetSize();
      const itk::SizeValueType angles = largestSize[SliceDimension];
      ElevationalSliceAnglesType elevationalSliceAngles( angles );
      itk::ExposeMetaData< ElevationalSliceAnglesType >( dictionary, "ElevationalSliceAngles", elevationalSliceAngles );
      typename SliceImageType::PointType sliceOrigin;
      for( unsigned int ii = 0; ii < angles; ++ii )
        {
        typename TransformType::Pointer transform = TransformType::New();
        transform->SetRotation( 0.0, elevationalSliceAngles[ii], 0.0 );
        image->SetSliceTransform( ii, transform );
        }
      }
  }
};


} // end anonymous namespace
#ifdef ITK_HAS_GCC_PRAGMA_DIAG_PUSHPOP
  ITK_GCC_PRAGMA_DIAG_POP()
#else
  ITK_GCC_PRAGMA_DIAG(warning "-Wattributes")
#endif


namespace itk
{

template < typename TOutputImage >
UltrasoundImageFileReader< TOutputImage >
::UltrasoundImageFileReader()
{
}


template < typename TOutputImage >
void
UltrasoundImageFileReader< TOutputImage >
::GenerateOutputInformation()
{
  Superclass::GenerateOutputInformation();
  UltrasoundImageFileReaderDispatch< TOutputImage >::ExtractMetaData( this->GetOutput() );
}

} // end namespace itk

#endif
