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
/*=========================================================================
 *
 *  Portions of this file are subject to the VTK Toolkit Version 3 copyright.
 *
 *  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 *
 *  For complete copyright, license and disclaimer of warranty information
 *  please refer to the NOTICE file at the top of the ITK source tree.
 *
 *=========================================================================*/
#ifndef itkCurvilinearArraySpecialCoordinatesImage_hxx
#define itkCurvilinearArraySpecialCoordinatesImage_hxx
#include "itkCurvilinearArraySpecialCoordinatesImage.h"

#include <complex>

namespace itk
{

template< typename TPixel, unsigned int VDimension >
void
CurvilinearArraySpecialCoordinatesImage< TPixel, VDimension >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent
     << "RadiusSampleSize = " << m_RadiusSampleSize
     << std::endl;
  os << indent
     << "LateralAngularSeparation = "
     << m_LateralAngularSeparation
     << std::endl;
  os << indent
     << "FirstSampleDistance = "
     << m_FirstSampleDistance
     << std::endl;
}


template< typename TPixel, unsigned int VDimension >
void
CurvilinearArraySpecialCoordinatesImage< TPixel, VDimension >
::CopyInformation(const DataObject *data)
{
  Superclass::CopyInformation( data );

  if ( data )
    {
    // Attempt to cast data to an CurvilinearArraySpecialCoordinatesImage
    const CurvilinearArraySpecialCoordinatesImage< TPixel, VDimension > * const curvilinearArray =
      dynamic_cast< const CurvilinearArraySpecialCoordinatesImage< TPixel, VDimension > * >( data );
    const CurvilinearArraySpecialCoordinatesImage< unsigned char, VDimension > * const curvilinearArrayUnsignedChar =
      dynamic_cast< const CurvilinearArraySpecialCoordinatesImage< unsigned char, VDimension > * >( data );
    const CurvilinearArraySpecialCoordinatesImage< signed char, VDimension > * const curvilinearArraySignedChar =
      dynamic_cast< const CurvilinearArraySpecialCoordinatesImage< signed char, VDimension > * >( data );
    const CurvilinearArraySpecialCoordinatesImage< unsigned short, VDimension > * const curvilinearArrayUnsignedShort =
      dynamic_cast< const CurvilinearArraySpecialCoordinatesImage< unsigned short, VDimension > * >( data );
    const CurvilinearArraySpecialCoordinatesImage< short, VDimension > * const curvilinearArrayShort =
      dynamic_cast< const CurvilinearArraySpecialCoordinatesImage< short, VDimension > * >( data );
    const CurvilinearArraySpecialCoordinatesImage< unsigned int, VDimension > * const curvilinearArrayUnsignedInt =
      dynamic_cast< const CurvilinearArraySpecialCoordinatesImage< unsigned int, VDimension > * >( data );
    const CurvilinearArraySpecialCoordinatesImage< int, VDimension > * const curvilinearArrayInt =
      dynamic_cast< const CurvilinearArraySpecialCoordinatesImage< int, VDimension > * >( data );
    const CurvilinearArraySpecialCoordinatesImage< unsigned long long, VDimension > * const curvilinearArrayUnsignedLongLong =
      dynamic_cast< const CurvilinearArraySpecialCoordinatesImage< unsigned long long, VDimension > * >( data );
    const CurvilinearArraySpecialCoordinatesImage< long long, VDimension > * const curvilinearArrayLongLong =
      dynamic_cast< const CurvilinearArraySpecialCoordinatesImage< long long, VDimension > * >( data );
    const CurvilinearArraySpecialCoordinatesImage< float, VDimension > * const curvilinearArrayFloat =
      dynamic_cast< const CurvilinearArraySpecialCoordinatesImage< float, VDimension > * >( data );
    const CurvilinearArraySpecialCoordinatesImage< double, VDimension > * const curvilinearArrayDouble =
      dynamic_cast< const CurvilinearArraySpecialCoordinatesImage< double, VDimension > * >( data );
    const CurvilinearArraySpecialCoordinatesImage< std::complex< float >, VDimension > * const curvilinearArrayComplexFloat =
      dynamic_cast< const CurvilinearArraySpecialCoordinatesImage< std::complex< float >, VDimension > * >( data );
    const CurvilinearArraySpecialCoordinatesImage< std::complex< double >, VDimension > * const curvilinearArrayComplexDouble =
      dynamic_cast< const CurvilinearArraySpecialCoordinatesImage< std::complex< double >, VDimension > * >( data );
    if ( curvilinearArray != ITK_NULLPTR )
      {
      this->SetLateralAngularSeparation( curvilinearArray->GetLateralAngularSeparation() );
      this->SetRadiusSampleSize( curvilinearArray->GetRadiusSampleSize() );
      this->SetFirstSampleDistance( curvilinearArray->GetFirstSampleDistance() );
      }
    else if ( curvilinearArray != ITK_NULLPTR )
      {
      this->SetLateralAngularSeparation( curvilinearArray->GetLateralAngularSeparation() );
      this->SetRadiusSampleSize( curvilinearArray->GetRadiusSampleSize() );
      this->SetFirstSampleDistance( curvilinearArray->GetFirstSampleDistance() );
      }
    else if ( curvilinearArrayUnsignedChar != ITK_NULLPTR )
      {
      this->SetLateralAngularSeparation( curvilinearArrayUnsignedChar->GetLateralAngularSeparation() );
      this->SetRadiusSampleSize( curvilinearArrayUnsignedChar->GetRadiusSampleSize() );
      this->SetFirstSampleDistance( curvilinearArrayUnsignedChar->GetFirstSampleDistance() );
      }
    else if ( curvilinearArraySignedChar != ITK_NULLPTR )
      {
      this->SetLateralAngularSeparation( curvilinearArraySignedChar->GetLateralAngularSeparation() );
      this->SetRadiusSampleSize( curvilinearArraySignedChar->GetRadiusSampleSize() );
      this->SetFirstSampleDistance( curvilinearArraySignedChar->GetFirstSampleDistance() );
      }
    else if ( curvilinearArrayUnsignedShort != ITK_NULLPTR )
      {
      this->SetLateralAngularSeparation( curvilinearArrayUnsignedShort->GetLateralAngularSeparation() );
      this->SetRadiusSampleSize( curvilinearArrayUnsignedShort->GetRadiusSampleSize() );
      this->SetFirstSampleDistance( curvilinearArrayUnsignedShort->GetFirstSampleDistance() );
      }
    else if ( curvilinearArrayShort != ITK_NULLPTR )
      {
      this->SetLateralAngularSeparation( curvilinearArrayShort->GetLateralAngularSeparation() );
      this->SetRadiusSampleSize( curvilinearArrayShort->GetRadiusSampleSize() );
      this->SetFirstSampleDistance( curvilinearArrayShort->GetFirstSampleDistance() );
      }
    else if ( curvilinearArrayUnsignedInt != ITK_NULLPTR )
      {
      this->SetLateralAngularSeparation( curvilinearArrayUnsignedInt->GetLateralAngularSeparation() );
      this->SetRadiusSampleSize( curvilinearArrayUnsignedInt->GetRadiusSampleSize() );
      this->SetFirstSampleDistance( curvilinearArrayUnsignedInt->GetFirstSampleDistance() );
      }
    else if ( curvilinearArrayInt != ITK_NULLPTR )
      {
      this->SetLateralAngularSeparation( curvilinearArrayInt->GetLateralAngularSeparation() );
      this->SetRadiusSampleSize( curvilinearArrayInt->GetRadiusSampleSize() );
      this->SetFirstSampleDistance( curvilinearArrayInt->GetFirstSampleDistance() );
      }
    else if ( curvilinearArrayUnsignedLongLong != ITK_NULLPTR )
      {
      this->SetLateralAngularSeparation( curvilinearArrayUnsignedLongLong->GetLateralAngularSeparation() );
      this->SetRadiusSampleSize( curvilinearArrayUnsignedLongLong->GetRadiusSampleSize() );
      this->SetFirstSampleDistance( curvilinearArrayUnsignedLongLong->GetFirstSampleDistance() );
      }
    else if ( curvilinearArrayLongLong != ITK_NULLPTR )
      {
      this->SetLateralAngularSeparation( curvilinearArrayLongLong->GetLateralAngularSeparation() );
      this->SetRadiusSampleSize( curvilinearArrayLongLong->GetRadiusSampleSize() );
      this->SetFirstSampleDistance( curvilinearArrayLongLong->GetFirstSampleDistance() );
      }
    else if ( curvilinearArrayFloat != ITK_NULLPTR )
      {
      this->SetLateralAngularSeparation( curvilinearArrayFloat->GetLateralAngularSeparation() );
      this->SetRadiusSampleSize( curvilinearArrayFloat->GetRadiusSampleSize() );
      this->SetFirstSampleDistance( curvilinearArrayFloat->GetFirstSampleDistance() );
      }
    else if ( curvilinearArrayDouble != ITK_NULLPTR )
      {
      this->SetLateralAngularSeparation( curvilinearArrayDouble->GetLateralAngularSeparation() );
      this->SetRadiusSampleSize( curvilinearArrayDouble->GetRadiusSampleSize() );
      this->SetFirstSampleDistance( curvilinearArrayDouble->GetFirstSampleDistance() );
      }
    else if ( curvilinearArrayComplexFloat != ITK_NULLPTR )
      {
      this->SetLateralAngularSeparation( curvilinearArrayComplexFloat->GetLateralAngularSeparation() );
      this->SetRadiusSampleSize( curvilinearArrayComplexFloat->GetRadiusSampleSize() );
      this->SetFirstSampleDistance( curvilinearArrayComplexFloat->GetFirstSampleDistance() );
      }
    else if ( curvilinearArrayComplexDouble != ITK_NULLPTR )
      {
      this->SetLateralAngularSeparation( curvilinearArrayComplexDouble->GetLateralAngularSeparation() );
      this->SetRadiusSampleSize( curvilinearArrayComplexDouble->GetRadiusSampleSize() );
      this->SetFirstSampleDistance( curvilinearArrayComplexDouble->GetFirstSampleDistance() );
      }
    else if ( std::string(data->GetNameOfClass()) == "Image" )
      {
      // No extra meta data to transfer
      }
    else
      {
      // pointer could not be cast back down
      itkExceptionMacro( << "itk::CurvilinearArraySpecialCoordinatesImage::CopyInformation() cannot cast "
                         << typeid( data ).name() << " to "
                         << typeid( const CurvilinearArraySpecialCoordinatesImage * ).name() );
      }
    }
}


template< typename TPixel, unsigned int VDimension >
void
CurvilinearArraySpecialCoordinatesImage< TPixel, VDimension >
::Graft(const DataObject *data)
{
  // call the superclass' implementation
  Superclass::Graft(data);

  if ( data )
    {
    // Attempt to cast data to a CurvilinearArraySpecialCoordinatesImage
    const Self * const imgData = dynamic_cast< const Self * >( data );

    if ( imgData )
      {
      // Now copy anything remaining that is needed
      this->SetPixelContainer( const_cast< PixelContainer * >
                               ( imgData->GetPixelContainer() ) );
      }
    else
      {
      // pointer could not be cast back down
      itkExceptionMacro( << "itk::Image::Graft() cannot cast "
                         << typeid( data ).name() << " to "
                         << typeid( const Self * ).name() );
      }
    }
}

} // end namespace itk

#endif
