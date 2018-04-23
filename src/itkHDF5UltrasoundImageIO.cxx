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
#include "itkVersion.h"
#include "itkHDF5UltrasoundImageIO.h"
#include "itkMetaDataObject.h"
#include "itkArray.h"
#include "itksys/SystemTools.hxx"
#include "itk_H5Cpp.h"

namespace itk
{

HDF5UltrasoundImageIO
::HDF5UltrasoundImageIO() : m_H5File(ITK_NULLPTR),
                            m_VoxelDataSet(ITK_NULLPTR)
{
}


HDF5UltrasoundImageIO
::~HDF5UltrasoundImageIO()
{
  if(this->m_VoxelDataSet != ITK_NULLPTR)
    {
    m_VoxelDataSet->close();
    delete m_VoxelDataSet;
    }
  this->CloseH5File();
}

namespace
{
template <typename TScalar>
H5::PredType GetType()
{
  itkGenericExceptionMacro(<< "Type not handled "
                           << "in HDF5 File: "
                           << typeid(TScalar).name());

}
#define GetH5TypeSpecialize(CXXType,H5Type) \
  template <>                             \
  H5::PredType GetType<CXXType>()         \
  {                                       \
    return H5Type;                        \
  }

GetH5TypeSpecialize(float,                  H5::PredType::NATIVE_FLOAT)
GetH5TypeSpecialize(double,                 H5::PredType::NATIVE_DOUBLE)

GetH5TypeSpecialize(char,                   H5::PredType::NATIVE_CHAR)
GetH5TypeSpecialize(unsigned char,          H5::PredType::NATIVE_UCHAR)

GetH5TypeSpecialize(short int,              H5::PredType::NATIVE_SHORT)
GetH5TypeSpecialize(short unsigned int,     H5::PredType::NATIVE_USHORT)

GetH5TypeSpecialize(int,                    H5::PredType::NATIVE_INT)
GetH5TypeSpecialize(unsigned int,           H5::PredType::NATIVE_UINT)

GetH5TypeSpecialize(long int,               H5::PredType::NATIVE_LONG)
GetH5TypeSpecialize(long unsigned int,      H5::PredType::NATIVE_ULONG)

#if defined(_MSC_VER) && defined(ITK_USE_64BITS_IDS) && ((ULLONG_MAX != ULONG_MAX) || (LLONG_MAX != LONG_MAX))
GetH5TypeSpecialize(long long int,          H5::PredType::NATIVE_LLONG)
GetH5TypeSpecialize(unsigned long long int, H5::PredType::NATIVE_ULLONG)
#endif

/* The following types are not implmented.  This comment serves
 * to indicate that the full complement of possible H5::PredType
 * types are not implemented int the ITK IO reader/writer
 * GetH5TypeSpecialize(bool,              H5::PredType::NATIVE_HBOOL)
*/

#undef GetH5TypeSpecialize

inline
ImageIOBase::IOComponentType
PredTypeToComponentType(const H5::DataType & type)
{
  if(type ==  H5::PredType::NATIVE_UCHAR)
    {
    return ImageIOBase::UCHAR;
    }
  else if(type ==  H5::PredType::NATIVE_CHAR)
    {
    return ImageIOBase::CHAR;
    }
  else if(type ==  H5::PredType::NATIVE_USHORT)
    {
    return ImageIOBase::USHORT;
    }
  else if(type ==  H5::PredType::NATIVE_SHORT)
    {
    return ImageIOBase::SHORT;
    }
  else if(type ==  H5::PredType::NATIVE_UINT)
    {
    return ImageIOBase::UINT;
    }
  else if(type ==  H5::PredType::NATIVE_INT)
    {
    return ImageIOBase::INT;
    }
  else if(type ==  H5::PredType::NATIVE_ULLONG)
    {
    return ImageIOBase::ULONG;
    }
  else if(type ==  H5::PredType::NATIVE_LLONG)
    {
    return ImageIOBase::LONG;
    }
  else if(type ==  H5::PredType::NATIVE_FLOAT)
    {
    return ImageIOBase::FLOAT;
    }
  else if(type ==  H5::PredType::NATIVE_DOUBLE)
    {
    return ImageIOBase::DOUBLE;
    }
  else if(type ==  H5::PredType::NATIVE_ULONG)
    {
    if(sizeof(unsigned int) == sizeof(unsigned long))
      {
      return ImageIOBase::UINT;
      }
    else if(sizeof(unsigned int) == sizeof(unsigned long long))
      {
      return ImageIOBase::ULONGLONG;
      }
    }
  else if(type ==  H5::PredType::NATIVE_LONG)
    {
    if(sizeof(int) == sizeof(long))
      {
      return ImageIOBase::INT;
      }
    else if(sizeof(int) == sizeof(long long))
      {
      return ImageIOBase::LONGLONG;
      }
    }
  itkGenericExceptionMacro(<< "unsupported data type "
                           << type.fromClass());
}

H5::PredType
ComponentToPredType(ImageIOBase::IOComponentType cType)
{
  switch ( cType )
    {
    case ImageIOBase::UCHAR:
      return H5::PredType::NATIVE_UCHAR;
    case ImageIOBase::CHAR:
      return H5::PredType::NATIVE_CHAR;
    case ImageIOBase::USHORT:
      return H5::PredType::NATIVE_USHORT;
    case ImageIOBase::SHORT:
      return H5::PredType::NATIVE_SHORT;
    case ImageIOBase::UINT:
      return H5::PredType::NATIVE_UINT;
    case ImageIOBase::INT:
      return H5::PredType::NATIVE_INT;
    case ImageIOBase::ULONG:
      return H5::PredType::NATIVE_ULONG;
    case ImageIOBase::ULONGLONG:
      return H5::PredType::NATIVE_ULLONG;
    case ImageIOBase::LONG:
      return H5::PredType::NATIVE_LONG;
    case ImageIOBase::LONGLONG:
      return H5::PredType::NATIVE_LLONG;
    case ImageIOBase::FLOAT:
      return H5::PredType::NATIVE_FLOAT;
    case ImageIOBase::DOUBLE:
      return H5::PredType::NATIVE_DOUBLE;
    case ImageIOBase::UNKNOWNCOMPONENTTYPE:
      itkGenericExceptionMacro(<< "unsupported IOComponentType"
                               << cType);
    }

    itkGenericExceptionMacro(<< "unsupported IOComponentType"
                             << cType);
}

std::string
ComponentToString(ImageIOBase::IOComponentType cType)
{
  std::string rval;
  switch ( cType )
    {
    case ImageIOBase::UCHAR:
      rval = "UCHAR";
      break;
    case ImageIOBase::CHAR:
      rval = "CHAR";
      break;
    case ImageIOBase::USHORT:
      rval = "USHORT";
      break;
    case ImageIOBase::SHORT:
      rval = "SHORT";
      break;
    case ImageIOBase::UINT:
      rval = "UINT";
      break;
    case ImageIOBase::INT:
      rval = "INT";
      break;
    case ImageIOBase::ULONG:
      rval = "ULONG";
      break;
    case ImageIOBase::LONG:
      rval = "LONG";
      break;
    case ImageIOBase::FLOAT:
      rval = "FLOAT";
      break;
    case ImageIOBase::DOUBLE:
      rval = "DOUBLE";
      break;
    default:
      itkGenericExceptionMacro(<< "unsupported IOComponentType"
                               << cType);
    }
  return rval;
}

// Function:    H5Object::doesAttrExist
///\brief       test for existence of attribut
///\param       name - IN: Name of the attribute
///\return      true if attribute exists, false otherwise
///\exception   none
// Programmer   Kent Williams 2011
//--------------------------------------------------------------------------
static bool doesAttrExist(const H5::H5Object &object, const char * const name )
{
  return( H5Aexists(object.getId(), name) > 0 ? true : false );
}

} // end anonymous namespace


template <typename TScalar>
std::vector<TScalar>
HDF5UltrasoundImageIO
::ReadVector(const std::string & dataSetName)
{
  std::vector< TScalar > result;
  H5::DataSet dataSet = this->m_H5File->openDataSet(dataSetName);
  H5::DataSpace space = dataSet.getSpace();

  std::vector<hsize_t> dim(space.getSimpleExtentNdims());
  space.getSimpleExtentDims( &dim[0], ITK_NULLPTR );
  result.resize( dim[0] );

  H5::PredType vecType = GetType<TScalar>();
  dataSet.read( &(result[0]), vecType );
  dataSet.close();
  return result;
}


// This method will only test if the header looks like an
// HDF5 Header.  Some code is redundant with ReadImageInformation
// a StateMachine could provide a better implementation
bool
HDF5UltrasoundImageIO
::CanReadFile(const char *fileNameToRead)
{
  //HDF5 is overly verbose in complaining that
  //     a file does not exist.
  if ( !itksys::SystemTools::FileExists(fileNameToRead) )
    {
    return false;
    }

  //Do not read if it is a MINC file.
  std::string filename(fileNameToRead);
  std::string::size_type mncPos = filename.rfind(".mnc");
  if ( (mncPos != std::string::npos)
       && (mncPos == filename.length() - 4) )
    {
    return false;
    }

  mncPos = filename.rfind(".MNC");
  if ( (mncPos != std::string::npos)
       && (mncPos == filename.length() - 4) )
    {
    return false;
    }

  mncPos = filename.rfind(".mnc2");
  if ( (mncPos != std::string::npos)
       && (mncPos == filename.length() - 5) )
    {
    return false;
    }

  mncPos = filename.rfind(".MNC2");
  if ( (mncPos != std::string::npos)
       && (mncPos == filename.length() - 5) )
    {
    return false;
    }


  // call standard method to determine HDF-ness
  bool rval;
  // HDF5 is so exception happy, we have to worry about
  // it throwing a wobbly here if the file doesn't exist
  // or has some other problem.
  try
    {
    rval = H5::H5File::isHdf5(fileNameToRead);
    }
  catch(...)
    {
    return false;
    }

  this->CloseH5File();
  this->m_H5File = new H5::H5File( fileNameToRead, H5F_ACC_RDONLY );
  // Do not try to read this file if it is an ITK HDF5 Image file.
  const herr_t status = H5Lexists( this->m_H5File->getId(), "/ITKImage", H5P_DEFAULT );
  if( status == 1 )
    {
    rval = false;
    }

  return rval;
}


void
HDF5UltrasoundImageIO
::CloseH5File()
{
  if(this->m_H5File != ITK_NULLPTR)
    {
    this->m_H5File->close();
    delete this->m_H5File;
    }
}


void
HDF5UltrasoundImageIO
::ReadImageInformation()
{
  try
    {
    this->CloseH5File();
    this->m_H5File = new H5::H5File(this->GetFileName(),
                                    H5F_ACC_RDONLY);

    this->SetNumberOfDimensions( 3 );

    const std::string axialPixelLocationsDataSet( "/axial" );
    typedef std::vector< double > AxialPixelLocationsType;
    AxialPixelLocationsType axialPixelLocations = this->ReadVector< double >( axialPixelLocationsDataSet );
    this->SetDimensions( 0, axialPixelLocations.size() );

    const std::string lateralPixelLocationsDataSet( "/lat" );
    typedef std::vector< double > LateralPixelLocationsType;
    LateralPixelLocationsType lateralPixelLocations = this->ReadVector< double >( lateralPixelLocationsDataSet );
    this->SetDimensions( 1, lateralPixelLocations.size() );

    const std::string elevationalSliceAngleDataSet( "/eleAngle" );
    typedef std::vector< double > ElevationalSliceAngleType;
    ElevationalSliceAngleType elevationalSliceAngles = this->ReadVector< double >( elevationalSliceAngleDataSet );
    const size_t angles = elevationalSliceAngles.size();
    for( size_t ii = 0; ii < angles; ++ii )
      {
      elevationalSliceAngles[ii] *= Math::pi_over_180;
      }
    this->SetDimensions( 2, angles );

    // set the ComponentType
    const std::string pixelDataName( "/bimg" );
    const H5::DataSet pixelDataSet = this->m_H5File->openDataSet( pixelDataName );
    const H5::DataType pixelDataType = pixelDataSet.getDataType();
    this->SetComponentType( PredTypeToComponentType( pixelDataType ) );


    // Read out metadata
    MetaDataDictionary & metaDataDict = this->GetMetaDataDictionary();
    // It is necessary to clear dict if ImageIO object is re-used
    metaDataDict.Clear();

    EncapsulateMetaData< std::string >( metaDataDict, "SliceType", "Image" );

    typedef Array< double > SliceSpacingType;
    SliceSpacingType sliceSpacing( 2 );
    sliceSpacing[0] = ( axialPixelLocations[axialPixelLocations.size() - 1] - axialPixelLocations[0] ) / ( axialPixelLocations.size() - 1 );
    sliceSpacing[1] = ( lateralPixelLocations[lateralPixelLocations.size() - 1] - lateralPixelLocations[0] ) / ( lateralPixelLocations.size() - 1 );
    EncapsulateMetaData< SliceSpacingType >( metaDataDict, "SliceSpacing", sliceSpacing );

    typedef Array< double > SliceOriginType;
    SliceOriginType sliceOrigin( 2 );
    sliceOrigin[0] = axialPixelLocations[0];
    sliceOrigin[1] = lateralPixelLocations[1];
    EncapsulateMetaData< SliceOriginType >( metaDataDict, "SliceOrigin", sliceOrigin );

    typedef Array< double > ElevationalSliceAnglesMetaDataType;
    ElevationalSliceAnglesMetaDataType elevationalSliceAnglesMetaData( &(elevationalSliceAngles[0]), elevationalSliceAngles.size() );
    EncapsulateMetaData< ElevationalSliceAnglesMetaDataType >( metaDataDict, "ElevationalSliceAngles", elevationalSliceAnglesMetaData );
    }
  // catch failure caused by the H5File operations
  catch( H5::AttributeIException & error )
    {
    itkExceptionMacro(<< error.getCDetailMsg());
    }
  catch( H5::FileIException & error )
    {
    itkExceptionMacro(<< error.getCDetailMsg());
    }
  // catch failure caused by the DataSet operations
  catch( H5::DataSetIException & error )
    {
    itkExceptionMacro(<< error.getCDetailMsg());
    }
  // catch failure caused by the DataSpace operations
  catch( H5::DataSpaceIException & error )
    {
    itkExceptionMacro(<< error.getCDetailMsg());
    }
  // catch failure caused by the DataSpace operations
  catch( H5::DataTypeIException & error )
    {
    itkExceptionMacro(<< error.getCDetailMsg());
    }
}

void
HDF5UltrasoundImageIO
::SetupStreaming(H5::DataSpace *imageSpace, H5::DataSpace *slabSpace)
{
  const ImageIORegion            regionToRead = this->GetIORegion();
  const ImageIORegion::SizeType  size = regionToRead.GetSize();
  const ImageIORegion::IndexType start = regionToRead.GetIndex();
  //
  const int numComponents = this->GetNumberOfComponents();

  const int HDFDim(this->GetNumberOfDimensions() +
                   (numComponents > 1 ? 1 : 0));

  hsize_t *offset = new hsize_t[HDFDim];
  hsize_t *HDFSize = new hsize_t[HDFDim];
  const int limit = regionToRead.GetImageDimension();
  //
  // fastest moving dimension is intra-voxel
  // index
  int i = 0;
  if(numComponents > 1)
    {
    offset[HDFDim - 1] = 0;
    HDFSize[HDFDim - 1] = numComponents;
    ++i;
    }

  // HDF5 dimensions listed slowest moving first, ITK are fastest
  // moving first.
  for(int j = 0; j < limit && i < HDFDim; ++i, ++j )
    {
      offset[HDFDim - i - 1] = start[j];
      HDFSize[HDFDim - i - 1] = size[j];
    }

  while ( i < HDFDim )
    {
    offset[HDFDim - i - 1] = 0;
    HDFSize[HDFDim - i - 1] = 1;
    ++i;
    }

  slabSpace->setExtentSimple(HDFDim,HDFSize);
  imageSpace->selectHyperslab(H5S_SELECT_SET,HDFSize,offset);
  delete[] HDFSize;
  delete[] offset;
}

void
HDF5UltrasoundImageIO
::Read(void *buffer)
{
  ImageIORegion            regionToRead = this->GetIORegion();
  ImageIORegion::SizeType  size = regionToRead.GetSize();
  ImageIORegion::IndexType start = regionToRead.GetIndex();

  std::string pixelDataName( "/bimg" );
  H5::DataSet pixelDataSet = this->m_H5File->openDataSet( pixelDataName );
  H5::DataType voxelType = pixelDataSet.getDataType();
  H5::DataSpace imageSpace = pixelDataSet.getSpace();

  H5::DataSpace dataSpace;
  this->SetupStreaming( &imageSpace, &dataSpace );
  pixelDataSet.read( buffer, voxelType, dataSpace, imageSpace);
}


bool
HDF5UltrasoundImageIO
::CanWriteFile(const char *FileNameToWrite)
{
  return false;
}


void
HDF5UltrasoundImageIO
::WriteImageInformation()
{
  return;
}

void
HDF5UltrasoundImageIO
::Write(const void *buffer)
{
  return;
}


ImageIOBase::SizeType
HDF5UltrasoundImageIO
::GetHeaderSize() const
{
  return 0;
}


void
HDF5UltrasoundImageIO
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  // just prints out the pointer value.
  os << indent << "H5File: " << this->m_H5File << std::endl;
}

} // end namespace itk
