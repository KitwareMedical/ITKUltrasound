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
#ifndef itkHDF5UltrasoundImageIO_h
#define itkHDF5UltrasoundImageIO_h
#include "UltrasoundExport.h"


// itk namespace first suppresses
// kwstyle error for the H5 namespace below
namespace itk
{
}
namespace H5
{
class H5File;
class DataSpace;
class DataSet;
}

#include "itkStreamingImageIOBase.h"

namespace itk
{
class MetaDataObjectBase;
class MetaDataDictionary;
/**
 * \class HDF5UltrasoundImageIO
 *
 * \brief Class that defines how to read ultrasound images stored in the HDF5 file format.
 *
 * Currently supports reading the format used by Duke University to read mechanically
 * rotated linear array volumes.
 *
 * \author Matt McCormick
 *
 * \ingroup Ultrasound
 *
 */

class Ultrasound_EXPORT HDF5UltrasoundImageIO: public StreamingImageIOBase
{
public:
  /** Standard class typedefs. */
  typedef HDF5UltrasoundImageIO Self;
  typedef StreamingImageIOBase  Superclass;
  typedef SmartPointer< Self >  Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(HDF5UltrasoundImageIO, StreamingImageIOBase);

  virtual bool CanReadFile(const char * fileNameToRead) override;

  virtual void ReadImageInformation() override;

  virtual void Read(void *buffer) override;

  virtual bool CanWriteFile(const char * fileNameToWrite) override;

  virtual void WriteImageInformation() override;

  virtual void Write(const void *buffer) override;

protected:
  HDF5UltrasoundImageIO();
  ~HDF5UltrasoundImageIO();

  virtual SizeType GetHeaderSize() const override;

  virtual void PrintSelf(std::ostream & os, Indent indent) const override;

private:
  HDF5UltrasoundImageIO(const Self &) ITK_DELETED_FUNCTION;
  void operator=(const Self &) ITK_DELETED_FUNCTION;

  template <typename TScalar>
  std::vector<TScalar> ReadVector(const std::string &dataSetName);

  void SetupStreaming(H5::DataSpace *imageSpace,
                      H5::DataSpace *slabSpace);

  void CloseH5File();

  H5::H5File  * m_H5File;
  H5::DataSet * m_VoxelDataSet;
};
} // end namespace itk

#endif // itkHDF5UltrasoundImageIO_h
