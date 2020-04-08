/*=========================================================================
 *
 *  Copyright NumFOCUS
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
#ifndef itkHDF5UltrasoundImageIOFactory_h
#define itkHDF5UltrasoundImageIOFactory_h
#include "UltrasoundExport.h"


#include "itkObjectFactoryBase.h"
#include "itkImageIOBase.h"

namespace itk
{
/** \class HDF5UltrasoundImageIOFactory
   * \brief Create instances of HDF5UltrasoundImageIO objects using an object
   * factory.
   *
   * \author Matt McCormick
   *
   * \ingroup Ultrasound
   */
class Ultrasound_EXPORT HDF5UltrasoundImageIOFactory: public ObjectFactoryBase
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(HDF5UltrasoundImageIOFactory);

  /** Standard class type alias. */
  using Self = HDF5UltrasoundImageIOFactory;
  using Superclass = ObjectFactoryBase;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Class methods used to interface with the registered factories. */
  virtual const char * GetITKSourceVersion() const override;

  virtual const char * GetDescription() const override;

  /** Method for class instantiation. */
  itkFactorylessNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(HDF5UltrasoundImageIOFactory, ObjectFactoryBase);

  /** Register one factory of this type  */
  static void RegisterOneFactory(void)
  {
    HDF5UltrasoundImageIOFactory::Pointer metaFactory = HDF5UltrasoundImageIOFactory::New();

    ObjectFactoryBase::RegisterFactoryInternal(metaFactory);
  }

protected:
  HDF5UltrasoundImageIOFactory();
  ~HDF5UltrasoundImageIOFactory();
  virtual void PrintSelf(std::ostream & os, Indent indent) const override;
};
} // end namespace itk

#endif
