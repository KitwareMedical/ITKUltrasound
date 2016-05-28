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
#include "itkHDF5UltrasoundImageIOFactory.h"
#include "itkHDF5UltrasoundImageIO.h"
#include "itkVersion.h"

namespace itk
{

HDF5UltrasoundImageIOFactory::HDF5UltrasoundImageIOFactory()
{
  this->RegisterOverride( "itkImageIOBase",
                          "itkHDF5UltrasoundImageIO",
                          "HDF5 Ultrasound Image IO",
                          1,
                          CreateObjectFunction< HDF5UltrasoundImageIO >::New() );
}


HDF5UltrasoundImageIOFactory::~HDF5UltrasoundImageIOFactory()
{
}


const char *
HDF5UltrasoundImageIOFactory::GetITKSourceVersion() const
{
  return ITK_SOURCE_VERSION;
}


const char *
HDF5UltrasoundImageIOFactory::GetDescription() const
{
  return "HDF5 ImageIO Factory, allows the loading of HDF5 ultrasound images into ITK";
}


void
HDF5UltrasoundImageIOFactory
::PrintSelf(std::ostream &, Indent) const
{
}


// Undocumented API used to register during static initialization.
// DO NOT CALL DIRECTLY.

static bool HDF5UltrasoundImageIOFactoryHasBeenRegistered;

void Ultrasound_EXPORT HDF5UltrasoundImageIOFactoryRegister__Private()
{
  if( ! HDF5UltrasoundImageIOFactoryHasBeenRegistered )
    {
    HDF5UltrasoundImageIOFactoryHasBeenRegistered = true;
    HDF5UltrasoundImageIOFactory::RegisterOneFactory();
    }
}

} // end namespace itk
