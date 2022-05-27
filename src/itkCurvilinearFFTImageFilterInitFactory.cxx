/*=========================================================================
 *
 *  Copyright NumFOCUS
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         https://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#include "UltrasoundExport.h"

#include "itkObjectFactoryBase.h"
#include "itkFFTImageFilterFactory.h"
#include "itkVnlComplexToComplex1DFFTImageFilter.h"
#include "itkVnlForward1DFFTImageFilter.h"
#include "itkVnlInverse1DFFTImageFilter.h"

#include "itkCurvilinearArraySpecialCoordinatesImage.h"
#include "itkCurvilinearFFTImageFilterInitFactory.h"

namespace itk
{
CurvilinearFFTImageFilterInitFactory::CurvilinearFFTImageFilterInitFactory()
{
  CurvilinearFFTImageFilterInitFactory::RegisterFactories();
}

void
CurvilinearFFTImageFilterInitFactory::RegisterFactories()
{
  // Curvilinear -> Curvilinear
  itk::ObjectFactoryBase::RegisterFactory(FFTImageFilterFactory<VnlComplexToComplex1DFFTImageFilter,
                                                                itk::CurvilinearArraySpecialCoordinatesImage,
                                                                itk::CurvilinearArraySpecialCoordinatesImage>::New());
  itk::ObjectFactoryBase::RegisterFactory(
    FFTImageFilterFactory<VnlForward1DFFTImageFilter, itk::CurvilinearArraySpecialCoordinatesImage, itk::CurvilinearArraySpecialCoordinatesImage>::New());
  itk::ObjectFactoryBase::RegisterFactory(
    FFTImageFilterFactory<VnlInverse1DFFTImageFilter, itk::CurvilinearArraySpecialCoordinatesImage, itk::CurvilinearArraySpecialCoordinatesImage>::New());
  // Curvilinear -> Image
  itk::ObjectFactoryBase::RegisterFactory(FFTImageFilterFactory<VnlForward1DFFTImageFilter, itk::CurvilinearArraySpecialCoordinatesImage, itk::Image>::New());
  // Image -> Curvilinear
  itk::ObjectFactoryBase::RegisterFactory(FFTImageFilterFactory<VnlInverse1DFFTImageFilter, itk::Image, itk::CurvilinearArraySpecialCoordinatesImage>::New());
}

// Undocumented API used to register during static initialization.
// DO NOT CALL DIRECTLY.
void Ultrasound_EXPORT
CurvilinearFFTImageFilterInitFactoryRegister__Private()
{
  CurvilinearFFTImageFilterInitFactory::RegisterFactories();
}

} // end namespace itk
