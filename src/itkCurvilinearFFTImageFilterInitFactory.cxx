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
#include "itkMacro.h"

#if ITK_VERSION_MAJOR >= 6
#  include "itkPocketFFTComplexToComplex1DFFTImageFilter.h"
#  include "itkPocketFFTForward1DFFTImageFilter.h"
#  include "itkPocketFFTInverse1DFFTImageFilter.h"
#  define ITKULTRASOUND_COMPLEX_TO_COMPLEX_1DFFT PocketFFTComplexToComplex1DFFTImageFilter
#  define ITKULTRASOUND_FORWARD_1DFFT PocketFFTForward1DFFTImageFilter
#  define ITKULTRASOUND_INVERSE_1DFFT PocketFFTInverse1DFFTImageFilter
#else
#  include "itkVnlComplexToComplex1DFFTImageFilter.h"
#  include "itkVnlForward1DFFTImageFilter.h"
#  include "itkVnlInverse1DFFTImageFilter.h"
#  define ITKULTRASOUND_COMPLEX_TO_COMPLEX_1DFFT VnlComplexToComplex1DFFTImageFilter
#  define ITKULTRASOUND_FORWARD_1DFFT VnlForward1DFFTImageFilter
#  define ITKULTRASOUND_INVERSE_1DFFT VnlInverse1DFFTImageFilter
#endif

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
  itk::ObjectFactoryBase::RegisterFactory(FFTImageFilterFactory<ITKULTRASOUND_COMPLEX_TO_COMPLEX_1DFFT,
                                                                itk::CurvilinearArraySpecialCoordinatesImage,
                                                                itk::CurvilinearArraySpecialCoordinatesImage>::New());
  itk::ObjectFactoryBase::RegisterFactory(FFTImageFilterFactory<ITKULTRASOUND_FORWARD_1DFFT,
                                                                itk::CurvilinearArraySpecialCoordinatesImage,
                                                                itk::CurvilinearArraySpecialCoordinatesImage>::New());
  itk::ObjectFactoryBase::RegisterFactory(FFTImageFilterFactory<ITKULTRASOUND_INVERSE_1DFFT,
                                                                itk::CurvilinearArraySpecialCoordinatesImage,
                                                                itk::CurvilinearArraySpecialCoordinatesImage>::New());
  // Curvilinear -> Image
  itk::ObjectFactoryBase::RegisterFactory(FFTImageFilterFactory<ITKULTRASOUND_FORWARD_1DFFT,
                                                                itk::CurvilinearArraySpecialCoordinatesImage,
                                                                itk::Image>::New());
  // Image -> Curvilinear
  itk::ObjectFactoryBase::RegisterFactory(FFTImageFilterFactory<ITKULTRASOUND_INVERSE_1DFFT,
                                                                itk::Image,
                                                                itk::CurvilinearArraySpecialCoordinatesImage>::New());
}

// Undocumented API used to register during static initialization.
// DO NOT CALL DIRECTLY.
void Ultrasound_EXPORT
CurvilinearFFTImageFilterInitFactoryRegister__Private()
{
  CurvilinearFFTImageFilterInitFactory::RegisterFactories();
}

} // end namespace itk
