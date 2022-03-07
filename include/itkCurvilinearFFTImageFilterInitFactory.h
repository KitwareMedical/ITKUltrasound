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
#ifndef itkCurvilinearFFTImageFilterInitFactory_h
#define itkCurvilinearFFTImageFilterInitFactory_h
#include "UltrasoundExport.h"

#include "itkLightObject.h"

namespace itk
{
/**
 *\class CurvilinearFFTImageFilterInitFactory
 * \brief Initialize FFT image filter factory backends.
 *
 * The purpose of CurvilinearFFTImageFilterInitFactory is to perform
 * one-time registration of factory objects that handle
 * creation of FFT image filter classes as applied to
 * CurvilinearArraySpecialCoordinatesImage objects
 * through the ITK object factory singleton mechanism.
 *
 * By default ITK registers Vnl FFT backends to operate on
 * rectangular `itk::Image` inputs and outputs. ITKUltrasound
 * expands default availability to apply FFTs to
 * curvilinear images.
 *
 * \ingroup FourierTransform
 * \ingroup ITKFFT
 * \ingroup ITKSystemObjects
 * \ingroup Ultrasound
 */
class Ultrasound_EXPORT CurvilinearFFTImageFilterInitFactory : public LightObject
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(CurvilinearFFTImageFilterInitFactory);

  /** Standard class type aliases. */
  using Self = CurvilinearFFTImageFilterInitFactory;
  using Superclass = LightObject;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  /** Method for class instantiation. */
  itkFactorylessNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(CurvilinearFFTImageFilterInitFactory, LightObject);

  /** Mimic factory interface for Python initialization  */
  static void
  RegisterOneFactory()
  {
    RegisterFactories();
  }

  /** Register all Curvilinear FFT factories */
  static void
  RegisterFactories();

protected:
  CurvilinearFFTImageFilterInitFactory();
  ~CurvilinearFFTImageFilterInitFactory() override = default;
};
} // end namespace itk

#endif // itkCurvilinearFFTImageFilterInitFactory_h
