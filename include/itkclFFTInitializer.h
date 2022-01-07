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
#if !defined(itkclFFTInitializer_h) && defined(ITKUltrasound_USE_clFFT)
#  define itkclFFTInitializer_h
#  include "UltrasoundExport.h"

#  define __CL_ENABLE_EXCEPTIONS
#  include "CL/cl.hpp"
#  include "clFFT.h"

namespace itk
{
class Ultrasound_EXPORT clFFTInitializer
{
public:
  clfftSetupData m_clFFTdefaults;

  clFFTInitializer();
  ~clFFTInitializer() override
  {
    clfftTeardown();
  }
};

// make sure clFFT has been initialized
Ultrasound_EXPORT clFFTInitializer &
                  clFFFInitialization();

} // namespace itk
#endif // itkclFFTInitializer_h
