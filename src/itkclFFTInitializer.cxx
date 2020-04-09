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
#include "itkclFFTInitializer.h"

namespace itk
{
clFFTInitializer::clFFTInitializer()
{
  clfftInitSetupData(&this->m_clFFTdefaults);
  clfftSetup(&this->m_clFFTdefaults);
}


clFFTInitializer &
clFFFInitialization()
{
  static clFFTInitializer clFFTInitializerGlobal;
  return clFFTInitializerGlobal;
}


bool
clFFFFactorization(size_t dimension)
{
  const unsigned char primeFactors[] = { 2, 3, 5, 7 };
  for (const auto & factor : primeFactors)
  {
    while (dimension % factor == 0)
    {
      dimension /= factor;
    }
  }
  return (dimension == 1);
}
} // namespace itk
