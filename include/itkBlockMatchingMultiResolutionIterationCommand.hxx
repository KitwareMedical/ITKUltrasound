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
#ifndef itkBlockMatchingMultiResolutionIterationCommand_hxx
#define itkBlockMatchingMultiResolutionIterationCommand_hxx


namespace itk
{
namespace BlockMatching
{

template <typename TMultiResolutionMethod>
void
MultiResolutionIterationCommand<TMultiResolutionMethod>::Execute(const Object * itkNotUsed(object), const EventObject & event)
{
  if (!(IterationEvent().CheckEvent(&event)))
  {
    return;
  }

  if (m_MultiResolutionMethod.GetPointer() == nullptr)
  {
    itkExceptionMacro(<< "The associated MultiResolutionMethod must be set.");
  }

  m_FixedImagePyramid = m_MultiResolutionMethod->GetFixedImagePyramid();
  m_MovingImagePyramid = m_MultiResolutionMethod->GetMovingImagePyramid();

  m_BlockRadiusCalculator = m_MultiResolutionMethod->GetBlockRadiusCalculator();

  m_SearchRegionImageSource = m_MultiResolutionMethod->GetSearchRegionImageSource();
}

} // end namespace BlockMatching
} // end namespace itk

#endif
