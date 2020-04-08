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
#ifndef itkBlockMatchingMultiResolutionFixedBlockRadiusCalculator_h
#define itkBlockMatchingMultiResolutionFixedBlockRadiusCalculator_h

#include "itkBlockMatchingMultiResolutionBlockRadiusCalculator.h"
#include "itkObjectFactory.h"

namespace itk
{
namespace BlockMatching
{
/** \class MultiResolutionFixedBlockRadiusCalculator
 *
 * \brief A fixed radius is used for every level.
 *
 * This class generates the fixed image matching kernel radius in a
 * BlockMatching::MultiResolutionImageRegistrationMethod.  A fixed block radius
 * is used, which means the size of the block in physical coordinates then
 * scales with the pyramid schedule.
 *
 * \sa MultiResolutionBlockRadiusCalculator
 *
 * \ingroup Ultrasound
 * */
template< class TFixedImage >
class ITK_TEMPLATE_EXPORT MultiResolutionFixedBlockRadiusCalculator :
  public  MultiResolutionBlockRadiusCalculator< TFixedImage >
{
public:
  /** Standard class type alias. */
  using Self = MultiResolutionFixedBlockRadiusCalculator;
  using Superclass = MultiResolutionBlockRadiusCalculator< TFixedImage >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( MultiResolutionFixedBlockRadiusCalculator, MultiResolutionBlockRadiusCalculator );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** ImageDimension enumeration. */
  itkStaticConstMacro( ImageDimension, unsigned int,
                      Superclass::ImageDimension );

  /** Type of the fixed image. */
  using FixedImageType = typename Superclass::FixedImageType;
  using RadiusType = typename Superclass::RadiusType;

  /** Set the radius used for all pyramid levels. */
  void SetRadius( const RadiusType& radius )
    { m_Radius = radius; }
  itkGetConstReferenceMacro( Radius, RadiusType );

  const RadiusType& Compute( unsigned long level ) override
    {
    return m_Radius;
    }

protected:
  MultiResolutionFixedBlockRadiusCalculator()
   {
   m_Radius.Fill( 1 );
   }

  RadiusType m_Radius;

private:
  MultiResolutionFixedBlockRadiusCalculator( const Self& );
  void operator=( const Self& );
};

} // end namespace itk
} // end namespace BlockMatching

#endif
