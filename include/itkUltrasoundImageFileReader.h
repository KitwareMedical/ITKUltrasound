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
#ifndef itkUltrasoundImageFileReader_h
#define itkUltrasoundImageFileReader_h

#include "itkImageFileReader.h"

#include "itkCurvilinearArraySpecialCoordinatesImage.h"

#ifdef ITK_HAS_GCC_PRAGMA_DIAG_PUSHPOP
  ITK_GCC_PRAGMA_DIAG_PUSH()
#endif
ITK_GCC_PRAGMA_DIAG(ignored "-Wattributes")
namespace itk
{

/**
 * \class UltrasoundImageFileReader
 *
 * Read an ultrasound itk::SpecialCoordinatesImage and populate its parameters
 * based on its itk::MetaDataDictionary entries.
 *
 * \ingroup Ultrasound
 *
 * \sa ImageFileReader
 * \sa SpecialCoordinatesImage
 * \sa CurvilinearArraySpecialCoordinatesImage
 *
 */
template < typename TOutputImage >
class ITK_TEMPLATE_EXPORT UltrasoundImageFileReader :
  public ImageFileReader< TOutputImage >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(UltrasoundImageFileReader);

  /** Standard class type alias.   */
  using Self = UltrasoundImageFileReader;
  using Superclass = ImageFileReader< TOutputImage >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(UltrasoundImageFileReader, ImageFileReader);

  using OutputImageType = TOutputImage;
  using OutputImagePixelType = typename Superclass::OutputImagePixelType;

  itkStaticConstMacro(ImageDimension, unsigned int, OutputImageType::ImageDimension);

protected:
  UltrasoundImageFileReader();
  ~UltrasoundImageFileReader() {}

  /** Prepare the allocation of the output image during the first back
   * propagation of the pipeline. */
  virtual void GenerateOutputInformation() override;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkUltrasoundImageFileReader.hxx"
#endif

#ifdef ITK_HAS_GCC_PRAGMA_DIAG_PUSHPOP
  ITK_GCC_PRAGMA_DIAG_POP()
#else
  ITK_GCC_PRAGMA_DIAG(warning "-Wattributes")
#endif

#endif // itkUltrasoundImageFileReader_h
