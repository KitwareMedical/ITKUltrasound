#==========================================================================
#
#   Copyright NumFOCUS
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#          https://www.apache.org/licenses/LICENSE-2.0.txt
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
#==========================================================================*/


import itk
import sys

itk.auto_progress(2)

PixelType = itk.VariableLengthVector[itk.F]
ImageDimension = 3
ImageType = itk.VectorImage[itk.F, ImageDimension]
ImageType2 = itk.VectorImage[itk.F, 2]

print("Reading input")
input = itk.imread(sys.argv[1], pixel_type=PixelType)
print("Reading reference")
# reference = itk.imread(sys.argv[2], pixel_type=PixelType) # This produces VIF3
reader = itk.ImageFileReader[ImageType2].New()
reader.SetFileName(sys.argv[2])
reader.Update()
reference = reader.GetOutput() # This VIF2, which is what we need

print("Running the filter")
result = itk.spectra1_d_normalize_image_filter(input, reference_image=reference)

print(f"Writing resulting image into file: {sys.argv[2]}")
itk.imwrite(result, sys.argv[2], compression=False)

print("Test finished")
