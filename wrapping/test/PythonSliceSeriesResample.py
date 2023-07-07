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

import math
import itk
import argparse

parser = argparse.ArgumentParser(description="Estimate three back-scatter coefficients.")
parser.add_argument("-i", "--input-image", required=True)
parser.add_argument("-o", "--output-image", required=True)
args = parser.parse_args()

itk.auto_progress(2)

# define the types we will use
pixel_type = itk.UC
dimension = 3
slice_type = itk.Image[pixel_type, dimension - 1]
transform_type = itk.Euler3DTransform[itk.D]
image_type = itk.SliceSeriesSpecialCoordinatesImage[slice_type, transform_type]

# read the image
reader = itk.UltrasoundImageFileReader[image_type].New()
reader.SetFileName(args.input_image)
reader.Update()
image = reader.GetOutput()
image.DisconnectPipeline()

print("Verify resampling works with SliceSeriesSpecialCoordinatesImage input")
output_size = [46, 55, 82]
output_spacing = [1.0] * dimension
output_origin = [0.0, -27.2693, -40.6222]
result = itk.resample_image_filter(image, size=output_size, output_spacing=output_spacing, output_origin=output_origin)
itk.imwrite(result, args.output_image)
print(f"Image written to {args.output_image}")
