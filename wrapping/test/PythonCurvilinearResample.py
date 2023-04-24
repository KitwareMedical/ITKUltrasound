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
parser.add_argument("-a", "--lateral-angular-separation", type=float, required=True)
parser.add_argument("-r", "--radius-sample-size", type=float, required=True)
parser.add_argument("-f", "--first-sample-distance", type=float, required=True)
args = parser.parse_args()

itk.auto_progress(2)

pixel_type = itk.UC
dimension = 3
image_type = itk.CurvilinearArraySpecialCoordinatesImage[pixel_type, dimension]
reader = itk.ImageFileReader[image_type].New()
reader.SetFileName(args.input_image)
reader.Update()
image = reader.GetOutput()
image.DisconnectPipeline()

image.SetLateralAngularSeparation(args.lateral_angular_separation)
image.SetFirstSampleDistance(args.first_sample_distance)
image.SetRadiusSampleSize(args.radius_sample_size)
print(image)

print("Verify resampling works with curvilinear image input")
output_size = [800, 800, 1]
output_spacing = [0.15] * dimension
output_origin = [output_size[0] * output_spacing[0] / -2.0,
                 args.first_sample_distance * math.cos(math.pi / 4.0), 0.0]
result = itk.resample_image_filter(image, size=output_size, output_spacing=output_spacing, output_origin=output_origin)
itk.imwrite(result, args.output_image)
