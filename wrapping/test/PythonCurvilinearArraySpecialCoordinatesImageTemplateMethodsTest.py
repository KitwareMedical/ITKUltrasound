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

# Verify itk.CurvilinearArraySpecialCoordinatesImage template methods
# are properly wrapped and behave correctly

import math
import itk
itk.auto_progress(2)

pixel_type = itk.F
dimension = 2
image_type = itk.CurvilinearArraySpecialCoordinatesImage[pixel_type, dimension]

# Set up image
image = image_type.New()
image_size = [2048, 241]
image.SetRegions(image_size)
image.Allocate()

# Metadata is taken from itkCurvilinearArraySpecialCoordinatesImageTest.cxx
lateral_angular_separation = (math.pi / 2.0 + 0.5) / (itk.size(image)[1] - 1)
radius_start = 26.4
radius_stop = 131.5
radius_sample_size = (radius_stop - radius_start) / (itk.size(image)[0] - 1)

image.SetLateralAngularSeparation(lateral_angular_separation)
image.SetFirstSampleDistance(radius_start)
image.SetRadiusSampleSize(radius_sample_size)
print(image)

# Try transforming first index [0,0]
continuous_index = itk.ContinuousIndex[itk.D, 2]()
print(f'Continuous index: {continuous_index}')

expected_point = itk.Point[itk.D,2]([-22.7057, 13.4688])
point = image.TransformContinuousIndexToPhysicalPoint(continuous_index)
print(f'Expected point: Transformed point: {point}')
assert all([abs(el1 - el2) < 1e-3 for el1, el2 in zip(list(point), list(expected_point))]), f'Output point {point} differs from expectation {expected_point}'

transformed_continuous_index = image.TransformPhysicalPointToContinuousIndex(point)
print(f'Transformed continuous index: {continuous_index}')
assert transformed_continuous_index == continuous_index, f'Output {transformed_continuous_index} differs from expectation {continuous_index}'

# Try transforming last index [2047, 240]
continuous_index[0] = 2047
continuous_index[1] = 240
print(f'Continuous index: {continuous_index}')

expected_point = itk.Point[itk.D,2]([113.099, 67.0891])
point = image.TransformContinuousIndexToPhysicalPoint(continuous_index)
print(f'Transformed point: {point}')
assert all([abs(el1 - el2) < 1e-3 for el1, el2 in zip(list(point), list(expected_point))]), f'Output point {point} differs from expectation {expected_point}'

transformed_continuous_index = image.TransformPhysicalPointToContinuousIndex(point)
print(f'Transformed continuous index: {continuous_index}')
assert transformed_continuous_index == continuous_index, f'Output {transformed_continuous_index} differs from expectation {continuous_index}'

# Verify resampling works with curvilinear image input
output_size = [800] * dimension
output_spacing = [0.15] * dimension
output_origin = [output_size[0] * output_spacing[0] / -2.0,
                 radius_start * math.cos(math.pi / 4.0)]
_ = itk.resample_image_filter(image, size=output_size, output_spacing=output_spacing, output_origin=output_origin)
