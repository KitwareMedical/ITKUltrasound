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
import argparse

parser = argparse.ArgumentParser(description="Estimate three back-scatter coefficients.")
parser.add_argument("-i", "--input-image", required=True)
parser.add_argument("-a", "--average-image", required=True)
parser.add_argument("-s", "--slope-image", required=True)
parser.add_argument("-n", "--intercept-image", required=True)
args = parser.parse_args()

itk.auto_progress(2)

input_image = itk.imread(args.input_image, pixel_type=itk.VariableLengthVector[itk.F])

# trying to use _mhz instead of _m_hz results in an error:
# unknown method name _Mhz, did you mean _MHz? (letter case difference, h vs H)
# so we use a quirk of snake_case_to_camel_case to work around it
average_image, slope_image, intercept_image = itk.backscatter_image_filter(
    input_image,
    number_of_work_units=100,
    sampling_frequency_m_hz=60,
    frequency_band_start_m_hz=5.0,
    frequency_band_end_m_hz=20.0,
)

itk.imwrite(average_image, args.average_image, compression=False)
itk.imwrite(slope_image, args.slope_image, compression=False)
itk.imwrite(intercept_image, args.intercept_image, compression=False)
