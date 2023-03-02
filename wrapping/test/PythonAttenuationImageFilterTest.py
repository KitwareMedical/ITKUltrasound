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

parser = argparse.ArgumentParser(description="Computes the estimated attentuation in dB/(MHz*cm).")
parser.add_argument("-i", "--spectra-image", help="Main input image", required=True)
parser.add_argument("-m", "--mask-image", help="Label map for regions of interest", required=True)
parser.add_argument("-a", "--attenuation-image", help="Main output image")
parser.add_argument("-o", "--attenuation-mask-image", help="Label map indicating computed pixels")
parser.add_argument("-w", "--number-of-work-units", nargs='?', const=100, type=int, default=100)
parser.add_argument("-d", "--fixed-estimation-depth-mm", nargs='?', const=0.0, type=float, default=0.0)
parser.add_argument("-n", "--consider-negative-attenuations", nargs='?', const=False, type=bool, default=False)
parser.add_argument("-c", "--computation-mode", nargs='?', const=0, type=int, default=0)
args = parser.parse_args()

itk.auto_progress(2)

spectra_image = itk.imread(args.spectra_image, pixel_type=itk.VariableLengthVector[itk.F])
mask_image = itk.imread(args.mask_image)

# trying to use _mhz instead of _m_hz, or _mm instead of _m_m results in an error:
# AttributeError: 'itkAttenuationImageFilterVIF3IF3' object has no attribute
# 'SetFixedEstimationDepthMm'. Did you mean: 'SetFixedEstimationDepthMM'?
# so we use a quirk of snake_case_to_camel_case to work around it
attenuation_image, attenuation_mask_image = itk.attenuation_image_filter(
    spectra_image,
    input_mask_image=mask_image,
    sampling_frequency_m_hz=60,
    frequency_band_start_m_hz=5.0,
    frequency_band_end_m_hz=20.0,
    number_of_work_units=args.number_of_work_units,
    fixed_estimation_depth_m_m=args.fixed_estimation_depth_mm,
    consider_negative_attenuations=args.consider_negative_attenuations,
    computation_mode=args.computation_mode,
)

itk.imwrite(attenuation_image, args.attenuation_image, compression=False)
itk.imwrite(attenuation_mask_image, args.attenuation_mask_image, compression=True)
