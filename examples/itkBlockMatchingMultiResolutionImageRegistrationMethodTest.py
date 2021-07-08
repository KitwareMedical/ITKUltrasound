#!/usr/bin/env python

# Copyright NumFOCUS
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0.txt
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
from urllib.request import urlretrieve

import itk

# Verify wrapping succeeded
assert 'MultiResolutionFixedBlockRadiusCalculator' in dir(itk.Ultrasound)

# Get image input for registration
FIXED_IMAGE_PATH = 'Input/rf_pre15.mha'
MOVING_IMAGE_PATH = 'Input/rf_post15.mha'

if not os.path.exists(FIXED_IMAGE_PATH):
    url = 'https://data.kitware.com/api/v1/file/58ee3b778d777f16d095fd8a/download'
    urlretrieve(url, FIXED_IMAGE_PATH)
    
if not os.path.exists(MOVING_IMAGE_PATH):
    url = 'https://data.kitware.com/api/v1/file/58ee3b758d777f16d095fd87/download'
    urlretrieve(url, MOVING_IMAGE_PATH)

fixed_image = itk.imread(FIXED_IMAGE_PATH, itk.D)
moving_image = itk.imread(MOVING_IMAGE_PATH, itk.D)

dimension = fixed_image.GetImageDimension()
displacement_image_type = itk.Image[itk.Vector[itk.D,dimension],dimension]

block_radius_calculator = itk.Ultrasound.MultiResolutionFixedBlockRadiusCalculator[type(fixed_image)].New()
block_radius_calculator.SetRadius([12,4])

# Create schedule for iterative registration
search_region_source = itk.Ultrasound.MultiResolutionFixedSearchRegionImageSource[type(fixed_image), 
                                                                                  type(fixed_image),
                                                                                  displacement_image_type].New()
pyramid_schedule = itk.Array2D[itk.UI]()
pyramid_schedule.SetSize(3,2)
pyramid_schedule.SetElement(0,0,3)
pyramid_schedule.SetElement(0,1,2)
pyramid_schedule.SetElement(1,0,2)
pyramid_schedule.SetElement(1,1,1)
pyramid_schedule.SetElement(2,0,1)
pyramid_schedule.SetElement(2,1,1)

search_region_source.SetPyramidSchedule(pyramid_schedule)
search_region_source.SetSearchRegionRadiusSchedule([50,6])
search_region_source.SetOverlapSchedule(1.0)

metric_image_filter = itk.Ultrasound.NormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter[type(fixed_image),
                                                                                                     type(fixed_image),
                                                                                                     type(fixed_image)].New()

level_registration_method = itk.Ultrasound.ImageRegistrationMethod[type(fixed_image),
                                                                   type(fixed_image),
                                                                   type(fixed_image),
                                                                   displacement_image_type,
                                                                   itk.D].New()
level_registration_method.SetMetricImageFilter(metric_image_filter)

# Set up the multi-resolution registration object
multi_res_registration_method = itk.Ultrasound.MultiResolutionImageRegistrationMethod[type(fixed_image),
                                                                   type(fixed_image),
                                                                   type(fixed_image),
                                                                   displacement_image_type,
                                                                   itk.D].New()
multi_res_registration_method.SetFixedImage(fixed_image)
multi_res_registration_method.SetMovingImage(moving_image)
multi_res_registration_method.SetBlockRadiusCalculator(block_radius_calculator)
multi_res_registration_method.SetSearchRegionImageSource(search_region_source)
multi_res_registration_method.SetSchedules(pyramid_schedule, pyramid_schedule)
multi_res_registration_method.SetImageRegistrationMethod(level_registration_method)

# Run the actual registration
multi_res_registration_method.Update()

# Write out results
itk.imwrite(multi_res_registration_method.GetOutput(), 'Output/rf_post15_registered.mha')
