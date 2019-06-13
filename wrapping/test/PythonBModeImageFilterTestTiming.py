#==========================================================================
#
#   Copyright Insight Software Consortium
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#          http://www.apache.org/licenses/LICENSE-2.0.txt
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
#==========================================================================*/

#
#  Test the performance of the BModeImageFilter using Python
#

import itk
from sys import argv

input_filename = argv[1]
output_filename = argv[2]
PixelType = itk.F
dim = 3
ImageType = itk.Image[PixelType, dim]
inputImg = itk.imread(input_filename, PixelType)
inputImg.DisconnectPipeline()

BModeFilterType = itk.BModeImageFilter[ImageType, ImageType]
bMode = BModeFilterType.New()
bMode.SetInput(inputImg)

WindowingType = itk.IntensityWindowingImageFilter[ImageType, ImageType]
window = WindowingType.New()
window.SetInput(bMode.GetOutput())

clock = itk.TimeProbe()

runs = 1000
for i in range(runs):
    bMode.Modified()
    clock.Start()
    window.Update()
    clock.Stop()

size = bMode.GetOutput().GetLargestPossibleRegion().GetSize()
frame_rate = size[2]/clock.GetMean()
print("Frame rate achieved over {} runs was {} fp{}.".format(clock.GetNumberOfStarts(), frame_rate, clock.GetUnit()))  

itk.imwrite(window.GetOutput(), output_filename)