/*=========================================================================
 *
 *  Copyright NumFOCUS
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         https://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "itkTestingMacros.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageToVideoFilter.h"

/* Example demonstration conversion of NRRD sequence file to a video stream.
 *
 * 3D Slicer writes out video frame sequences in NRRD list format, while ITK
 * interprets NRRD lists as pixel channels. This leads to the undesirable situation
 * where video frames load into ITK with the time domain along the pixel channel
 * direction, which is not easy to work with for ITK filtering pipelines or general
 * analysis.
 *
 * This example demonstrates reading in a sequence of 2D video frames captured from
 * an ultrasound phantom as a 2D VectorImage with each pixel representing a vector
 * of measurements in time. We manually allocate an itk::Image and manually iterate over
 * the VectorImage getting pixel values to produce a 3D volume representing two spatial
 * dimensions and one temporal dimension, which we save to disk. We then convert
 * the volume to an itk::VideoStream of 2D ultrasound images in order to carry out
 * subsequent analysis with time as a first-class citizen of the data object.
 *
 * Underlying metadata (spacing, origin, directionality) and assumptions (LPS orientation)
 * are propagated through each procedural step.
 */
int
itkNrrdSequenceToVideoStreamTest(int argc, char * argv[])
{
  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0];
    std::cerr << " <inputImage.seq.nrrd> <outputFrame0.mha>";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }
  const char * inputImageFileName = argv[1];
  const char * outputFrameFileName = argv[2];

  using PixelType = float;
  const unsigned int Dimension = 3;
  using VectorImageType = itk::VectorImage<PixelType, Dimension>;

  // Read in sequence data.
  // Our example data has 3 spatial axes of size 512x512x1
  // and a temporal axis of size 118.
  // We must read in all axes to preserve metadata from the NRRD file.
  using ReaderType = itk::ImageFileReader<VectorImageType>;
  ReaderType::Pointer reader = ReaderType::New();
  ITK_TRY_EXPECT_NO_EXCEPTION(reader->SetFileName(inputImageFileName));
  ITK_TRY_EXPECT_NO_EXCEPTION(reader->Update());
  VectorImageType::Pointer input = reader->GetOutput();

  // Allocate the 3D destination volume.
  // In ITK the 0th accessor index is the fastest-moving index in memory,
  // so we will let the 0th index correspond to the pixel component (time) dimension
  // as the remaining 1,..,n indices correspond to the 0,...,n-1 spatial indices
  // in the input image.
  using ImageType = itk::Image<PixelType, Dimension>;
  typename ImageType::RegionType region;
  region.SetSize(itk::MakeSize(input->GetNumberOfComponentsPerPixel(),
                               input->GetLargestPossibleRegion().GetSize(0),
                               input->GetLargestPossibleRegion().GetSize(1)));
  region.SetIndex(itk::MakeIndex(0, 0, 0));

  ImageType::Pointer scalarImage = ImageType::New();
  ITK_TRY_EXPECT_NO_EXCEPTION(scalarImage->SetRegions(region));
  ITK_TRY_EXPECT_NO_EXCEPTION(scalarImage->Allocate());

  // Propagate metadata.
  // The Image class requires isotropic spacing, though this may be
  // unrealistic for video data acquired with a variable frame rate.
  // Here we specify that our data is collected at approximately 1 frame / 30 ms.
  ImageType::SpacingType imageSpacing;
  imageSpacing[0] = 0.03f;
  imageSpacing[1] = input->GetSpacing()[0];
  imageSpacing[2] = input->GetSpacing()[1];
  scalarImage->SetSpacing(imageSpacing);

  ImageType::PointType imageOrigin;
  imageOrigin[0] = 0.0f;
  imageOrigin[1] = input->GetOrigin()[0];
  imageOrigin[2] = input->GetOrigin()[1];
  ITK_TRY_EXPECT_NO_EXCEPTION(scalarImage->SetOrigin(imageOrigin));

  // We will never rotate along the time axis, so we set the identity
  // direction for time and propagate the remaining direction submatrix.
  // We intentionally discard the elevational axis.
  ImageType::DirectionType imageDirection;
  imageDirection.SetIdentity();
  for (itk::IndexValueType r = 0; r < input->GetImageDimension() - 1; ++r)
  {
    for (itk::IndexValueType c = 0; c < input->GetImageDimension() - 1; ++c)
    {
      imageDirection(r + 1, c + 1) = input->GetDirection()(r, c);
    }
  }
  ITK_TRY_EXPECT_NO_EXCEPTION(scalarImage->SetDirection(imageDirection));

  // Iteratively copy input channel values to the corresponding pixel in the 3D volume.
  // In this process the elevational axis of size 1 is intentionally discarded so that
  // we are left with axes of size > 1.
  using ScalarIteratorType = typename itk::ImageRegionIteratorWithIndex<ImageType>;
  ScalarIteratorType it(scalarImage, scalarImage->GetLargestPossibleRegion());
  it.GoToBegin();
  while (!it.IsAtEnd())
  {
    auto scalarIndex = it.GetIndex();
    it.Set(input->GetPixel(itk::MakeIndex(scalarIndex[1], scalarIndex[2], 0)).GetElement(scalarIndex[0]));
    ++it;
  }
  ITK_TEST_EXPECT_EQUAL(scalarImage->GetPixel(itk::MakeIndex(0, 0, 0)),
                       input->GetPixel(itk::MakeIndex(0, 0, 0)).GetElement(0));

  // Convert the 3D volume to a VideoStream of 2D frames for analysis.
  using VideoFilterType = typename itk::ImageToVideoFilter<ImageType>;
  VideoFilterType::Pointer videoFilter = VideoFilterType::New();
  videoFilter->SetInput(scalarImage);
  videoFilter->SetFrameAxis(0);
  ITK_TRY_EXPECT_NO_EXCEPTION(videoFilter->Update());

  for (size_t axis = 0; axis < VideoFilterType::OutputFrameType::ImageDimension; ++axis)
  {
    ITK_TEST_EXPECT_EQUAL(videoFilter->GetOutput()->GetFrame(0)->GetLargestPossibleRegion().GetSize()[axis],
                          input->GetLargestPossibleRegion().GetSize()[axis]);
    ITK_TEST_EXPECT_EQUAL(videoFilter->GetOutput()->GetFrame(0)->GetSpacing()[axis], input->GetSpacing()[axis]);
    ITK_TEST_EXPECT_EQUAL(videoFilter->GetOutput()->GetFrame(0)->GetOrigin()[axis], input->GetOrigin()[axis]);
    for (size_t axis2 = 0; axis2 < VideoFilterType::OutputFrameType::ImageDimension; ++axis2)
    {
      ITK_TEST_EXPECT_EQUAL(videoFilter->GetOutput()->GetFrame(0)->GetDirection()(axis,axis2), input->GetDirection()(axis,axis2));
    }
  }

  // Save the first frame to disk for verification.
  ITK_TRY_EXPECT_NO_EXCEPTION(itk::WriteImage(videoFilter->GetOutput()->GetFrame(0), outputFrameFileName));

  return EXIT_SUCCESS;
}
