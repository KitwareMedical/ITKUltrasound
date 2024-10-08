itk_module_test()

set(UltrasoundTests
  itkAnalyticSignalImageFilterTest.cxx
  itkAttenuationImageFilterTest.cxx
  itkBackscatterImageFilterTest.cxx
  itkBModeImageFilterTestTiming.cxx
  itkCurvilinearArraySpecialCoordinatesImageTest.cxx
  itkCurvilinearArrayUltrasoundImageFileReaderTest.cxx
  itkFFT1DImageFilterTest.cxx
  itkHDF5BModeUltrasoundImageFileReaderTest.cxx
  itkHDF5UltrasoundImageIOTest.cxx
  itkHDF5UltrasoundImageIOCanReadITKImageTest.cxx
  itkInverseScanConvertPhasedArray3DSpecialCoordinatesImageTest.cxx
  itkReplaceNonFiniteImageFilterTest.cxx
  itkScanConvertPhasedArray3DSpecialCoordinatesImageTest.cxx
  itkSliceSeriesSpecialCoordinatesImageTest.cxx
  itkSpeckleReducingAnisotropicDiffusionImageFilterTest.cxx
  itkSpectra1DImageFilterTest.cxx
  itkSpectra1DSupportWindowImageFilterTest.cxx
  itkSpectra1DSupportWindowToMaskImageFilterTest.cxx
  itkSpectra1DAveragingImageFilterTest.cxx
  itkSpectra1DNormalizeImageFilterTest.cxx
  itkTimeGainCompensationImageFilterTest.cxx
  itkInverse1DFFTImageFilterTest.cxx
  itkComplexToComplex1DFFTImageFilterTest.cxx
  itkForward1DFFTImageFilterTest.cxx
  itkBlockMatchingNormalizedCrossCorrelationFFTMetricImageFilterTest.cxx
  itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilterTest.cxx
  itkBlockMatchingBayesianRegularizationDisplacementCalculatorTest.cxx
  itkBlockMatchingImageRegistrationMethodTest.cxx
  itkBlockMatchingMultiResolutionImageRegistrationMethodTest.cxx
  itkButterworthBandpass1DFilterTest.cxx
  itkNrrdSequenceToVideoStreamTest.cxx
  )
if(ITKUltrasound_USE_VTK)
  list(APPEND UltrasoundTests
    itkSpecialCoordinatesImageToVTKStructuredGridFilterTest.cxx
    itkSpecialCoordinatesImageToVTKStructuredGridFilterSliceSeriesTest.cxx
    )
endif()

CreateTestDriver(Ultrasound "${Ultrasound-Test_LIBRARIES}" "${UltrasoundTests}")

itk_add_test(NAME itkAnalyticSignalImageFilterTest
  COMMAND UltrasoundTestDriver
  --compare
    ${CMAKE_CURRENT_SOURCE_DIR}/Baseline/itkAnalyticSignalImageFilterReal.mhd
    ${ITK_TEST_OUTPUT_DIR}/itkAnalyticSignalImageFilterTestOutputReal.mha
  --compare
    ${CMAKE_CURRENT_SOURCE_DIR}/Baseline/itkAnalyticSignalImageFilterImaginary.mhd
    ${ITK_TEST_OUTPUT_DIR}/itkAnalyticSignalImageFilterTestOutputImaginary.mha
  itkAnalyticSignalImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
    ${ITK_TEST_OUTPUT_DIR}/itkAnalyticSignalImageFilterTestOutput
    )
itk_add_test(NAME itkAttenuationImageFilterTest
  COMMAND UltrasoundTestDriver
  --compareIntensityTolerance 0.01
  --compare
    DATA{Baseline/itkAttenuationImageFilterTestOutputBaseline.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkAttenuationImageFilterTestOutput.mha
  --compare
    DATA{Baseline/itkAttenuationImageFilterTestMaskOutput.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkAttenuationImageFilterTestMaskOutput.mha
  itkAttenuationImageFilterTest
    DATA{Input/UnfusedRF-a0-spectra-cropped.mhd,UnfusedRF-a0-spectra-cropped.raw}
    DATA{Input/UnfusedRF-a0-spectra-label_1-cropped.mhd,UnfusedRF-a0-spectra-label_1-cropped.raw}
    ${ITK_TEST_OUTPUT_DIR}/itkAttenuationImageFilterTestOutput.mha
    ${ITK_TEST_OUTPUT_DIR}/itkAttenuationImageFilterTestMaskOutput.mha
    1 # work units
    0.0 # Distance between pixels sampled for attenuation estimation; 0.0 -> inclusion length
    0   # consider negative attenuations
    0   # computation mode
    0.273958 # median attenuation
  )
itk_add_test(NAME itkAttenuationImageFilterMultithreadedTest
  COMMAND UltrasoundTestDriver
  --compareIntensityTolerance 0.01
  --compare
    DATA{Baseline/itkAttenuationImageFilterTestOutputBaseline.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkAttenuationImageFilterMultithreadedTestOutput.mha
  --compare
    DATA{Baseline/itkAttenuationImageFilterTestMaskOutput.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkAttenuationImageFilterTestMaskOutput.mha
  itkAttenuationImageFilterTest
    DATA{Input/UnfusedRF-a0-spectra-cropped.mhd,UnfusedRF-a0-spectra-cropped.raw}
    DATA{Input/UnfusedRF-a0-spectra-label_1-cropped.mhd,UnfusedRF-a0-spectra-label_1-cropped.raw}
    ${ITK_TEST_OUTPUT_DIR}/itkAttenuationImageFilterMultithreadedTestOutput.mha
    ${ITK_TEST_OUTPUT_DIR}/itkAttenuationImageFilterTestMaskOutput.mha
    5 # work units
    0.0 # Distance between pixels sampled for attenuation estimation; 0.0 -> inclusion length
    0   # consider negative attenuations
    0   # computation mode
    0.273958 # median attenuation
  )
itk_add_test(NAME itkAttenuationImageFilterFixedEstimationDistanceTest
  COMMAND UltrasoundTestDriver
  --compareIntensityTolerance 0.01
  --compare
    DATA{Baseline/itkAttenuationImageFilterTestFixedEstimationDepthBaseline.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkAttenuationImageFilterTestFixedEstimationDepthOutput.mha
  --compare
    DATA{Baseline/itkAttenuationImageFilterTestMaskOutput.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkAttenuationImageFilterTestMaskOutput.mha
  itkAttenuationImageFilterTest
    DATA{Input/UnfusedRF-a0-spectra-cropped.mhd,UnfusedRF-a0-spectra-cropped.raw}
    DATA{Input/UnfusedRF-a0-spectra-label_1-cropped.mhd,UnfusedRF-a0-spectra-label_1-cropped.raw}
    ${ITK_TEST_OUTPUT_DIR}/itkAttenuationImageFilterTestFixedEstimationDepthOutput.mha
    ${ITK_TEST_OUTPUT_DIR}/itkAttenuationImageFilterTestMaskOutput.mha
    5   # work units
    3.0 # Distance between pixels sampled for attenuation estimation
    0   # consider negative attenuations
    0   # computation mode
    0.273958 # median attenuation
  )
itk_add_test(NAME itkAttenuationImageFilterFixedDistanceModeTest
  COMMAND UltrasoundTestDriver
  --compareIntensityTolerance 0.01
  --compare
    DATA{Baseline/itkAttenuationImageFilterFixedDistanceModeTestOutput.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkAttenuationImageFilterFixedDistanceModeTestOutput.mha
  --compare
    DATA{Baseline/itkAttenuationImageFilterFixedDistanceModeTestMask.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkAttenuationImageFilterFixedDistanceModeTestMask.mha
  itkAttenuationImageFilterTest
    DATA{Input/UnfusedRF-a0-spectra-cropped.mhd,UnfusedRF-a0-spectra-cropped.raw}
    DATA{Input/UnfusedRF-a0-spectra-label_1-cropped.mhd,UnfusedRF-a0-spectra-label_1-cropped.raw}
    ${ITK_TEST_OUTPUT_DIR}/itkAttenuationImageFilterFixedDistanceModeTestOutput.mha
    ${ITK_TEST_OUTPUT_DIR}/itkAttenuationImageFilterFixedDistanceModeTestMask.mha
    600 # work units
    5.0 # Distance between pixels sampled for attenuation estimation
    1   # consider negative attenuations
    2   # computation mode
    -3.72332 # median attenuation
  )
itk_add_test(NAME itkAttenuationImageFilterPairModeTest
  COMMAND UltrasoundTestDriver
  --compareIntensityTolerance 0.01
  --compare
    DATA{Baseline/itkAttenuationImageFilterPairModeTestOutput.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkAttenuationImageFilterPairModeTestOutput.mha
  --compare
    DATA{Baseline/itkAttenuationImageFilterPairModeTestMask.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkAttenuationImageFilterPairModeTestMask.mha
  itkAttenuationImageFilterTest
    DATA{Input/UnfusedRF-a0-spectra-cropped.mhd,UnfusedRF-a0-spectra-cropped.raw}
    DATA{Input/UnfusedRF-a0-spectra-label_1-cropped.mhd,UnfusedRF-a0-spectra-label_1-cropped.raw}
    ${ITK_TEST_OUTPUT_DIR}/itkAttenuationImageFilterPairModeTestOutput.mha
    ${ITK_TEST_OUTPUT_DIR}/itkAttenuationImageFilterPairModeTestMask.mha
    100 # work units
    0.0 # Distance between pixels sampled for attenuation estimation
    1   # consider negative attenuations
    1   # computation mode
    0.0693445 # median attenuation
  )
itk_add_test(NAME itkAttenuationImageFilterWholeImageTest
  COMMAND UltrasoundTestDriver
  --compareIntensityTolerance 0.01
  --compare
    DATA{Baseline/itkAttenuationImageFilterWholeImageTestOutput.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkAttenuationImageFilterWholeImageTestOutput.mha
  --compare
    DATA{Baseline/itkAttenuationImageFilterWholeImageTestMask.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkAttenuationImageFilterWholeImageTestMask.mha
  itkAttenuationImageFilterTest
    DATA{Input/UnfusedRF-a0-spectra-cropped.mhd,UnfusedRF-a0-spectra-cropped.raw}
    nul
    ${ITK_TEST_OUTPUT_DIR}/itkAttenuationImageFilterWholeImageTestOutput.mha
    ${ITK_TEST_OUTPUT_DIR}/itkAttenuationImageFilterWholeImageTestMask.mha
    100 # work units
    0.0 # Distance between pixels sampled for attenuation estimation
    1   # consider negative attenuations
    0   # computation mode
    0.458781 # median attenuation
  )
itk_add_test(NAME itkBackscatterImageFilterTest100
  COMMAND UltrasoundTestDriver
  --compareIntensityTolerance 0.01
  --compare
    DATA{Baseline/itkBackscatterImageFilterTest0Output.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkBackscatterImageFilterTest0_100Output.mha
  --compare
    DATA{Baseline/itkBackscatterImageFilterTest1Output.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkBackscatterImageFilterTest1_100Output.mha
  --compare
    DATA{Baseline/itkBackscatterImageFilterTest2Output.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkBackscatterImageFilterTest2_100Output.mha
  itkBackscatterImageFilterTest
    DATA{Input/UnfusedRF-a0-spectra-cropped.mhd,UnfusedRF-a0-spectra-cropped.raw}
    ${ITK_TEST_OUTPUT_DIR}/itkBackscatterImageFilterTest0_100Output.mha
    ${ITK_TEST_OUTPUT_DIR}/itkBackscatterImageFilterTest1_100Output.mha
    ${ITK_TEST_OUTPUT_DIR}/itkBackscatterImageFilterTest2_100Output.mha
    100 # work units
  )
itk_add_test(NAME itkBackscatterImageFilterTest16
  COMMAND UltrasoundTestDriver
  --compareIntensityTolerance 0.01
  --compare
    DATA{Baseline/itkBackscatterImageFilterTest0Output.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkBackscatterImageFilterTest0_16Output.mha
  --compare
    DATA{Baseline/itkBackscatterImageFilterTest1Output.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkBackscatterImageFilterTest1_16Output.mha
  --compare
    DATA{Baseline/itkBackscatterImageFilterTest2Output.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkBackscatterImageFilterTest2_16Output.mha
  itkBackscatterImageFilterTest
    DATA{Input/UnfusedRF-a0-spectra-cropped.mhd,UnfusedRF-a0-spectra-cropped.raw}
    ${ITK_TEST_OUTPUT_DIR}/itkBackscatterImageFilterTest0_16Output.mha
    ${ITK_TEST_OUTPUT_DIR}/itkBackscatterImageFilterTest1_16Output.mha
    ${ITK_TEST_OUTPUT_DIR}/itkBackscatterImageFilterTest2_16Output.mha
    16 # work units
  )
itk_add_test(NAME itkBModeImageFilterTestTiming
  COMMAND UltrasoundTestDriver
  --compare
    DATA{Baseline/itkBModeImageFilterTestTimingBaseline.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkBModeImageFilterTestTiming.mha
  itkBModeImageFilterTestTiming
    DATA{Input/uniform_phantom_8.9_MHz.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkBModeImageFilterTestTiming.mha
    )
itk_add_test(NAME itkCurvilinearArraySpecialCoordinatesImageTest1
  COMMAND UltrasoundTestDriver
  --compare
    DATA{Baseline/itkCurvilinearArraySpecialCoordinatesImageTest1.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkCurvilinearArraySpecialCoordinatesImageTest1.mha
  itkCurvilinearArraySpecialCoordinatesImageTest
    DATA{Input/curvilinear_envelope.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkCurvilinearArraySpecialCoordinatesImageTest1.mha
    )
itk_add_test(NAME itkCurvilinearArraySpecialCoordinatesImageTest2
  COMMAND UltrasoundTestDriver
  --compare
    DATA{Baseline/itkCurvilinearArraySpecialCoordinatesImageTest2.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkCurvilinearArraySpecialCoordinatesImageTest2.mha
  itkCurvilinearArraySpecialCoordinatesImageTest
    DATA{Input/curvilinear_envelope_multiframe.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkCurvilinearArraySpecialCoordinatesImageTest2.mha
    )
itk_add_test(NAME itkCurvilinearArrayUltrasoundImageFileReaderTest
  COMMAND UltrasoundTestDriver
  itkCurvilinearArrayUltrasoundImageFileReaderTest
    DATA{Input/curvilinear_envelope_with_metadata.mha}
    )

itk_add_test(NAME itkSliceSeriesSpecialCoordinatesImageTest
  COMMAND UltrasoundTestDriver
  --compareIntensityTolerance 10.0
  --compare
    DATA{Baseline/itkSliceSeriesSpecialCoordinatesImageTest.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkSliceSeriesSpecialCoordinatesImageTest.mha
  itkSliceSeriesSpecialCoordinatesImageTest
    DATA{Input/curvilinear_envelope_multiframe.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkSliceSeriesSpecialCoordinatesImageTest.mha
    )
itk_add_test(NAME itkSpectra1DImageFilterTest
  COMMAND UltrasoundTestDriver
  itkSpectra1DImageFilterTest
    DATA{Input/rf_voltage_15_freq_0005000000_2017-5-31_12-36-44.nrrd}
    DATA{Input/rf_voltage_15_freq_0005000000_2017-5-31_12-36-44_ReferenceSpectrum_side_lines_03_fft1d_size_128.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkSpectra1DImageFilterTest.mha
    )
itk_add_test(NAME itkSpectra1DSupportWindowImageFilterTest
  COMMAND UltrasoundTestDriver
  itkSpectra1DSupportWindowImageFilterTest
    DATA{Input/uniform_phantom_8.9_MHz.mha}
    )
itk_add_test(NAME itkSpectra1DSupportWindowToMaskImageFilterTest1
  COMMAND UltrasoundTestDriver
  --compare
    DATA{Baseline/itkSpectra1DSupportWindowToMaskImageFilterTest1.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkSpectra1DSupportWindowToMaskImageFilterTest1.mha
  itkSpectra1DSupportWindowToMaskImageFilterTest
    DATA{Input/uniform_phantom_8.9_MHz.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkSpectra1DSupportWindowToMaskImageFilterTest1.mha
    0 0
    )
itk_add_test(NAME itkSpectra1DSupportWindowToMaskImageFilterTest2
  COMMAND UltrasoundTestDriver
  --compare
    DATA{Baseline/itkSpectra1DSupportWindowToMaskImageFilterTest2.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkSpectra1DSupportWindowToMaskImageFilterTest2.mha
  itkSpectra1DSupportWindowToMaskImageFilterTest
    DATA{Input/uniform_phantom_8.9_MHz.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkSpectra1DSupportWindowToMaskImageFilterTest2.mha
    100 100
    )
itk_add_test(NAME itkSpectra1DAveragingImageFilterTest
  COMMAND UltrasoundTestDriver
    --compare
    DATA{Input/spectra_pixel_image.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkSpectra1DAveragingImageFilterTest.mha
  itkSpectra1DAveragingImageFilterTest
    ${ITK_TEST_OUTPUT_DIR}/itkSpectra1DAveragingImageFilterTest.mha
    DATA{Input/spectra_pixel_image.mha}
    DATA{Input/spectra_pixel_image.mha}
    DATA{Input/spectra_pixel_image.mha}
    )
itk_add_test(NAME itkSpectra1DAveragingImageFilterTest2
  COMMAND UltrasoundTestDriver
    --compare
    DATA{Baseline/itkSpectra1DAveragingImageFilterTest2Baseline.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkSpectra1DAveragingImageFilterTest2.mha
  itkSpectra1DAveragingImageFilterTest
    ${ITK_TEST_OUTPUT_DIR}/itkSpectra1DAveragingImageFilterTest2.mha
    DATA{Input/UnfusedRF-a0-spectra-cropped.mhd,UnfusedRF-a0-spectra-cropped.raw}
    DATA{Input/UnfusedRF-a0-spectra-cropped.mhd,UnfusedRF-a0-spectra-cropped.raw}
    )
itk_add_test(NAME itkSpectra1DNormalizeImageFilterTest
  COMMAND UltrasoundTestDriver
    --compare
    DATA{Baseline/itkSpectra1DNormalizeImageFilterTestBaseline.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkSpectra1DNormalizeImageFilterTest.mha
  itkSpectra1DNormalizeImageFilterTest
    DATA{Input/spectra_pixel_image.mha}
    DATA{Input/spectra_pixel_image.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkSpectra1DNormalizeImageFilterTest.mha
    )
itk_add_test(NAME itkSpectra1DNormalizeImageFilterTest2
  COMMAND UltrasoundTestDriver
    --compare
    DATA{Baseline/itkSpectra1DNormalizeImageFilterTest2Baseline.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkSpectra1DNormalizeImageFilterTest2.mha
  itkSpectra1DNormalizeImageFilterTest
    DATA{Input/UnfusedRF-a0-spectra-cropped.mhd,UnfusedRF-a0-spectra-cropped.raw}
    DATA{Baseline/itkSpectra1DAveragingImageFilterTest2Baseline.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkSpectra1DNormalizeImageFilterTest2.mha
    )
itk_add_test(NAME itkTimeGainCompensationImageFilterTest
  COMMAND UltrasoundTestDriver
  --compare
    DATA{Baseline/itkTimeGainCompensationImageFilterTest.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkTimeGainCompensationImageFilterTest.mha
  itkTimeGainCompensationImageFilterTest
    DATA{Input/EyeRoundOverChuck_5.0MHz_20V.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkTimeGainCompensationImageFilterTest.mha
    )
itk_add_test(NAME itkInverseScanConvertPhasedArray3DSpecialCoordinatesImageTest
  COMMAND UltrasoundTestDriver
  --compare
    DATA{Baseline/itkInverseScanConvertPhasedArray3DSpecialCoordinatesImageTestOutput.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkInverseScanConvertPhasedArray3DSpecialCoordinatesImageTestOutput.mha
  itkInverseScanConvertPhasedArray3DSpecialCoordinatesImageTest
    ${ITK_TEST_OUTPUT_DIR}/itkInverseScanConvertPhasedArray3DSpecialCoordinatesImageTestOutput.mha
    )
itk_add_test(NAME itkScanConvertPhasedArray3DSpecialCoordinatesImageTest
  COMMAND UltrasoundTestDriver
  --compare
    DATA{Baseline/itkScanConvertPhasedArray3DSpecialCoordinatesImageTestOutput.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkScanConvertPhasedArray3DSpecialCoordinatesImageTestOutput.mha
  itkScanConvertPhasedArray3DSpecialCoordinatesImageTest
    DATA{Input/itkInverseScanConvertPhasedArray3DSpecialCoordinatesImageTestOutput.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkScanConvertPhasedArray3DSpecialCoordinatesImageTestOutput.mha
    )
itk_add_test(NAME itkHDF5UltrasoundImageIOTest
  COMMAND UltrasoundTestDriver
  itkHDF5UltrasoundImageIOTest
    DATA{Input/bmode_p59.hdf5}
    )
itk_add_test(NAME itkHDF5UltrasoundImageIOCanReadITKImageTest
  COMMAND UltrasoundTestDriver
  itkHDF5UltrasoundImageIOCanReadITKImageTest
    DATA{Input/ITKImage.hdf5}
    )
itk_add_test(NAME itkHDF5BModeUltrasoundImageFileReaderTest
  COMMAND UltrasoundTestDriver
  itkHDF5BModeUltrasoundImageFileReaderTest
    DATA{Input/bmode_p59.hdf5}
  )
itk_add_test(NAME itkReplaceNonFiniteImageFilterTest
  COMMAND UltrasoundTestDriver
  --compare
    DATA{Baseline/itkReplaceNonFiniteImageFilterTestOutput.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkReplaceNonFiniteImageFilterTestOutput.mha
  itkReplaceNonFiniteImageFilterTest
    ${ITK_TEST_OUTPUT_DIR}/itkReplaceNonFiniteImageFilterTestOutput.mha
  )
itk_add_test(NAME itkBlockMatchingNormalizedCrossCorrelationFFTMetricImageFilterTest
  COMMAND UltrasoundTestDriver
  --compare
    DATA{Baseline/itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilterTestBaseline.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkBlockMatchingNormalizedCrossCorrelationFFTMetricImageFilterTestOutput.mha
  itkBlockMatchingNormalizedCrossCorrelationFFTMetricImageFilterTest
    DATA{Input/rf_pre15.mha}
    DATA{Input/rf_post15.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkBlockMatchingNormalizedCrossCorrelationFFTMetricImageFilterTestOutput.mha
  )
itk_add_test(NAME itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilterTest
  COMMAND UltrasoundTestDriver
  --compare
    DATA{Baseline/itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilterTestBaseline.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilterTestOutput.mha
  itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilterTest
    DATA{Input/rf_pre15.mha}
    DATA{Input/rf_post15.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilterTestOutput.mha
  )
itk_add_test(NAME itkBlockMatchingBayesianRegularizationDisplacementCalculatorTest
  COMMAND UltrasoundTestDriver
  --compare
    DATA{Baseline/itkBlockMatchingBayesianRegularizationDisplacementCalculatorTestBaselineComponent0.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkBlockMatchingBayesianRegularizationDisplacementCalculatorTestOutputComponent0.mha
  --compare
    DATA{Baseline/itkBlockMatchingBayesianRegularizationDisplacementCalculatorTestBaselineComponent1.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkBlockMatchingBayesianRegularizationDisplacementCalculatorTestOutputComponent1.mha
  itkBlockMatchingBayesianRegularizationDisplacementCalculatorTest
    DATA{Input/rf_pre15.mha}
    DATA{Input/rf_post15.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkBlockMatchingBayesianRegularizationDisplacementCalculatorTestOutput
  )
itk_add_test(NAME itkBlockMatchingImageRegistrationMethodTest
  COMMAND UltrasoundTestDriver
  --compare
    DATA{Baseline/itkBlockMatchingImageRegistrationMethodTestBaseline.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkBlockMatchingImageRegistrationMethodTestOutput.mha
  itkBlockMatchingImageRegistrationMethodTest
    DATA{Input/rf_pre15.mha}
    DATA{Input/rf_post15.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkBlockMatchingImageRegistrationMethodTestOutput.mha
  )
itk_add_test(NAME itkBlockMatchingMultiResolutionImageRegistrationMethodTest
  COMMAND UltrasoundTestDriver
  --compare
    DATA{Baseline/itkBlockMatchingMultiResolutionImageRegistrationMethodTestBaseline.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkBlockMatchingMultiResolutionImageRegistrationMethodTestOutput.mha
  itkBlockMatchingMultiResolutionImageRegistrationMethodTest
    DATA{Input/rf_pre15.mha}
    DATA{Input/rf_post15.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkBlockMatchingMultiResolutionImageRegistrationMethodTestOutput.mha
  )
itk_add_test(NAME itkButterworthBandpass1DFilterTest
  COMMAND UltrasoundTestDriver
  --compare
    DATA{Baseline/itkButterworthBandpass1DFilterTestOutputBaseline.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkButterworthBandpass1DFilterTestOutput.mha
  itkButterworthBandpass1DFilterTest
    DATA{Input/rf_voltage_15_freq_0005000000_2017-5-31_12-36-44.nrrd}
    ${ITK_TEST_OUTPUT_DIR}/itkButterworthBandpass1DFilterTestOutput.mha
  )
    
itk_add_test(NAME itkSpeckleReducingAnisotropicDiffusionImageFilterTest
  COMMAND UltrasoundTestDriver
  --compareIntensityTolerance 0.01
  --compare
  DATA{Baseline/itkSpeckleReducingAnisotropicDiffusionImageFilterTestOutputBaseline.mha}
  ${ITK_TEST_OUTPUT_DIR}/itkSpeckleReducingAnisotropicDiffusionImageFilterTestOutput.mha
  itkSpeckleReducingAnisotropicDiffusionImageFilterTest
  DATA{Input/PhantomRFFrame0.mha}
  ${ITK_TEST_OUTPUT_DIR}/itkSpeckleReducingAnisotropicDiffusionImageFilterTestOutput.mha
  )

itk_add_test(NAME itkNrrdSequenceToVideoStreamTest
  COMMAND UltrasoundTestDriver
  --compare
  DATA{Baseline/itkNrrdSequenceToVideoStreamTestOutputBaseline.mha}
  ${ITK_TEST_OUTPUT_DIR}/itkNrrdSequenceToVideoStreamTestOutput.mha
  itkNrrdSequenceToVideoStreamTest
  DATA{Input/SeqDemoBMode.seq.nrrd}
  ${ITK_TEST_OUTPUT_DIR}/itkNrrdSequenceToVideoStreamTestOutput.mha)

if(ITKUltrasound_USE_clFFT)
  itk_add_test(NAME itkOpenCLForward1DFFTImageFilterTest
    COMMAND UltrasoundTestDriver
    --compare
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaselineRealFull.mhd
      ${ITK_TEST_OUTPUT_DIR}/itkOpenCLForward1DFFTImageFilterTestOutputReal.mha
    --compare
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaselineImaginaryFull.mhd
      ${ITK_TEST_OUTPUT_DIR}/itkOpenCLForward1DFFTImageFilterTestOutputImaginary.mha
    itkForward1DFFTImageFilterTest
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
      ${ITK_TEST_OUTPUT_DIR}/itkOpenCLForward1DFFTImageFilterTestOutput
      3
      )
  itk_add_test(NAME itkOpenCLInverse1DFFTImageFilterTest
    COMMAND UltrasoundTestDriver
    --compare
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
      ${ITK_TEST_OUTPUT_DIR}/itkOpenCLInverse1DFFTImageFilterTestOutput.mha
    itkInverse1DFFTImageFilterTest
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaseline
      ${ITK_TEST_OUTPUT_DIR}/itkOpenCLInverse1DFFTImageFilterTestOutput.mha
      3
      )
  itk_add_test(NAME itkOpenCLComplexToComplex1DFFTImageFilterTest
    COMMAND UltrasoundTestDriver
    --compare
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
      ${ITK_TEST_OUTPUT_DIR}/itkOpenCLComplexToComplex1DFFTImageFilterTestOutput.mha
    itkComplexToComplex1DFFTImageFilterTest
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaseline
      ${ITK_TEST_OUTPUT_DIR}/itkOpenCLComplexToComplex1DFFTImageFilterTestOutput.mha
      3
      )
  itk_add_test(NAME itkOpenCLFFT1DImageFilterTest
    COMMAND UltrasoundTestDriver
    --compare
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
      ${ITK_TEST_OUTPUT_DIR}/itkOpenCLFFT1DImageFilterTestOutput.mha
    itkFFT1DImageFilterTest
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
      ${ITK_TEST_OUTPUT_DIR}/itkOpenCLFFT1DImageFilterTestOutput.mha
      3
      )
endif()

if(ITKUltrasound_USE_VTK)
  itk_add_test(NAME itkSpecialCoordinatesImageToVTKStructuredGridFilterTest
    COMMAND UltrasoundTestDriver
    itkSpecialCoordinatesImageToVTKStructuredGridFilterTest
      DATA{Input/curvilinear_envelope_with_metadata.mha}
      ${ITK_TEST_OUTPUT_DIR}/itkSpecialCoordinatesImageToVTKStructuredGridTestOutput.vtk
      )
    itk_add_test(NAME itkSpecialCoordinatesImageToVTKStructuredGridFilterSliceSeriesTest
    COMMAND UltrasoundTestDriver
    itkSpecialCoordinatesImageToVTKStructuredGridFilterSliceSeriesTest
      DATA{Input/bmode_p59.hdf5}
      ${ITK_TEST_OUTPUT_DIR}/itkSpecialCoordinatesImageToVTKStructuredGridSliceSeriesTestOutput.vtk
      )
endif()
