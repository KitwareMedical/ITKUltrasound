itk_module_test()

if(ITK_USE_FFTWF OR ITK_USE_FFTWD)
  set(use_fftw 1)
else()
  set(use_fftw 0)
endif()

set(UltrasoundTests
  itkAnalyticSignalImageFilterTest.cxx
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
  itkSpectra1DImageFilterTest.cxx
  itkSpectra1DSupportWindowImageFilterTest.cxx
  itkSpectra1DSupportWindowToMaskImageFilterTest.cxx
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
itk_add_test(NAME itkComplexToComplex1DFFTImageFilterTest
  COMMAND UltrasoundTestDriver
  --compare
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
    ${ITK_TEST_OUTPUT_DIR}/itkComplexToComplex1DFFTImageFilterTestOutput.mha
  itkComplexToComplex1DFFTImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaseline
    ${ITK_TEST_OUTPUT_DIR}/itkComplexToComplex1DFFTImageFilterTestOutput.mha
    )
# use a forward fft, then inverse fft for this instead because of the FullMatrix
# issue
itk_add_test(NAME itkFFT1DImageFilterTest
  COMMAND UltrasoundTestDriver
  --compare
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
    ${ITK_TEST_OUTPUT_DIR}/itkFFT1DImageFilterTestOutput.mha
  itkFFT1DImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
    ${ITK_TEST_OUTPUT_DIR}/itkFFT1DImageFilterTestOutput.mha
    )
itk_add_test(NAME itkVnlComplexToComplex1DFFTImageFilterTest
  COMMAND UltrasoundTestDriver
  --compare
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
    ${ITK_TEST_OUTPUT_DIR}/itkVnlComplexToComplex1DFFTImageFilterTestOutput.mhd
  itkComplexToComplex1DFFTImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaseline
    ${ITK_TEST_OUTPUT_DIR}/itkVnlComplexToComplex1DFFTImageFilterTestOutput.mhd
    1
    )
itk_add_test(NAME itkVnlForward1DFFTImageFilterTest
  COMMAND UltrasoundTestDriver
  --compare
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaselineRealFull.mhd
    ${ITK_TEST_OUTPUT_DIR}/itkVnlForward1DFFTImageFilterTestOutputReal.mha
  --compare
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaselineImaginaryFull.mhd
    ${ITK_TEST_OUTPUT_DIR}/itkVnlForward1DFFTImageFilterTestOutputImaginary.mha
  itkForward1DFFTImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
    ${ITK_TEST_OUTPUT_DIR}/itkVnlForward1DFFTImageFilterTestOutput
    1
    )
itk_add_test( NAME itkVnlInverse1DFFTImageFilterTest
  COMMAND UltrasoundTestDriver
  --compare
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
    ${ITK_TEST_OUTPUT_DIR}/itkVnlInverse1DFFTImageFilterTestOutput.mhd
  itkInverse1DFFTImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaseline
    ${ITK_TEST_OUTPUT_DIR}/itkVnlInverse1DFFTImageFilterTestOutput.mhd
    1
    )
itk_add_test(NAME itkVnlFFT1DImageFilterTest
  COMMAND UltrasoundTestDriver
  --compare
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
    ${ITK_TEST_OUTPUT_DIR}/itkVnlFFT1DImageFilterTestOutput.mha
  itkFFT1DImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
    ${ITK_TEST_OUTPUT_DIR}/itkVnlFFT1DImageFilterTestOutput.mha
    1
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


itk_add_test(NAME itkForward1DFFTImageFilterTest
  COMMAND UltrasoundTestDriver
  --compare
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaselineRealFull.mhd
    ${ITK_TEST_OUTPUT_DIR}/itkForward1DFFTImageFilterTestOutputReal.mha
  --compare
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaselineImaginaryFull.mhd
    ${ITK_TEST_OUTPUT_DIR}/itkForward1DFFTImageFilterTestOutputImaginary.mha
  itkForward1DFFTImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
    ${ITK_TEST_OUTPUT_DIR}/itkForward1DFFTImageFilterTestOutput
    )
itk_add_test(NAME itkInverse1DFFTImageFilterTest
  COMMAND UltrasoundTestDriver
  --compare
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
    ${ITK_TEST_OUTPUT_DIR}/itkInverse1DFFTImageFilterTestOutput.mha
  itkInverse1DFFTImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaseline
    ${ITK_TEST_OUTPUT_DIR}/itkInverse1DFFTImageFilterTestOutput.mha
    )
    
if(use_fftw)
  itk_add_test(NAME itkFFTWForward1DFFTImageFilterTest
    COMMAND UltrasoundTestDriver
    --compare
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaselineRealFull.mhd
      ${ITK_TEST_OUTPUT_DIR}/itkFFTWForward1DFFTImageFilterTestOutputReal.mha
    --compare
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaselineImaginaryFull.mhd
      ${ITK_TEST_OUTPUT_DIR}/itkFFTWForward1DFFTImageFilterTestOutputImaginary.mha
    itkForward1DFFTImageFilterTest
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
      ${ITK_TEST_OUTPUT_DIR}/itkFFTWForward1DFFTImageFilterTestOutput
      2
      )
  itk_add_test(NAME itkFFTWInverse1DFFTImageFilterTest
    COMMAND UltrasoundTestDriver
    --compare
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
      ${ITK_TEST_OUTPUT_DIR}/itkFFTWInverse1DFFTImageFilterTestOutput.mha
    itkInverse1DFFTImageFilterTest
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaseline
      ${ITK_TEST_OUTPUT_DIR}/itkFFTWInverse1DFFTImageFilterTestOutput.mha
      2
      )
  itk_add_test(NAME itkFFTWComplexToComplex1DFFTImageFilterTest
    COMMAND UltrasoundTestDriver
    --compare
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
      ${ITK_TEST_OUTPUT_DIR}/itkFFTWComplexToComplex1DFFTImageFilterTestOutput.mha
    itkComplexToComplex1DFFTImageFilterTest
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaseline
      ${ITK_TEST_OUTPUT_DIR}/itkFFTWComplexToComplex1DFFTImageFilterTestOutput.mha
      2
      )
  itk_add_test(NAME itkFFTW1DImageFilterTest
    COMMAND UltrasoundTestDriver
    --compare
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
      ${ITK_TEST_OUTPUT_DIR}/itkFFTW1DImageFilterTestOutput.mha
    itkFFT1DImageFilterTest
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
      ${ITK_TEST_OUTPUT_DIR}/itkFFTW1DImageFilterTestOutput.mha
      2
      )
endif()

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
