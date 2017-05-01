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
  itkComplexToComplex1DFFTImageFilterTest.cxx
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
  itkVnlInverse1DFFTImageFilterTest.cxx
  itkVnlComplexToComplex1DFFTImageFilterTest.cxx
  itkVnlForward1DFFTImageFilterTest.cxx
  itkBlockMatchingNormalizedCrossCorrelationFFTMetricImageFilterTest.cxx
  itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilterTest.cxx
  itkBlockMatchingBayesianRegularizationDisplacementCalculatorTest.cxx
  )
if(use_fftw)
  list(APPEND UltrasoundTests
    itkFFTWForward1DFFTImageFilterTest.cxx
    itkForward1DFFTImageFilterTest.cxx
    itkFFTWInverse1DFFTImageFilterTest.cxx
    itkInverse1DFFTImageFilterTest.cxx
    itkFFTWComplexToComplex1DFFTImageFilterTest.cxx
    )
endif()
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
    ${ITK_TEST_OUTPUT_DIR}/itkBModeImageFilterTestTiming.mha
    DATA{Baseline/itkBModeImageFilterTestTimingBaseline.mha}
  itkBModeImageFilterTestTiming
    DATA{Input/uniform_phantom_8.9_MHz.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkBModeImageFilterTestTiming.mha
    )
itk_add_test(NAME itkCurvilinearArraySpecialCoordinatesImageTest1
  COMMAND UltrasoundTestDriver
  --compare
    ${ITK_TEST_OUTPUT_DIR}/itkCurvilinearArraySpecialCoordinatesImageTest1.mha
    DATA{Baseline/itkCurvilinearArraySpecialCoordinatesImageTest1.mha}
  itkCurvilinearArraySpecialCoordinatesImageTest
    DATA{Input/curvilinear_envelope.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkCurvilinearArraySpecialCoordinatesImageTest1.mha
    )
itk_add_test(NAME itkCurvilinearArraySpecialCoordinatesImageTest2
  COMMAND UltrasoundTestDriver
  --compare
    ${ITK_TEST_OUTPUT_DIR}/itkCurvilinearArraySpecialCoordinatesImageTest2.mha
    DATA{Baseline/itkCurvilinearArraySpecialCoordinatesImageTest2.mha}
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
  itkVnlComplexToComplex1DFFTImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaseline
    ${ITK_TEST_OUTPUT_DIR}/itkVnlComplexToComplex1DFFTImageFilterTestOutput.mhd
    )
itk_add_test(NAME itkVnlForward1DFFTImageFilterTest
  COMMAND UltrasoundTestDriver
  --compare
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaselineRealFull.mhd
    ${ITK_TEST_OUTPUT_DIR}/itkVnlForward1DFFTImageFilterTestOutputReal.mha
  --compare
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaselineImaginaryFull.mhd
    ${ITK_TEST_OUTPUT_DIR}/itkVnlForward1DFFTImageFilterTestOutputImaginary.mha
  itkVnlForward1DFFTImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
    ${ITK_TEST_OUTPUT_DIR}/itkVnlForward1DFFTImageFilterTestOutput
    )
itk_add_test( NAME itkVnlInverse1DFFTImageFilterTest
  COMMAND UltrasoundTestDriver
  --compare
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
    ${ITK_TEST_OUTPUT_DIR}/itkVnlInverse1DFFTImageFilterTestOutput.mhd
  itkVnlInverse1DFFTImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaseline
    ${ITK_TEST_OUTPUT_DIR}/itkVnlInverse1DFFTImageFilterTestOutput.mhd
    )
  itk_add_test(NAME itkSliceSeriesSpecialCoordinatesImageTest
  COMMAND UltrasoundTestDriver
  --compareIntensityTolerance 10.0
  --compare
    ${ITK_TEST_OUTPUT_DIR}/itkSliceSeriesSpecialCoordinatesImageTest.mha
    DATA{Baseline/itkSliceSeriesSpecialCoordinatesImageTest.mha}
  itkSliceSeriesSpecialCoordinatesImageTest
    DATA{Input/curvilinear_envelope_multiframe.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkSliceSeriesSpecialCoordinatesImageTest.mha
    )
itk_add_test(NAME itkSpectra1DImageFilterTest
  COMMAND UltrasoundTestDriver
  itkSpectra1DImageFilterTest
    DATA{Input/uniform_phantom_8.9_MHz.mha}
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
    ${ITK_TEST_OUTPUT_DIR}/itkSpectra1DSupportWindowToMaskImageFilterTest1.mha
    DATA{Baseline/itkSpectra1DSupportWindowToMaskImageFilterTest1.mha}
  itkSpectra1DSupportWindowToMaskImageFilterTest
    DATA{Input/uniform_phantom_8.9_MHz.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkSpectra1DSupportWindowToMaskImageFilterTest1.mha
    0 0
    )
itk_add_test(NAME itkSpectra1DSupportWindowToMaskImageFilterTest2
  COMMAND UltrasoundTestDriver
  --compare
    ${ITK_TEST_OUTPUT_DIR}/itkSpectra1DSupportWindowToMaskImageFilterTest2.mha
    DATA{Baseline/itkSpectra1DSupportWindowToMaskImageFilterTest2.mha}
  itkSpectra1DSupportWindowToMaskImageFilterTest
    DATA{Input/uniform_phantom_8.9_MHz.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkSpectra1DSupportWindowToMaskImageFilterTest2.mha
    100 100
    )
itk_add_test(NAME itkTimeGainCompensationImageFilterTest
  COMMAND UltrasoundTestDriver
  --compare
    ${ITK_TEST_OUTPUT_DIR}/itkTimeGainCompensationImageFilterTest.mha
    DATA{Baseline/itkTimeGainCompensationImageFilterTest.mha}
  itkTimeGainCompensationImageFilterTest
    DATA{Input/EyeRoundOverChuck_5.0MHz_20V.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkTimeGainCompensationImageFilterTest.mha
    )
itk_add_test(NAME itkInverseScanConvertPhasedArray3DSpecialCoordinatesImageTest
  COMMAND UltrasoundTestDriver
  --compare
    ${ITK_TEST_OUTPUT_DIR}/itkInverseScanConvertPhasedArray3DSpecialCoordinatesImageTestOutput.mha
    DATA{Baseline/itkInverseScanConvertPhasedArray3DSpecialCoordinatesImageTestOutput.mha}
  itkInverseScanConvertPhasedArray3DSpecialCoordinatesImageTest
    ${ITK_TEST_OUTPUT_DIR}/itkInverseScanConvertPhasedArray3DSpecialCoordinatesImageTestOutput.mha
    )
itk_add_test(NAME itkScanConvertPhasedArray3DSpecialCoordinatesImageTest
  COMMAND UltrasoundTestDriver
  --compare
    ${ITK_TEST_OUTPUT_DIR}/itkScanConvertPhasedArray3DSpecialCoordinatesImageTestOutput.mha
    DATA{Baseline/itkScanConvertPhasedArray3DSpecialCoordinatesImageTestOutput.mha}
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
    ${ITK_TEST_OUTPUT_DIR}/itkReplaceNonFiniteImageFilterTestOutput.mha
    DATA{Baseline/itkReplaceNonFiniteImageFilterTestOutput.mha}
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

if(use_fftw)
  itk_add_test(NAME itkFFTWForward1DFFTImageFilterTest
    COMMAND UltrasoundTestDriver
    --compare
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaselineRealFull.mhd
      ${ITK_TEST_OUTPUT_DIR}/itkFFTWForward1DFFTImageFilterTestOutputReal.mha
    --compare
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaselineImaginaryFull.mhd
      ${ITK_TEST_OUTPUT_DIR}/itkFFTWForward1DFFTImageFilterTestOutputImaginary.mha
    itkFFTWForward1DFFTImageFilterTest
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
      ${ITK_TEST_OUTPUT_DIR}/itkFFTWForward1DFFTImageFilterTestOutput
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
  itk_add_test(NAME itkFFTWInverse1DFFTImageFilterTest
    COMMAND UltrasoundTestDriver
    --compare
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
      ${ITK_TEST_OUTPUT_DIR}/itkFFTWInverse1DFFTImageFilterTestOutput.mha
    itkFFTWInverse1DFFTImageFilterTest
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaseline
      ${ITK_TEST_OUTPUT_DIR}/itkFFTWInverse1DFFTImageFilterTestOutput.mha
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
  itk_add_test(NAME itkFFTWComplexToComplex1DFFTImageFilterTest
    COMMAND UltrasoundTestDriver
    --compare
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
      ${ITK_TEST_OUTPUT_DIR}/itkFFTWComplexToComplex1DFFTImageFilterTestOutput.mha
    itkFFTWComplexToComplex1DFFTImageFilterTest
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaseline
      ${ITK_TEST_OUTPUT_DIR}/itkFFTWComplexToComplex1DFFTImageFilterTestOutput.mha
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
