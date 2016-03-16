itk_module_test()

if( ITK_USE_FFTWF OR ITK_USE_FFTWD )
  set( USE_FFTW 1 )
else()
  set( USE_FFTW 0 )
endif()

set( UltrasoundTests
  itkAnalyticSignalImageFilterTest.cxx
  itkBModeImageFilterTestTiming.cxx
  itkCurvilinearArraySpecialCoordinatesImageTest.cxx
  itkFFT1DImageFilterTest.cxx
  itkComplexToComplex1DFFTImageFilterTest.cxx
  itkSpectra1DImageFilterTest.cxx
  itkSpectra1DSupportWindowImageFilterTest.cxx
  itkSpectra1DSupportWindowToMaskImageFilterTest.cxx
  itkTimeGainCompensationImageFilterTest.cxx
  itkVnlFFT1DComplexConjugateToRealImageFilterTest.cxx
  itkVnlComplexToComplex1DFFTImageFilterTest.cxx
  itkVnlForward1DFFTImageFilterTest.cxx
  )
if( USE_FFTW )
  list( APPEND UltrasoundTests
    itkFFTWForward1DFFTImageFilterTest.cxx
    itkForward1DFFTImageFilterTest.cxx
    itkFFTW1DComplexConjugateToRealImageFilterTest.cxx
    itkFFT1DComplexConjugateToRealImageFilterTest.cxx
    itkFFTWComplexToComplex1DFFTImageFilterTest.cxx
    )
endif()

CreateTestDriver( Ultrasound "${Ultrasound-Test_LIBRARIES}" "${UltrasoundTests}" )

itk_add_test( NAME itkAnalyticSignalImageFilterTest
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
itk_add_test( NAME itkBModeImageFilterTestTiming
  COMMAND UltrasoundTestDriver
  --compare
    ${ITK_TEST_OUTPUT_DIR}/itkBModeImageFilterTestTiming.mha
    DATA{Baseline/itkBModeImageFilterTestTimingBaseline.mha}
  itkBModeImageFilterTestTiming
    DATA{Input/uniform_phantom_8.9_MHz.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkBModeImageFilterTestTiming.mha
    )
itk_add_test( NAME itkCurvilinearArraySpecialCoordinatesImageTest1
  COMMAND UltrasoundTestDriver
  --compare
    ${ITK_TEST_OUTPUT_DIR}/itkCurvilinearArraySpecialCoordinatesImageTest1.mha
    DATA{Baseline/itkCurvilinearArraySpecialCoordinatesImageTest1.mha}
  itkCurvilinearArraySpecialCoordinatesImageTest
    DATA{Input/curvilinear_envelope.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkCurvilinearArraySpecialCoordinatesImageTest1.mha
    )
itk_add_test( NAME itkCurvilinearArraySpecialCoordinatesImageTest2
  COMMAND UltrasoundTestDriver
  --compare
    ${ITK_TEST_OUTPUT_DIR}/itkCurvilinearArraySpecialCoordinatesImageTest2.mha
    DATA{Baseline/itkCurvilinearArraySpecialCoordinatesImageTest2.mha}
  itkCurvilinearArraySpecialCoordinatesImageTest
    DATA{Input/curvilinear_envelope_multiframe.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkCurvilinearArraySpecialCoordinatesImageTest2.mha
    )
itk_add_test( NAME itkComplexToComplex1DFFTImageFilterTest
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
itk_add_test( NAME itkFFT1DImageFilterTest
  COMMAND UltrasoundTestDriver
  --compare
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
    ${ITK_TEST_OUTPUT_DIR}/itkFFT1DImageFilterTestOutput.mha
  itkFFT1DImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
    ${ITK_TEST_OUTPUT_DIR}/itkFFT1DImageFilterTestOutput.mha
    )
itk_add_test( NAME itkVnlComplexToComplex1DFFTImageFilterTest
  COMMAND UltrasoundTestDriver
  --compare
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
    ${ITK_TEST_OUTPUT_DIR}/itkVnlComplexToComplex1DFFTImageFilterTestOutput.mhd
  itkVnlComplexToComplex1DFFTImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaseline
    ${ITK_TEST_OUTPUT_DIR}/itkVnlComplexToComplex1DFFTImageFilterTestOutput.mhd
    )
itk_add_test( NAME itkVnlForward1DFFTImageFilterTest
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
itk_add_test( NAME itkVnlFFT1DComplexConjugateToRealImageFilterTest
  COMMAND UltrasoundTestDriver
  --compare
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
    ${ITK_TEST_OUTPUT_DIR}/itkVnlFFT1DComplexConjugateToRealImageFilterTestOutput.mhd
  itkVnlFFT1DComplexConjugateToRealImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaseline
    ${ITK_TEST_OUTPUT_DIR}/itkVnlFFT1DComplexConjugateToRealImageFilterTestOutput.mhd
    )
itk_add_test( NAME itkSpectra1DImageFilterTest
  COMMAND UltrasoundTestDriver
  itkSpectra1DImageFilterTest
    DATA{Input/uniform_phantom_8.9_MHz.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkSpectra1DImageFilterTest.mha
    )
itk_add_test( NAME itkSpectra1DSupportWindowImageFilterTest
  COMMAND UltrasoundTestDriver
  itkSpectra1DSupportWindowImageFilterTest
    DATA{Input/uniform_phantom_8.9_MHz.mha}
    )
itk_add_test( NAME itkSpectra1DSupportWindowToMaskImageFilterTest1
  COMMAND UltrasoundTestDriver
  --compare
    ${ITK_TEST_OUTPUT_DIR}/itkSpectra1DSupportWindowToMaskImageFilterTest1.mha
    DATA{Baseline/itkSpectra1DSupportWindowToMaskImageFilterTest1.mha}
  itkSpectra1DSupportWindowToMaskImageFilterTest
    DATA{Input/uniform_phantom_8.9_MHz.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkSpectra1DSupportWindowToMaskImageFilterTest1.mha
    0 0
    )
itk_add_test( NAME itkSpectra1DSupportWindowToMaskImageFilterTest2
  COMMAND UltrasoundTestDriver
  --compare
    ${ITK_TEST_OUTPUT_DIR}/itkSpectra1DSupportWindowToMaskImageFilterTest2.mha
    DATA{Baseline/itkSpectra1DSupportWindowToMaskImageFilterTest2.mha}
  itkSpectra1DSupportWindowToMaskImageFilterTest
    DATA{Input/uniform_phantom_8.9_MHz.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkSpectra1DSupportWindowToMaskImageFilterTest2.mha
    100 100
    )
itk_add_test( NAME itkTimeGainCompensationImageFilterTest
  COMMAND UltrasoundTestDriver
  --compare
    ${ITK_TEST_OUTPUT_DIR}/itkTimeGainCompensationImageFilterTest.mha
    DATA{Baseline/itkTimeGainCompensationImageFilterTest.mha}
  itkTimeGainCompensationImageFilterTest
    DATA{Input/EyeRoundOverChuck_5.0MHz_20V.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkTimeGainCompensationImageFilterTest.mha
    )

if( USE_FFTW )
  itk_add_test( NAME itkFFTWForward1DFFTImageFilterTest
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
  itk_add_test( NAME itkForward1DFFTImageFilterTest
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
  itk_add_test( NAME itkFFTW1DComplexConjugateToRealImageFilterTest
    COMMAND UltrasoundTestDriver
    --compare
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
      ${ITK_TEST_OUTPUT_DIR}/itkFFTW1DComplexConjugateToRealImageFilterTestOutput.mha
    itkFFTW1DComplexConjugateToRealImageFilterTest
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaseline
      ${ITK_TEST_OUTPUT_DIR}/itkFFTW1DComplexConjugateToRealImageFilterTestOutput.mha
      )
  itk_add_test( NAME itkFFT1DComplexConjugateToRealImageFilterTest
    COMMAND UltrasoundTestDriver
    --compare
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
      ${ITK_TEST_OUTPUT_DIR}/itkFFT1DComplexConjugateToRealImageFilterTestOutput.mha
    itkFFT1DComplexConjugateToRealImageFilterTest
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaseline
      ${ITK_TEST_OUTPUT_DIR}/itkFFT1DComplexConjugateToRealImageFilterTestOutput.mha
      )
  itk_add_test( NAME itkFFTWComplexToComplex1DFFTImageFilterTest
    COMMAND UltrasoundTestDriver
    --compare
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
      ${ITK_TEST_OUTPUT_DIR}/itkFFTWComplexToComplex1DFFTImageFilterTestOutput.mha
    itkFFTWComplexToComplex1DFFTImageFilterTest
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkForward1DFFTImageFilterTestBaseline
      ${ITK_TEST_OUTPUT_DIR}/itkFFTWComplexToComplex1DFFTImageFilterTestOutput.mha
      )
endif()
