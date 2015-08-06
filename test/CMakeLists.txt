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
  itkFFT1DComplexToComplexImageFilterTest.cxx
  itkSpectra1DImageFilterTest.cxx
  itkSpectra1DSupportWindowImageFilterTest.cxx
  itkSpectra1DSupportWindowToMaskImageFilterTest.cxx
  itkTimeGainCompensationImageFilterTest.cxx
  itkVnlFFT1DComplexConjugateToRealImageFilterTest.cxx
  itkVnlFFT1DComplexToComplexImageFilterTest.cxx
  itkVnlFFT1DRealToComplexConjugateImageFilterTest.cxx
  )
if( USE_FFTW )
  list( APPEND UltrasoundTests
    itkFFTW1DRealToComplexConjugateImageFilterTest.cxx
    itkFFT1DRealToComplexConjugateImageFilterTest.cxx
    itkFFTW1DComplexConjugateToRealImageFilterTest.cxx
    itkFFTW1DComplexToComplexImageFilterTest.cxx
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
itk_add_test( NAME itkCurvilinearArraySpecialCoordinatesImageTest
  COMMAND UltrasoundTestDriver
  --compare
    ${ITK_TEST_OUTPUT_DIR}/itkCurvilinearArraySpecialCoordinatesImageTest.mha
    DATA{Baseline/itkCurvilinearArraySpecialCoordinatesImageTest.mha}
  itkCurvilinearArraySpecialCoordinatesImageTest
    DATA{Input/curvilinear_envelope.mha}
    ${ITK_TEST_OUTPUT_DIR}/itkCurvilinearArraySpecialCoordinatesImageTest.mha
    )
itk_add_test( NAME itkFFT1DComplexToComplexImageFilterTest
  COMMAND UltrasoundTestDriver
  --compare
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
    ${ITK_TEST_OUTPUT_DIR}/itkFFT1DComplexToComplexImageFilterTestOutput.mha
  itkFFT1DComplexToComplexImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkFFT1DRealToComplexConjugateImageFilterTestBaseline
    ${ITK_TEST_OUTPUT_DIR}/itkFFT1DComplexToComplexImageFilterTestOutput.mha
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
itk_add_test( NAME itkVnlFFT1DComplexToComplexImageFilterTest
  COMMAND UltrasoundTestDriver
  --compare
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
    ${ITK_TEST_OUTPUT_DIR}/itkVnlFFT1DComplexToComplexImageFilterTestOutput.mhd
  itkVnlFFT1DComplexToComplexImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkFFT1DRealToComplexConjugateImageFilterTestBaseline
    ${ITK_TEST_OUTPUT_DIR}/itkVnlFFT1DComplexToComplexImageFilterTestOutput.mhd
    )
itk_add_test( NAME itkVnlFFT1DRealToComplexConjugateImageFilterTest
  COMMAND UltrasoundTestDriver
  --compare
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkFFT1DRealToComplexConjugateImageFilterTestBaselineRealFull.mhd
    ${ITK_TEST_OUTPUT_DIR}/itkVnlFFT1DRealToComplexConjugateImageFilterTestOutputReal.mha
  --compare
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkFFT1DRealToComplexConjugateImageFilterTestBaselineImaginaryFull.mhd
    ${ITK_TEST_OUTPUT_DIR}/itkVnlFFT1DRealToComplexConjugateImageFilterTestOutputImaginary.mha
  itkVnlFFT1DRealToComplexConjugateImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
    ${ITK_TEST_OUTPUT_DIR}/itkVnlFFT1DRealToComplexConjugateImageFilterTestOutput
    )
itk_add_test( NAME itkVnlFFT1DComplexConjugateToRealImageFilterTest
  COMMAND UltrasoundTestDriver
  --compare
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
    ${ITK_TEST_OUTPUT_DIR}/itkVnlFFT1DComplexConjugateToRealImageFilterTestOutput.mhd
  itkVnlFFT1DComplexConjugateToRealImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkFFT1DRealToComplexConjugateImageFilterTestBaseline
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
  itk_add_test( NAME itkFFTW1DRealToComplexConjugateImageFilterTest
    COMMAND UltrasoundTestDriver
    --compare
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkFFT1DRealToComplexConjugateImageFilterTestBaselineRealFull.mhd
      ${ITK_TEST_OUTPUT_DIR}/itkFFTW1DRealToComplexConjugateImageFilterTestOutputReal.mha
    --compare
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkFFT1DRealToComplexConjugateImageFilterTestBaselineImaginaryFull.mhd
      ${ITK_TEST_OUTPUT_DIR}/itkFFTW1DRealToComplexConjugateImageFilterTestOutputImaginary.mha
    itkFFTW1DRealToComplexConjugateImageFilterTest
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
      ${ITK_TEST_OUTPUT_DIR}/itkFFTW1DRealToComplexConjugateImageFilterTestOutput
      )
  itk_add_test( NAME itkFFT1DRealToComplexConjugateImageFilterTest
    COMMAND UltrasoundTestDriver
    --compare
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkFFT1DRealToComplexConjugateImageFilterTestBaselineRealFull.mhd
      ${ITK_TEST_OUTPUT_DIR}/itkFFT1DRealToComplexConjugateImageFilterTestOutputReal.mha
    --compare
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkFFT1DRealToComplexConjugateImageFilterTestBaselineImaginaryFull.mhd
      ${ITK_TEST_OUTPUT_DIR}/itkFFT1DRealToComplexConjugateImageFilterTestOutputImaginary.mha
    itkFFT1DRealToComplexConjugateImageFilterTest
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
      ${ITK_TEST_OUTPUT_DIR}/itkFFT1DRealToComplexConjugateImageFilterTestOutput
      )
  itk_add_test( NAME itkFFTW1DComplexConjugateToRealImageFilterTest
    COMMAND UltrasoundTestDriver
    --compare
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
      ${ITK_TEST_OUTPUT_DIR}/itkFFTW1DComplexConjugateToRealImageFilterTestOutput.mha
    itkFFTW1DComplexConjugateToRealImageFilterTest
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkFFT1DRealToComplexConjugateImageFilterTestBaseline
      ${ITK_TEST_OUTPUT_DIR}/itkFFTW1DComplexConjugateToRealImageFilterTestOutput.mha
      )
  itk_add_test( NAME itkFFTW1DComplexToComplexImageFilterTest
    COMMAND UltrasoundTestDriver
    --compare
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
      ${ITK_TEST_OUTPUT_DIR}/itkFFTW1DComplexToComplexImageFilterTestOutput.mha
    itkFFTW1DComplexToComplexImageFilterTest
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkFFT1DRealToComplexConjugateImageFilterTestBaseline
      ${ITK_TEST_OUTPUT_DIR}/itkFFTW1DComplexToComplexImageFilterTestOutput.mha
      )
endif()
