itk_module_test()

if( ITK_USE_FFTWF OR ITK_USE_FFTWD )
  set( USE_FFTW 1 )
else()
  set( USE_FFTW 0 )
endif()

set( FFT1DTests
  itkAnalyticSignalImageFilterTest.cxx
  itkFFT1DImageFilterTest.cxx
  itkFFT1DComplexToComplexImageFilterTest.cxx
  itkSpectra1DImageFilterTest.cxx
  itkSpectra1DSupportWindowImageFilterTest.cxx
  itkVnlFFT1DComplexConjugateToRealImageFilterTest.cxx
  itkVnlFFT1DComplexToComplexImageFilterTest.cxx
  itkVnlFFT1DRealToComplexConjugateImageFilterTest.cxx
  )
if( USE_FFTW )
  list( APPEND FFT1DTests
    itkFFTW1DRealToComplexConjugateImageFilterTest.cxx
    itkFFTW1DComplexConjugateToRealImageFilterTest.cxx
    itkFFTW1DComplexToComplexImageFilterTest.cxx
    )
endif()

CreateTestDriver( FFT1D "${FFT1D-Test_LIBRARIES}" "${FFT1DTests}" )

itk_add_test( NAME itkAnalyticSignalImageFilterTest
  COMMAND FFT1DTestDriver
  --compare ${CMAKE_CURRENT_SOURCE_DIR}/Baseline/itkAnalyticSignalImageFilterReal.mhd
    itkAnalyticSignalImageFilterTestOutputReal.mha
  --compare ${CMAKE_CURRENT_SOURCE_DIR}/Baseline/itkAnalyticSignalImageFilterImaginary.mhd
    itkAnalyticSignalImageFilterTestOutputImaginary.mha
  itkAnalyticSignalImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
    itkAnalyticSignalImageFilterTestOutput
    )
itk_add_test( NAME itkFFT1DComplexToComplexImageFilterTest
  COMMAND FFT1DTestDriver
  --compare ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
            itkFFT1DComplexToComplexImageFilterTestOutput.mha
  itkFFT1DComplexToComplexImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkFFT1DRealToComplexConjugateImageFilterTestBaseline
    itkFFT1DComplexToComplexImageFilterTestOutput.mha
    )
# use a forward fft, then inverse fft for this instead because of the FullMatrix
# issue
itk_add_test( NAME itkFFT1DImageFilterTest
  COMMAND FFT1DTestDriver
  --compare ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
            itkFFT1DImageFilterTestOutput.mha
  itkFFT1DImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
    itkFFT1DImageFilterTestOutput.mha
    )
itk_add_test( NAME itkVnlFFT1DComplexToComplexImageFilterTest
  COMMAND FFT1DTestDriver
  --compare ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
            itkVnlFFT1DComplexToComplexImageFilterTestOutput.mhd
  itkVnlFFT1DComplexToComplexImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkFFT1DRealToComplexConjugateImageFilterTestBaseline
    itkVnlFFT1DComplexToComplexImageFilterTestOutput.mhd
    )
itk_add_test( NAME itkVnlFFT1DRealToComplexConjugateImageFilterTest
  COMMAND FFT1DTestDriver
  --compare ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkFFT1DRealToComplexConjugateImageFilterTestBaselineRealFull.mhd
            itkVnlFFT1DRealToComplexConjugateImageFilterTestOutputReal.mha
  --compare ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkFFT1DRealToComplexConjugateImageFilterTestBaselineImaginaryFull.mhd
            itkVnlFFT1DRealToComplexConjugateImageFilterTestOutputImaginary.mha
  itkVnlFFT1DRealToComplexConjugateImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
    itkVnlFFT1DRealToComplexConjugateImageFilterTestOutput
    )
itk_add_test( NAME itkVnlFFT1DComplexConjugateToRealImageFilterTest
  COMMAND FFT1DTestDriver
  --compare ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
            itkVnlFFT1DComplexConjugateToRealImageFilterTestOutput.mhd
  itkVnlFFT1DComplexConjugateToRealImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkFFT1DRealToComplexConjugateImageFilterTestBaseline
    itkVnlFFT1DComplexConjugateToRealImageFilterTestOutput.mhd
    )
itk_add_test( NAME itkSpectra1DImageFilterTest
  COMMAND FFT1DTestDriver
  itkSpectra1DImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkSpectra1DImageFilterTestBaseline
    itkSpectra1DImageFilterTest
    )
itk_add_test( NAME itkSpectra1DSupportWindowImageFilterTest
  COMMAND FFT1DTestDriver
  itkSpectra1DSupportWindowImageFilterTest
    ${CMAKE_CURRENT_SOURCE_DIR}/Input/uniform_phantom_8.9_MHz.mha
    itkSpectra1DSupportWindowImageFilterTestOutput.mha
    )

if( USE_FFTW )
  itk_add_test( NAME itkFFTW1DRealToComplexConjugateImageFilterTest
    COMMAND FFT1DTestDriver
    --compare ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkFFT1DRealToComplexConjugateImageFilterTestBaselineRealFull.mhd
              itkFFTW1DRealToComplexConjugateImageFilterTestOutputReal.mha
    --compare ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkFFT1DRealToComplexConjugateImageFilterTestBaselineImaginaryFull.mhd
              itkFFTW1DRealToComplexConjugateImageFilterTestOutputImaginary.mha
    itkFFTW1DRealToComplexConjugateImageFilterTest
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
      itkFFTW1DRealToComplexConjugateImageFilterTestOutput
      )
  itk_add_test( NAME itkFFTW1DComplexConjugateToRealImageFilterTest
    COMMAND FFT1DTestDriver
    --compare ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
              itkFFTW1DComplexConjugateToRealImageFilterTestOutput.mha
    itkFFTW1DComplexConjugateToRealImageFilterTest
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkFFT1DRealToComplexConjugateImageFilterTestBaseline
      itkFFTW1DComplexConjugateToRealImageFilterTestOutput.mha
      )
  itk_add_test( NAME itkFFTW1DComplexToComplexImageFilterTest
    COMMAND FFT1DTestDriver
    --compare ${CMAKE_CURRENT_SOURCE_DIR}/Input/TreeBarkTexture.png
              itkFFTW1DComplexToComplexImageFilterTestOutput.mha
    itkFFTW1DComplexToComplexImageFilterTest
      ${CMAKE_CURRENT_SOURCE_DIR}/Input/itkFFT1DRealToComplexConjugateImageFilterTestBaseline
      itkFFTW1DComplexToComplexImageFilterTestOutput.mha
      )
endif()
