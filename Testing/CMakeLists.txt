set( DATA_INPUT ${ITK_FFT_EXTENSIONS_SOURCE_DIR}/Testing/Data/Input )
set( DATA_BASELINE ${ITK_FFT_EXTENSIONS_SOURCE_DIR}/Testing/Data/Baseline )

set( _test_srcs
  itkAnalyticSignalImageFilterTest.cxx
  itkFFT1DImageFilterTest.cxx
  itkFFT1DComplexToComplexImageFilterTest.cxx
  itkFFTW1DComplexToComplexImageFilterTest.cxx
  itkVnlFFT1DComplexConjugateToRealImageFilterTest.cxx
  itkVnlFFT1DComplexToComplexImageFilterTest.cxx
  itkVnlFFT1DRealToComplexConjugateImageFilterTest.cxx
  )
if( USE_FFTW )
  list( APPEND _test_srcs
    itkFFTW1DRealToComplexConjugateImageFilterTest.cxx
    itkFFTW1DComplexConjugateToRealImageFilterTest.cxx
    )
endif()
set(CMAKE_TESTDRIVER_BEFORE_TESTMAIN "#include \"itkTestDriverBeforeTest.inc\"")
set(CMAKE_TESTDRIVER_AFTER_TESTMAIN "#include \"itkTestDriverAfterTest.inc\"")
create_test_sourcelist( _test_driver_srcs ITKFFTExtensionsTestDriver.cxx
  ${_test_srcs}
  EXTRA_INCLUDE itkTestDriverIncludeRequiredIOFactories.h
  FUNCTION  ProcessArgumentsAndRegisterRequiredFactories
  )
add_executable( ITKFFTExtensionsTestDriver ITKFFTExtensionsTestDriver.cxx
  ${_test_driver_srcs} )
target_link_libraries( ITKFFTExtensionsTestDriver ${ITK_LIBRARIES} )

add_test( NAME itkAnalyticSignalImageFilterTest
  COMMAND ITKFFTExtensionsTestDriver
  --compare ${DATA_BASELINE}/itkAnalyticSignalImageFilterReal.mhd
    itkAnalyticSignalImageFilterTestOutputReal.mha
  --compare ${DATA_BASELINE}/itkAnalyticSignalImageFilterImaginary.mhd
    itkAnalyticSignalImageFilterTestOutputImaginary.mha
  itkAnalyticSignalImageFilterTest
    ${DATA_INPUT}/TreeBarkTexture.png
    itkAnalyticSignalImageFilterTestOutput
    )
add_test( NAME itkFFT1DComplexToComplexImageFilterTest
  COMMAND ITKFFTExtensionsTestDriver
  --compare ${DATA_INPUT}/TreeBarkTexture.png
            itkFFT1DComplexToComplexImageFilterTestOutput.mha
  itkFFT1DComplexToComplexImageFilterTest
    ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaseline
    itkFFT1DComplexToComplexImageFilterTestOutput.mha
    )
# use a forward fft, then inverse fft for this instead because of the FullMatrix
# issue
add_test( NAME itkFFT1DImageFilterTest
  COMMAND ITKFFTExtensionsTestDriver
  --compare ${DATA_INPUT}/TreeBarkTexture.png
            itkFFT1DImageFilterTestOutput.mha
  itkFFT1DImageFilterTest
    ${DATA_INPUT}/TreeBarkTexture.png
    itkFFT1DImageFilterTestOutput.mha
    )
add_test( NAME itkVnlFFT1DComplexToComplexImageFilterTest
  COMMAND ITKFFTExtensionsTestDriver
  --compare ${DATA_INPUT}/TreeBarkTexture.png
            itkVnlFFT1DComplexToComplexImageFilterTestOutput.mhd
  itkVnlFFT1DComplexToComplexImageFilterTest
    ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaseline
    itkVnlFFT1DComplexToComplexImageFilterTestOutput.mhd
    )
add_test( NAME itkVnlFFT1DRealToComplexConjugateImageFilterTest
  COMMAND ITKFFTExtensionsTestDriver
  --compare ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaselineRealFull.mhd
            itkVnlFFT1DRealToComplexConjugateImageFilterTestOutputReal.mha
  --compare ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaselineImaginaryFull.mhd
            itkVnlFFT1DRealToComplexConjugateImageFilterTestOutputImaginary.mha
  itkVnlFFT1DRealToComplexConjugateImageFilterTest
    ${DATA_INPUT}/TreeBarkTexture.png
    itkVnlFFT1DRealToComplexConjugateImageFilterTestOutput
    )
add_test( NAME itkVnlFFT1DComplexConjugateToRealImageFilterTest
  COMMAND ITKFFTExtensionsTestDriver
  --compare ${DATA_INPUT}/TreeBarkTexture.png
            itkVnlFFT1DComplexConjugateToRealImageFilterTestOutput.mhd
  itkVnlFFT1DComplexConjugateToRealImageFilterTest
    ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaseline
    itkVnlFFT1DComplexConjugateToRealImageFilterTestOutput.mhd
    )

if( USE_FFTW )
  add_test( NAME itkFFTW1DRealToComplexConjugateImageFilterTest
    COMMAND ITKFFTExtensionsTestDriver
    --compare ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaselineRealFull.mhd
              itkFFTW1DRealToComplexConjugateImageFilterTestOutputReal.mha
    --compare ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaselineImaginaryFull.mhd
              itkFFTW1DRealToComplexConjugateImageFilterTestOutputImaginary.mha
    itkFFTW1DRealToComplexConjugateImageFilterTest
      ${DATA_INPUT}/TreeBarkTexture.png
      itkFFTW1DRealToComplexConjugateImageFilterTestOutput
      )
  add_test( NAME itkFFTW1DComplexConjugateToRealImageFilterTest
    COMMAND ITKFFTExtensionsTestDriver
    --compare ${DATA_INPUT}/TreeBarkTexture.png
              itkFFTW1DComplexConjugateToRealImageFilterTestOutput.mha
    itkFFTW1DComplexConjugateToRealImageFilterTest
      ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaseline
      itkFFTW1DComplexConjugateToRealImageFilterTestOutput.mha
      )
  add_test( NAME itkFFTW1DComplexToComplexImageFilterTest
    COMMAND ITKFFTExtensionsTestDriver
    --compare ${DATA_INPUT}/TreeBarkTexture.png
              itkFFTW1DComplexToComplexImageFilterTestOutput.mha
    itkFFTW1DComplexToComplexImageFilterTest
      ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaseline
      itkFFTW1DComplexToComplexImageFilterTestOutput.mha
      )
endif()

#if( USE_OPENCL )
  #target_link_libraries( itkFFT1DImageFilterTest
    #OpenCL_FFT
    #)
  #target_link_libraries( itkFFT1DComplexToComplexImageFilterTest
    #OpenCL_FFT
    #)

  #add_executable( itkOpenCL1DRealToComplexConjugateImageFilterTest
    #itkOpenCL1DRealToComplexConjugateImageFilterTest.cxx
    #)
  #target_link_libraries( itkOpenCL1DRealToComplexConjugateImageFilterTest
    #ITKCommon
    #ITKIO
    #ITKAlgorithms
    #OpenCL_FFT
    #${OPENCL_LIBRARIES}
    #)
  ## gdb-pretty http://gitorious.com/gdb-pretty
  #add_custom_command (TARGET itkOpenCL1DRealToComplexConjugateImageFilterTest
    #POST_BUILD
    #COMMAND cmake -E copy
    #${CMAKE_CURRENT_SOURCE_DIR}/itkOpenCL1DRealToComplexConjugateImageFilterTest-gdb.py
    #${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
    #)
  #add_test( itkOpenCL1DRealToComplexConjugateImageFilterTest
    #${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkOpenCL1DRealToComplexConjugateImageFilterTest
    #--compare ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaselineReal128.mhd
	      #itkOpenCL1DRealToComplexConjugateImageFilterTestOutputReal.mha
    #--compare ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaselineImaginary128.mhd
	      #itkOpenCL1DRealToComplexConjugateImageFilterTestOutputImaginary.mha
    #itkOpenCL1DRealToComplexConjugateImageFilterTest
      #${DATA_INPUT}/TreeBarkTexture.png
      #itkOpenCL1DRealToComplexConjugateImageFilterTestOutput
      #)
  #add_executable( itkOpenCL1DComplexConjugateToRealImageFilterTest
    #itkOpenCL1DComplexConjugateToRealImageFilterTest.cxx
    #)
  #target_link_libraries( itkOpenCL1DComplexConjugateToRealImageFilterTest
    #ITKCommon
    #ITKIO
    #ITKAlgorithms
    #OpenCL_FFT
    #${OPENCL_LIBRARIES}
    #)
  #add_test( itkOpenCL1DComplexConjugateToRealImageFilterTest
    #${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkOpenCL1DComplexConjugateToRealImageFilterTest
    #--compare ${DATA_INPUT}/TreeBarkTexture.png
	      #itkOpenCL1DComplexConjugateToRealImageFilterTestOutput.mha
    #itkOpenCL1DComplexConjugateToRealImageFilterTest
      #${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaseline
      #itkOpenCL1DComplexConjugateToRealImageFilterTestOutput.mha
      #)
  #add_executable( itkOpenCL1DComplexToComplexImageFilterTest
    #itkOpenCL1DComplexToComplexImageFilterTest.cxx
    #)
  #target_link_libraries( itkOpenCL1DComplexToComplexImageFilterTest
    #ITKCommon
    #ITKIO
    #ITKAlgorithms
    #OpenCL_FFT
    #${OPENCL_LIBRARIES}
    #)
  #add_test( itkOpenCL1DComplexToComplexImageFilterTest
    #${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkOpenCL1DComplexToComplexImageFilterTest
    #--compare ${DATA_INPUT}/TreeBarkTexture.png
	      #itkOpenCL1DComplexToComplexImageFilterTestOutput.mha
    #itkOpenCL1DComplexToComplexImageFilterTest
      #${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaseline
      #itkOpenCL1DComplexToComplexImageFilterTestOutput.mha
      #)
#endif()

#if( USE_OPENCL )
  #target_link_libraries( itkAnalyticSignalImageFilterTest
    #OpenCL_FFT
    #${OPENCL_LIBRARIES}
    #)
#endif()
