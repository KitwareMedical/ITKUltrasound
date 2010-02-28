set( DATA_INPUT ${ITK_FFT_EXTENSIONS_SOURCE_DIR}/Testing/Data/Input )
set( DATA_BASELINE ${ITK_FFT_EXTENSIONS_SOURCE_DIR}/Testing/Data/Baseline )

# use a forward fft, then inverse fft for this instead because of the FullMatrix
# issue
add_executable( itkFFT1DImageFilterTest
  itkFFT1DImageFilterTest.cxx
  )
target_link_libraries( itkFFT1DImageFilterTest
  ITKCommon
  ITKIO
  ITKAlgorithms
  )
add_test( itkFFT1DImageFilterTest
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkFFT1DImageFilterTest
  --compare ${ITK_DATA_ROOT}/Input/TreeBarkTexture.png
	    itkFFT1DImageFilterTestOutput.mha
  itkFFT1DImageFilterTest
    ${ITK_DATA_ROOT}/Input/TreeBarkTexture.png
    itkFFT1DImageFilterTestOutput.mha
    )
add_executable( itkFFT1DComplexToComplexImageFilterTest
  itkFFT1DComplexToComplexImageFilterTest.cxx
  )
target_link_libraries( itkFFT1DComplexToComplexImageFilterTest
  ITKCommon
  ITKIO
  ITKAlgorithms
  )
add_test( itkFFT1DComplexToComplexImageFilterTest
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkFFT1DComplexToComplexImageFilterTest
  --compare ${ITK_DATA_ROOT}/Input/TreeBarkTexture.png
	    itkFFT1DComplexToComplexImageFilterTestOutput.mha
  itkFFT1DComplexToComplexImageFilterTest
    ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaseline
    itkFFT1DComplexToComplexImageFilterTestOutput.mha
    )

add_executable( itkVnlFFT1DRealToComplexConjugateImageFilterTest
  itkVnlFFT1DRealToComplexConjugateImageFilterTest.cxx
  )
target_link_libraries( itkVnlFFT1DRealToComplexConjugateImageFilterTest
  ITKCommon
  ITKIO
  )
add_test( itkVnlFFT1DRealToComplexConjugateImageFilterTest
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkVnlFFT1DRealToComplexConjugateImageFilterTest
  --compare ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaselineRealFull.mhd
	    itkVnlFFT1DRealToComplexConjugateImageFilterTestOutputReal.mha
  --compare ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaselineImaginaryFull.mhd
	    itkVnlFFT1DRealToComplexConjugateImageFilterTestOutputImaginary.mha
  itkVnlFFT1DRealToComplexConjugateImageFilterTest
    ${ITK_DATA_ROOT}/Input/TreeBarkTexture.png
    itkVnlFFT1DRealToComplexConjugateImageFilterTestOutput
    )
add_executable( itkVnlFFT1DComplexConjugateToRealImageFilterTest
  itkVnlFFT1DComplexConjugateToRealImageFilterTest.cxx
  )
target_link_libraries( itkVnlFFT1DComplexConjugateToRealImageFilterTest
  ITKCommon
  ITKIO
  )
add_test( itkVnlFFT1DComplexConjugateToRealImageFilterTest
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkVnlFFT1DComplexConjugateToRealImageFilterTest
  --compare ${ITK_DATA_ROOT}/Input/TreeBarkTexture.png
	    itkVnlFFT1DComplexConjugateToRealImageFilterTestOutput.mhd
  itkVnlFFT1DComplexConjugateToRealImageFilterTest
    ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaseline
    itkVnlFFT1DComplexConjugateToRealImageFilterTestOutput.mhd
    )
add_executable( itkVnlFFT1DComplexToComplexImageFilterTest
  itkVnlFFT1DComplexToComplexImageFilterTest.cxx
  )
target_link_libraries( itkVnlFFT1DComplexToComplexImageFilterTest
  ITKCommon
  ITKIO
  )
add_test( itkVnlFFT1DComplexToComplexImageFilterTest
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkVnlFFT1DComplexToComplexImageFilterTest
  --compare ${ITK_DATA_ROOT}/Input/TreeBarkTexture.png
	    itkVnlFFT1DComplexToComplexImageFilterTestOutput.mhd
  itkVnlFFT1DComplexToComplexImageFilterTest
    ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaseline
    itkVnlFFT1DComplexToComplexImageFilterTestOutput.mhd
    )

if( USE_FFTW )
  add_executable( itkFFTW1DRealToComplexConjugateImageFilterTest
    itkFFTW1DRealToComplexConjugateImageFilterTest.cxx
    )
  target_link_libraries( itkFFTW1DRealToComplexConjugateImageFilterTest
    ITKCommon
    ITKIO
    ITKAlgorithms
    )
  add_test( itkFFTW1DRealToComplexConjugateImageFilterTest
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkFFTW1DRealToComplexConjugateImageFilterTest
    --compare ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaselineRealFull.mhd
	      itkFFTW1DRealToComplexConjugateImageFilterTestOutputReal.mha
    --compare ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaselineImaginaryFull.mhd
	      itkFFTW1DRealToComplexConjugateImageFilterTestOutputImaginary.mha
    itkFFTW1DRealToComplexConjugateImageFilterTest
      ${ITK_DATA_ROOT}/Input/TreeBarkTexture.png
      itkFFTW1DRealToComplexConjugateImageFilterTestOutput
      )
  add_executable( itkFFTW1DComplexConjugateToRealImageFilterTest
    itkFFTW1DComplexConjugateToRealImageFilterTest.cxx
    )
  target_link_libraries( itkFFTW1DComplexConjugateToRealImageFilterTest
    ITKCommon
    ITKIO
    ITKAlgorithms
    )
  add_test( itkFFTW1DComplexConjugateToRealImageFilterTest
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkFFTW1DComplexConjugateToRealImageFilterTest
    --compare ${ITK_DATA_ROOT}/Input/TreeBarkTexture.png
	      itkFFTW1DComplexConjugateToRealImageFilterTestOutput.mha
    itkFFTW1DComplexConjugateToRealImageFilterTest
      ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaseline
      itkFFTW1DComplexConjugateToRealImageFilterTestOutput.mha
      )
  add_executable( itkFFTW1DComplexToComplexImageFilterTest
    itkFFTW1DComplexToComplexImageFilterTest.cxx
    )
  target_link_libraries( itkFFTW1DComplexToComplexImageFilterTest
    ITKCommon
    ITKIO
    ITKAlgorithms
    )
  add_test( itkFFTW1DComplexToComplexImageFilterTest
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkFFTW1DComplexToComplexImageFilterTest
    --compare ${ITK_DATA_ROOT}/Input/TreeBarkTexture.png
	      itkFFTW1DComplexToComplexImageFilterTestOutput.mha
    itkFFTW1DComplexToComplexImageFilterTest
      ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaseline
      itkFFTW1DComplexToComplexImageFilterTestOutput.mha
      )
endif()

if( USE_OPENCL )
  target_link_libraries( itkFFT1DImageFilterTest
    OpenCL_FFT
    )
  target_link_libraries( itkFFT1DComplexToComplexImageFilterTest
    OpenCL_FFT
    )

  add_executable( itkOpenCL1DRealToComplexConjugateImageFilterTest
    itkOpenCL1DRealToComplexConjugateImageFilterTest.cxx
    )
  target_link_libraries( itkOpenCL1DRealToComplexConjugateImageFilterTest
    ITKCommon
    ITKIO
    ITKAlgorithms
    OpenCL_FFT
    ${OPENCL_LIBRARIES}
    )
  # gdb-pretty http://gitorious.com/gdb-pretty
  add_custom_command (TARGET itkOpenCL1DRealToComplexConjugateImageFilterTest
    POST_BUILD
    COMMAND cmake -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/itkOpenCL1DRealToComplexConjugateImageFilterTest-gdb.py
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
    )
  add_test( itkOpenCL1DRealToComplexConjugateImageFilterTest
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkOpenCL1DRealToComplexConjugateImageFilterTest
    --compare ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaselineReal128.mhd
	      itkOpenCL1DRealToComplexConjugateImageFilterTestOutputReal.mha
    --compare ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaselineImaginary128.mhd
	      itkOpenCL1DRealToComplexConjugateImageFilterTestOutputImaginary.mha
    itkOpenCL1DRealToComplexConjugateImageFilterTest
      ${ITK_DATA_ROOT}/Input/TreeBarkTexture.png
      itkOpenCL1DRealToComplexConjugateImageFilterTestOutput
      )
  add_executable( itkOpenCL1DComplexConjugateToRealImageFilterTest
    itkOpenCL1DComplexConjugateToRealImageFilterTest.cxx
    )
  target_link_libraries( itkOpenCL1DComplexConjugateToRealImageFilterTest
    ITKCommon
    ITKIO
    ITKAlgorithms
    OpenCL_FFT
    ${OPENCL_LIBRARIES}
    )
  add_test( itkOpenCL1DComplexConjugateToRealImageFilterTest
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkOpenCL1DComplexConjugateToRealImageFilterTest
    --compare ${ITK_DATA_ROOT}/Input/TreeBarkTexture.png
	      itkOpenCL1DComplexConjugateToRealImageFilterTestOutput.mha
    itkOpenCL1DComplexConjugateToRealImageFilterTest
      ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaseline
      itkOpenCL1DComplexConjugateToRealImageFilterTestOutput.mha
      )
  add_executable( itkOpenCL1DComplexToComplexImageFilterTest
    itkOpenCL1DComplexToComplexImageFilterTest.cxx
    )
  target_link_libraries( itkOpenCL1DComplexToComplexImageFilterTest
    ITKCommon
    ITKIO
    ITKAlgorithms
    OpenCL_FFT
    ${OPENCL_LIBRARIES}
    )
  add_test( itkOpenCL1DComplexToComplexImageFilterTest
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkOpenCL1DComplexToComplexImageFilterTest
    --compare ${ITK_DATA_ROOT}/Input/TreeBarkTexture.png
	      itkOpenCL1DComplexToComplexImageFilterTestOutput.mha
    itkOpenCL1DComplexToComplexImageFilterTest
      ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaseline
      itkOpenCL1DComplexToComplexImageFilterTestOutput.mha
      )
endif()

add_executable( itkAnalyticSignalImageFilterTest
  itkAnalyticSignalImageFilterTest.cxx
  )
target_link_libraries( itkAnalyticSignalImageFilterTest
  ITKCommon
  ITKIO
  ITKAlgorithms
  )
if( USE_OPENCL )
  target_link_libraries( itkAnalyticSignalImageFilterTest
    OpenCL_FFT
    ${OPENCL_LIBRARIES}
    )
endif()
add_test( itkAnalyticSignalImageFilterTest
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkAnalyticSignalImageFilterTest
  --compare ${DATA_BASELINE}/itkAnalyticSignalImageFilterReal.mhd
    itkAnalyticSignalImageFilterTestOutputReal.mha
  --compare ${DATA_BASELINE}/itkAnalyticSignalImageFilterImaginary.mhd
    itkAnalyticSignalImageFilterTestOutputImaginary.mha
  itkAnalyticSignalImageFilterTest
    ${ITK_DATA_ROOT}/Input/TreeBarkTexture.png
    itkAnalyticSignalImageFilterTestOutput
    )
# gdb-pretty http://gitorious.com/gdb-pretty
add_custom_command (TARGET itkAnalyticSignalImageFilterTest
  POST_BUILD
  COMMAND cmake -E copy
  ${CMAKE_CURRENT_SOURCE_DIR}/itkAnalyticSignalImageFilterTest-gdb.py
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
  )
