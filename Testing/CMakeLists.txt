include_directories( ${ITK_FFT_EXTENSIONS_SOURCE_DIR}/Code )

set( DATA_INPUT ${ITK_FFT_EXTENSIONS_SOURCE_DIR}/Testing/Data/Input )

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
    --compare ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaselineRealNotFull.mhd
	      itkFFTW1DRealToComplexConjugateImageFilterTestOutputReal.mha
    --compare ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaselineImaginaryNotFull.mhd
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
endif()
