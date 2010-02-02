include_directories( ${ITK_FFT_EXTENSIONS_SOURCE_DIR}/Code )

set( DATA_INPUT ${ITK_FFT_EXTENSIONS_SOURCE_DIR}/Testing/Data/Input )

# use a forward fft, then inverse fft for this instead because of the FullMatrix
# issue
#add_executable( itkFFT1DRealToComplexConjugateImageFilterTest
  #itkFFT1DRealToComplexConjugateImageFilterTest.cxx
  #)
#target_link_libraries( itkFFT1DRealToComplexConjugateImageFilterTest
  #ITKCommon
  #ITKIO
  #ITKAlgorithms
  #)
#add_test( itkFFT1DRealToComplexConjugateImageFilterTest
  #${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkFFT1DRealToComplexConjugateImageFilterTest
  #--compare ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaselineReal.mhd
	    #itkFFT1DRealToComplexConjugateImageFilterTestOutputReal.mha
  #--compare ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaselineImaginary.mhd
	    #itkFFT1DRealToComplexConjugateImageFilterTestOutputImaginary.mha
  #itkFFT1DRealToComplexConjugateImageFilterTest
    #${ITK_DATA_ROOT}/Input/TreeBarkTexture.png
    #itkFFT1DRealToComplexConjugateImageFilterTestOutput
    #)

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
endif()
