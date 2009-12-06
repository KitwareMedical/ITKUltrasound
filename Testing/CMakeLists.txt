include_directories( ${ITK_FFT_EXTENSIONS_SOURCE_DIR}/Code )

add_executable( itkFFT1DRealToComplexConjugateImageFilterTest
  itkFFT1DRealToComplexConjugateImageFilterTest.cxx
  )
target_link_libraries( itkFFT1DRealToComplexConjugateImageFilterTest
  ITKCommon
  ITKIO
  )
add_test( itkFFT1DRealToComplexConjugateImageFilterTest
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkFFT1DRealToComplexConjugateImageFilterTest
  --compare itkFFT1DRealToComplexConjugateImageFilterTestBaselineReal.mhd
	    itkFFT1DRealToComplexConjugateImageFilterTestOutputReal.mha
	    itkFFT1DRealToComplexConjugateImageFilterTestBaselineImaginary.mhd
	    itkFFT1DRealToComplexConjugateImageFilterTestOutputImaginary.mha
  itkFFT1DRealToComplexConjugateImageFilterTest
    ${ITK_DATA_ROOT}/Input/TreeBarkTexture.png
    itkFFT1DRealToComplexConjugateImageFilterTestOutput
    )
