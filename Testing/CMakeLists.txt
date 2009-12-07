include_directories( ${ITK_FFT_EXTENSIONS_SOURCE_DIR}/Code )

set( DATA_INPUT ${ITK_FFT_EXTENSIONS_SOURCE_DIR}/Testing/Data/Input )


add_executable( itkFFT1DRealToComplexConjugateImageFilterTest
  itkFFT1DRealToComplexConjugateImageFilterTest.cxx
  )
target_link_libraries( itkFFT1DRealToComplexConjugateImageFilterTest
  ITKCommon
  ITKIO
  )
add_test( itkFFT1DRealToComplexConjugateImageFilterTest
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/itkFFT1DRealToComplexConjugateImageFilterTest
  --compare ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaselineReal.mhd
	    itkFFT1DRealToComplexConjugateImageFilterTestOutputReal.mha
  --compare ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaselineImaginary.mhd
	    itkFFT1DRealToComplexConjugateImageFilterTestOutputImaginary.mha
  itkFFT1DRealToComplexConjugateImageFilterTest
    ${ITK_DATA_ROOT}/Input/TreeBarkTexture.png
    itkFFT1DRealToComplexConjugateImageFilterTestOutput
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
  --compare ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaselineReal.mhd
	    itkVnlFFT1DRealToComplexConjugateImageFilterTestOutputReal.mha
  --compare ${DATA_INPUT}/itkFFT1DRealToComplexConjugateImageFilterTestBaselineImaginary.mhd
	    itkVnlFFT1DRealToComplexConjugateImageFilterTestOutputImaginary.mha
  itkVnlFFT1DRealToComplexConjugateImageFilterTest
    ${ITK_DATA_ROOT}/Input/TreeBarkTexture.png
    itkVnlFFT1DRealToComplexConjugateImageFilterTestOutput
    )
