install( FILES
  itkVnl1DRealToComplexConjugateImageFilter.h
  itkVnl1DRealToComplexConjugateImageFilter.txx
  DESTINATION ${ITK_INSTALL_PREFIX}/include/InsightToolkit/Algorithms
  )

if( USE_FFTWD OR USE_FFTWF )
  install( FILES
  itkFFT1DRealToComplexConjugateImageFilter.h
  itkFFT1DRealToComplexConjugateImageFilter.txx
  DESTINATION ${ITK_INSTALL_PREFIX}/include/InsightToolkit/Algorithms
  )
endif()
