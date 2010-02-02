install( FILES
  itkVnl1DRealToComplexConjugateImageFilter.h
  itkVnl1DRealToComplexConjugateImageFilter.txx
  DESTINATION ${ITK_INSTALL_PREFIX}/include/InsightToolkit/Algorithms
  )

if( USE_FFTW )
  install( FILES
  itkFFT1DRealToComplexConjugateImageFilter.h
  itkFFT1DRealToComplexConjugateImageFilter.txx
  DESTINATION ${ITK_INSTALL_PREFIX}/include/InsightToolkit/Algorithms
  )
endif()
