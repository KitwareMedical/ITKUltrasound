install( FILES
  itkFFT1DRealToComplexConjugateImageFilter.h
  itkFFT1DRealToComplexConjugateImageFilter.txx
  itkFFT1DComplexConjugateToRealImageFilter.h
  itkFFT1DComplexConjugateToRealImageFilter.txx
  itkVnl1DRealToComplexConjugateImageFilter.h
  itkVnl1DRealToComplexConjugateImageFilter.txx
  itkVnl1DComplexConjugateToRealImageFilter.h
  itkVnl1DComplexConjugateToRealImageFilter.txx
  DESTINATION ${ITK_INSTALL_PREFIX}/include/InsightToolkit/Algorithms
  )

if( USE_FFTW )
  install( FILES
  itkFFTW1DRealToComplexConjugateImageFilter.h
  itkFFTW1DRealToComplexConjugateImageFilter.txx
  DESTINATION ${ITK_INSTALL_PREFIX}/include/InsightToolkit/Algorithms
  )
endif()
