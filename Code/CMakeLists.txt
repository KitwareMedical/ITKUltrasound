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
    itkFFTW1DComplexConjugateToRealImageFilter.h
    itkFFTW1DComplexConjugateToRealImageFilter.txx
  DESTINATION ${ITK_INSTALL_PREFIX}/include/InsightToolkit/Algorithms
  )
endif()

if( USE_OPENCL )
  install( FILES
    itkOpenCL1DRealToComplexConjugateImageFilter.h
    itkOpenCL1DRealToComplexConjugateImageFilter.txx
    itkOpenCL1DComplexConjugateToRealImageFilter.h
    itkOpenCL1DComplexConjugateToRealImageFilter.txx
  DESTINATION ${ITK_INSTALL_PREFIX}/include/InsightToolkit/Algorithms
  )
endif()

