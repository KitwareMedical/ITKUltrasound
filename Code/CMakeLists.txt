install( FILES
  itkFFT1DRealToComplexConjugateImageFilter.h
  itkFFT1DRealToComplexConjugateImageFilter.hxx
  itkFFT1DComplexConjugateToRealImageFilter.h
  itkFFT1DComplexConjugateToRealImageFilter.hxx
  itkFFT1DComplexToComplexImageFilter.h
  itkFFT1DComplexToComplexImageFilter.hxx
  itkVnlFFT1DRealToComplexConjugateImageFilter.h
  itkVnlFFT1DRealToComplexConjugateImageFilter.hxx
  itkVnlFFT1DComplexConjugateToRealImageFilter.h
  itkVnlFFT1DComplexConjugateToRealImageFilter.hxx
  itkVnlFFT1DComplexToComplexImageFilter.h
  itkVnlFFT1DComplexToComplexImageFilter.hxx
  itkAnalyticSignalImageFilter.h
  itkAnalyticSignalImageFilter.hxx
  DESTINATION include/InsightToolkit/Algorithms
  )

if( USE_FFTW )
  install( FILES
    itkFFTW1DRealToComplexConjugateImageFilter.h
    itkFFTW1DRealToComplexConjugateImageFilter.hxx
    itkFFTW1DComplexConjugateToRealImageFilter.h
    itkFFTW1DComplexConjugateToRealImageFilter.hxx
    itkFFTW1DComplexToComplexImageFilter.h
    itkFFTW1DComplexToComplexImageFilter.hxx
    itkFFTWCommonExtended.h
  DESTINATION include/InsightToolkit/Algorithms
  )
endif()
