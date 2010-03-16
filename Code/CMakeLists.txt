install( FILES
  itkFFT1DRealToComplexConjugateImageFilter.h
  itkFFT1DRealToComplexConjugateImageFilter.txx
  itkFFT1DComplexConjugateToRealImageFilter.h
  itkFFT1DComplexConjugateToRealImageFilter.txx
  itkFFT1DComplexToComplexImageFilter.h
  itkFFT1DComplexToComplexImageFilter.txx
  itkVnlFFT1DRealToComplexConjugateImageFilter.h
  itkVnlFFT1DRealToComplexConjugateImageFilter.txx
  itkVnlFFT1DComplexConjugateToRealImageFilter.h
  itkVnlFFT1DComplexConjugateToRealImageFilter.txx
  itkVnlFFT1DComplexToComplexImageFilter.h
  itkVnlFFT1DComplexToComplexImageFilter.txx
  itkAnalyticSignalImageFilter.h
  itkAnalyticSignalImageFilter.txx
  DESTINATION include/InsightToolkit/Algorithms
  )

if( USE_FFTW )
  install( FILES
    itkFFTW1DRealToComplexConjugateImageFilter.h
    itkFFTW1DRealToComplexConjugateImageFilter.txx
    itkFFTW1DComplexConjugateToRealImageFilter.h
    itkFFTW1DComplexConjugateToRealImageFilter.txx
    itkFFTW1DComplexToComplexImageFilter.h
    itkFFTW1DComplexToComplexImageFilter.txx
    itkFFTWCommonExtended.h
  DESTINATION include/InsightToolkit/Algorithms
  )
endif()

if( USE_OPENCL )
  install( FILES
    itkOpenCL1DRealToComplexConjugateImageFilter.h
    itkOpenCL1DRealToComplexConjugateImageFilter.txx
    itkOpenCL1DComplexConjugateToRealImageFilter.h
    itkOpenCL1DComplexConjugateToRealImageFilter.txx
    itkOpenCL1DComplexToComplexImageFilter.h
    itkOpenCL1DComplexToComplexImageFilter.txx
  DESTINATION include/InsightToolkit/Algorithms
  )
endif()

