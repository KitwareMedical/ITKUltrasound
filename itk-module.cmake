set(DOCUMENTATION "The modules provides filters may be particularly useful for
ultrasound image reconstruction and analysis.")

# A library is only created when FFTW support is ON
set(_fft_depends ITKFFT)
if(ITK_USE_FFTWF OR ITK_USE_FFTWD)
  set(_fft_depends ITKFFT)
endif()

itk_module(Ultrasound
  DEPENDS
    ITKImageCompose
    ITKImageIntensity
    ITKIOImageBase
    ITKTransform
    ${_fft_depends}
  PRIVATE_DEPENDS
    ITKHDF5
  TEST_DEPENDS
    ITKTransform
    ITKTestKernel
    ITKImageSources
    ${_fft_depends}
  EXCLUDE_FROM_DEFAULT
  ENABLE_SHARED
  DESCRIPTION
    "${DOCUMENTATION}"
  )
