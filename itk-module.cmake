set(DOCUMENTATION "The modules provides filters may be particularly useful for
ultrasound image reconstruction and analysis.")

# A library is only created when FFTW support is ON
set(_fft_depends COMPILE_DEPENDS ITKFFT)
if(ITK_USE_FFTWF OR ITK_USE_FFTWD)
  set(_fft_depends ITKFFT)
endif()

itk_module(Ultrasound
  DEPENDS
    ITKImageCompose
    ITKImageIntensity
    ITKIOImageBase
    ${_fft_depends}
  TEST_DEPENDS
    ITKTestKernel
    ITKImageSources
  EXCLUDE_FROM_DEFAULT
  DESCRIPTION
    "${DOCUMENTATION}"
  )
