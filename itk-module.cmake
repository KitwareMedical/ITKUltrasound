set(DOCUMENTATION "The modules provides filters may be particularly useful for
ultrasound image reconstruction and analysis.")

set(_fft_depends ITKFFT)
if(ITK_USE_GPU)
  list(APPEND _fft_depends ITKGPUCommon) # also add OpenCL_FFT, and replace by clFFT
endif()

itk_module(Ultrasound
  DEPENDS
    ITKImageCompose
    ITKImageIntensity
    ITKIOImageBase
    ITKTransform
    ITKRegistrationCommon
    ${_fft_depends}
  COMPILE_DEPENDS
    SplitComponents
    BSplineGradient
    HigherOrderAccurateGradient
    Strain
  PRIVATE_DEPENDS
    ITKHDF5
    ITKTestKernel
    ITKImageSources
  EXCLUDE_FROM_DEFAULT
  ENABLE_SHARED
  DESCRIPTION
    "${DOCUMENTATION}"
  )
