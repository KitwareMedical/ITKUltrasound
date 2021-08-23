set(DOCUMENTATION "The modules provides filters may be particularly useful for
ultrasound image reconstruction and analysis.")

set(_fft_depends ITKFFT)
if(ITKUltrasound_USE_clFFT)
  list(APPEND _fft_depends ITKGPUCommon)
endif()

itk_module(Ultrasound
  DEPENDS
    ITKAnisotropicSmoothing
    ITKImageCompose
    ITKImageFunction
    ITKImageIntensity
    ITKIOImageBase
    ITKTransform
    ITKRegistrationCommon
    ITKVideoCore
    ITKVideoIO
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
