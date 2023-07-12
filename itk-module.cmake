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
    ${_fft_depends}
  COMPILE_DEPENDS
    SplitComponents
    MeshToPolyData
    BSplineGradient
    HigherOrderAccurateGradient
    Strain
  PRIVATE_DEPENDS
    ITKHDF5
    ITKImageSources
    ITKEigen3
  TEST_DEPENDS
    ITKTestKernel
    ITKVideoCore
    ITKVideoIO
  FACTORY_NAMES
    "FFTImageFilterInit::Curvilinear"
    ImageIO::HDF5Ultrasound
  EXCLUDE_FROM_DEFAULT
  ENABLE_SHARED
  DESCRIPTION
    "${DOCUMENTATION}"
  )
