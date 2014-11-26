set(DOCUMENTATION "The modules provides filters may be particularly useful for
ultrasound image reconstruction and analysis.")

itk_module(Ultrasound
  DEPENDS
    ITKFFT
    ITKImageCompose
    ITKImageIntensity
  TEST_DEPENDS
    ITKTestKernel
  EXCLUDE_FROM_DEFAULT
  DESCRIPTION
    "${DOCUMENTATION}"
  )
