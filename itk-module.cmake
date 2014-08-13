set(DOCUMENTATION "The modules provides FFT filters that operate only in one
direction of an image.  They may be particularly useful for ultrasound.")

itk_module(FFT1D
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
