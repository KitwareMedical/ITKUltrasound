set(Ultrasound_SRCS
  itkHDF5UltrasoundImageIOFactory.cxx
  itkHDF5UltrasoundImageIO.cxx
  itkTextProgressBarCommand.cxx
  )

if(ITKUltrasound_USE_clFFT)
  list(APPEND Ultrasound_SRCS itkclFFTInitializer.cxx)
endif()

itk_module_add_library(Ultrasound ${Ultrasound_SRCS})
