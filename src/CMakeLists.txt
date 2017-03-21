set(Ultrasound_SRCS
  itkHDF5UltrasoundImageIOFactory.cxx
  itkHDF5UltrasoundImageIO.cxx
  )

itk_module_add_library(Ultrasound ${Ultrasound_SRCS})
