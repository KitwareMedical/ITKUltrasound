set(Ultrasound_SRCS
  itkHDF5UltrasoundImageIOFactory.cxx
  itkHDF5UltrasoundImageIO.cxx
  )

add_library(Ultrasound ${ITK_LIBRARY_BUILD_TYPE} ${Ultrasound_SRCS})
itk_module_link_dependencies()
itk_module_target(Ultrasound)
