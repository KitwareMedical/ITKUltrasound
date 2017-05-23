itk_wrap_module(Ultrasound)

set(WRAPPER_LIBRARY_GROUPS
  itkSpectra1DSupportWindowImageFilter
  itkCurvilinearArraySpecialCoordinatesImage
  itkFrequencyDomain1DFilterFunction
  itkFrequencyDomain1DImageFilter
  )
itk_auto_load_submodules()

itk_end_wrap_module()
