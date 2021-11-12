if(${ITK_WRAP_double} AND NOT ${ITK_WRAP_complex_double})
  message(WARNING "ITK_WRAP_double is on but ITK_WRAP_complex_double is not. Some filters with double input / outputs depend on complex filters. Set ITK_WRAP_complex_double to on" )
endif()

if(${ITK_WRAP_float} AND NOT ${ITK_WRAP_complex_float})
  message(WARNING "ITK_WRAP_float is on but ITK_WRAP_complex_float is not. Some filters with float input / outputs depend on complex filters. Set ITK_WRAP_complex_float to on" )
endif()

itk_wrap_module(Ultrasound)
set(WRAPPER_LIBRARY_GROUPS
  itkBlockMatchingMetricImageFilter
  itkBlockMatchingNormalizedCrossCorrelationMetricImageFilter
  itkBlockMatchingNormalizedCrossCorrelationFFTMetricImageFilter
  itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter
  itkSpectra1DSupportWindowImageFilter
  itkCurvilinearArraySpecialCoordinatesImage
  itkFrequencyDomain1DFilterFunction
  itkFrequencyDomain1DImageFilter
  itkForward1DFFTImageFilter
  itkInverse1DFFTImageFilter
  itkRegionFromReferenceImageFilter
  itkAnalyticSignalImageFilter
  itkBModeImageFilter
  itkSpeckleReducingAnisotropicDiffusionImageFilter
  )
itk_auto_load_submodules()
itk_end_wrap_module()
