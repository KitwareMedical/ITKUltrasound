if(${ITK_WRAP_double} AND NOT ${ITK_WRAP_complex_double})
  message(WARNING "ITK_WRAP_double is on but ITK_WRAP_complex_double is not. Some filters with double input / outputs depend on complex filters. Set ITK_WRAP_complex_double to on" )
endif()

if(${ITK_WRAP_float} AND NOT ${ITK_WRAP_complex_float})
  message(WARNING "ITK_WRAP_float is on but ITK_WRAP_complex_float is not. Some filters with float input / outputs depend on complex filters. Set ITK_WRAP_complex_float to on" )
endif()

itk_wrap_module(Ultrasound)
set(WRAPPER_SUBMODULE_ORDER
  itkSpecialCoordinatesImage
  itkCurvilinearArraySpecialCoordinatesImage
  itkCurvilinearArraySpecialCoordinatesImageFilters
  itkSliceSeriesSpecialCoordinatesImage
  itkSliceSeriesSpecialCoordinatesImageFilters
  itkAttenuationImageFilter
  itkBlockMatchingMetricImageFilter
  itkBlockMatchingNormalizedCrossCorrelationMetricImageFilter
  itkBlockMatchingNormalizedCrossCorrelationFFTMetricImageFilter
  itkBlockMatchingNormalizedCrossCorrelationNeighborhoodIteratorMetricImageFilter
  itkImageUltrasound # must come before Spectra1DSupportWindow
  itkSpectra1DSupportWindowImageFilter
  itkSpectra1DSupportWindowToMaskImageFilter  # must come after CASCI SSSCI
  itkFrequencyDomain1DFilterFunction
  itkFrequencyDomain1DImageFilter
  itkRegionFromReferenceImageFilter
  itkAnalyticSignalImageFilter
  itkBModeImageFilter
  itkSpeckleReducingAnisotropicDiffusionImageFilter
  itkImageToImageFilterUltrasound # Extra wrappings for Spectra1D filters
  itkSpectra1DAveragingImageFilter
  itkSpectra1DNormalizeImageFilter
)

itk_auto_load_submodules()
itk_end_wrap_module()
