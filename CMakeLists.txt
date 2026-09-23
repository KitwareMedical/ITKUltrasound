cmake_minimum_required(VERSION 3.16.3)
project(Ultrasound)

if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()

set(Ultrasound_LIBRARIES Ultrasound)

if(NOT ITK_SOURCE_DIR)
  find_package(ITK REQUIRED)
  list(APPEND CMAKE_MODULE_PATH ${ITK_CMAKE_DIR})
endif()

option(ITKUltrasound_USE_clFFT "Build ITKUltrasound 1D FFT filters with clFFT as the back-end.\n\
  https://github.com/clMathLibraries/clFFT" ${ITK_USE_GPU})
if(ITKUltrasound_USE_clFFT)
  find_package(clFFT REQUIRED)
  list(APPEND Ultrasound_LIBRARIES ${CLFFT_LIBRARIES})

  list(APPEND Ultrasound_SYSTEM_INCLUDE_DIRS ${CLFFT_INCLUDE_DIRS})
  get_target_property(Ultrasound_clFFT_LIBRARY_DIRS clFFT LINK_DIRECTORIES)  
  set(Ultrasound_SYSTEM_LIBRARY_DIRS ${Ultrasound_clFFT_LIBRARY_DIRS})

  set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS ITKUltrasound_USE_clFFT)
  set(Ultrasound_EXPORT_CODE_INSTALL "
  set_property(DIRECTORY APPEND PROPERTY INCLUDE_DIRECTORIES ${CLFFT_INCLUDE_DIRS})
  set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS ITKUltrasound_USE_clFFT)
  ")
  set(Ultrasound_EXPORT_CODE_BUILD "
  if(NOT ITK_BINARY_DIR)
    set_property(DIRECTORY APPEND PROPERTY INCLUDE_DIRECTORIES ${CLFFT_INCLUDE_DIRS})
    set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS ITKUltrasound_USE_clFFT)
  endif()
  ")
endif()

# Optional VTK support. Keep version in sync with itk-module-init.cmake.
set(ITKUltrasound_VTK_MINIMUM_VERSION 9.1)
find_package(VTK ${ITKUltrasound_VTK_MINIMUM_VERSION} NO_MODULE QUIET)
set(itkultrasound_use_vtk_default OFF)
if(VTK_FOUND)
  set(itkultrasound_use_vtk_default ON)
endif()
option(ITKUltrasound_USE_VTK "Build ITKUltrasound filters that require VTK." ${itkultrasound_use_vtk_default})
mark_as_advanced(ITKUltrasound_USE_VTK)

if(ITKUltrasound_USE_VTK)
  find_package(VTK ${ITKUltrasound_VTK_MINIMUM_VERSION} NO_MODULE REQUIRED
    COMPONENTS
      CommonCore
      CommonDataModel
      IOLegacy # vtkStructuredGridWriter, used by the tests
    )

  # itkSpecialCoordinatesImageToVTKStructuredGridFilter is header-only, so
  # consumers of this module need VTK's usage requirements.
  list(APPEND Ultrasound_LIBRARIES VTK::CommonCore VTK::CommonDataModel)

  # Make the VTK targets available to projects that use this module.
  # Append, so that the clFFT-related export code above is preserved.
  string(APPEND Ultrasound_EXPORT_CODE_INSTALL "
  if(NOT TARGET VTK::CommonDataModel)
    set(VTK_DIR \"${VTK_DIR}\")
    find_package(VTK ${ITKUltrasound_VTK_MINIMUM_VERSION} NO_MODULE REQUIRED COMPONENTS CommonCore CommonDataModel)
  endif()
  ")
  string(APPEND Ultrasound_EXPORT_CODE_BUILD "
  if(NOT ITK_BINARY_DIR AND NOT TARGET VTK::CommonDataModel)
    set(VTK_DIR \"${VTK_DIR}\")
    find_package(VTK ${ITKUltrasound_VTK_MINIMUM_VERSION} NO_MODULE REQUIRED COMPONENTS CommonCore CommonDataModel)
  endif()
  ")
endif() # ITKUltrasound_USE_VTK

if(NOT ITK_SOURCE_DIR)
  include(ITKModuleExternal)
else()
  itk_module_impl()
endif()
