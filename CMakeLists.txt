project(Ultrasound)
cmake_minimum_required(VERSION 3.16.3)

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

# required VTK version
find_package(VTK NO_MODULE QUIET)
set(itkultrasound_use_vtk_default OFF)
if(VTK_FOUND)
  set(itkultrasound_use_vtk_default ON)
endif()
option(ITKUltrasound_USE_VTK "Build ITKUltrasound filters that require VTK." ${itkultrasound_use_vtk_default})
mark_as_advanced(ITKUltrasound_USE_VTK)

if(ITKUltrasound_USE_VTK)
  find_package(VTK NO_MODULE REQUIRED)

  # Add compiler flags needed to use VTK.
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${VTK_REQUIRED_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${VTK_REQUIRED_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${VTK_REQUIRED_EXE_LINKER_FLAGS}")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${VTK_REQUIRED_SHARED_LINKER_FLAGS}")
  set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${VTK_REQUIRED_MODULE_LINKER_FLAGS}")

  vtk_module_config(Ultrasound_VTK
    vtkCommonDataModel
    vtkIOLegacy
    )

  list(APPEND Ultrasound_LIBRARIES ${VTK_LIBRARIES})

  # The VTK DICOMParser and vtkmetaio includes conflict with the ITK
  # versions. Here we remove them from the include directories.
  #
  string(REGEX REPLACE "[^;]*MetaIO;"
           "" Ultrasound_VTK_INCLUDE_DIRS "${Ultrasound_VTK_INCLUDE_DIRS}")
  string(REGEX REPLACE "[^;]*vtkmetaio;"
           "" Ultrasound_VTK_INCLUDE_DIRS "${Ultrasound_VTK_INCLUDE_DIRS}")
  string(REGEX REPLACE "[^;]*DICOMParser;"
           "" Ultrasound_VTK_INCLUDE_DIRS "${Ultrasound_VTK_INCLUDE_DIRS}")

  # Set up the include and library dirs needed for the module
  # Use append syntax to avoid overwriting possible clFFT-related content
  list(APPEND Ultrasound_SYSTEM_INCLUDE_DIRS ${Ultrasound_VTK_INCLUDE_DIRS})
  list(APPEND Ultrasound_SYSTEM_LIBRARY_DIRS ${Ultrasound_VTK_LIBRARY_DIRS})

  # Add preprocessor definitions needed to use VTK.
  set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS ${Ultrasound_VTK_DEFINITIONS})

  set(Ultrasound_EXPORT_CODE_INSTALL "
  set(VTK_DIR \"${VTK_DIR}\")

  find_package(VTK NO_MODULE REQUIRED)

  vtk_module_config(Ultrasound_VTK
    vtkCommonDataModel
    )

  set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS \${Ultrasound_VTK_DEFINITIONS})
  ")
  set(Ultrasound_EXPORT_CODE_BUILD "
  if(NOT ITK_BINARY_DIR)
    set(VTK_DIR \"${VTK_DIR}\")

    find_package(VTK NO_MODULE REQUIRED)

    vtk_module_config(Ultrasound_VTK
      vtkCommonDataModel
      )

    set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS \${Ultrasound_VTK_DEFINITIONS})
  endif()
  ")
endif() # ITKUltrasound_USE_VTK

if(NOT ITK_SOURCE_DIR)
  include(ITKModuleExternal)
else()
  itk_module_impl()
endif()
