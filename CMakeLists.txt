cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)
# If all the submodules are installed and 'git submodule init && git submodule
# update --recursive' is not called.
cmake_policy( SET CMP0014 OLD )


project( ITK_FFT_EXTENSIONS )

include_directories( ${ITK_FFT_EXTENSIONS_SOURCE_DIR}/Code )

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/CMake )

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )
if( NOT ITK_USE_REVIEW AND NOT ITKReview_LOADED )
  message( SEND_ERROR "ITK must be configured with Module_ITKReview set to ON." )
endif()


if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb3 -fno-inline")

  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -DNDEBUG -ffast-math")

  set(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} -g")
  set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} -O")
endif(CMAKE_COMPILER_IS_GNUCXX)

# where all the executables go
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY
  ${ITK_FFT_EXTENSIONS_BINARY_DIR}/bin )

if( ITK_USE_FFTWF OR ITK_USE_FFTWD )
  set( USE_FFTW ON )
else()
  set( USE_FFTW OFF )
endif()

add_subdirectory( Code )


include(CTest)
if(BUILD_TESTING)
  enable_testing()
  add_subdirectory( Testing )
endif(BUILD_TESTING)
