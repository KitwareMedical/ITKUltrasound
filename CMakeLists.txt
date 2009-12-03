cmake_minimum_required(VERSION 2.8)
if(COMMAND cmake_policy)
  cmake_policy(VERSION 2.8)
endif(COMMAND cmake_policy)


project( ITK_FFT_EXTENSIONS )


find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )


if(CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb3 -fno-inline")

  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -DNDEBUG -ffast-math")

  set(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} -g")
  set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} -O")
endif(CMAKE_COMPILER_IS_GNUCXX)

# where all the executables go
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY
  ${STREAMING_RESAMPLE_IMAGE_BINARY_DIR}/bin )


add_subdirectory( Code )


include(CTest)
if(BUILD_TESTING)
  find_path( ITK_DATA_ROOT ITKData.readme
    $ENV{ITK_DATA_ROOT}
    DOC "Path corresponding to to ${ITK_SOURCE_DIR}/Testing/Data."
    )
  if( NOT ITK_DATA_ROOT )
    message( SEND_ERROR "ITK_DATA_ROOT must be defined to perform tests." )
  endif()
  enable_testing()
  add_subdirectory( Testing )
endif(BUILD_TESTING)
