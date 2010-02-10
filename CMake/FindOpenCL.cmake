# - Find OpenCL cross-platform parallel programming library
#
# This module defines the following non-cached variables:
#  OPENCL_FOUND        - TRUE if OpenCL was found
#  OPENCL_INCLUDE_DIRS - Include directories of OpenCL (not cached)
#  OPENCL_LIBRARIES    - Libraries to link against (not cached)
#
# The following cached variables are also defined, but are not intended for
# general use:
#  OPENCL_INCLUDE_DIR  - Directory containing CL/cl.h or OpenCL/cl.h (cached)
#  OPENCL_LIBRARY      - The OpenCL library (cached)
#

# Copyright (C) 2009 Michael Wild <themiwi@users.sf.net>

find_path(OPENCL_INCLUDE_DIR NAMES CL/cl.h OpenCL/cl.h
  PATHS ENV OPENCL_DIR
  PATH_SUFFIXES include
  )

# if not an Apple framework, provide a search hint for the
# library based on the include path and a path suffix based on whether the
# system is 32 or 64 bit.
if(NOT OPENCL_INCLUDE_DIR MATCHES "\\.framework/")
  get_filename_component(__opencl_lib_dir_hint "${OPENCL_INCLUDE_DIR}" PATH)
  set(__opencl_lib_dir_hint HINTS "${__opencl_lib_dir_hint}/lib")
  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(__opencl_libdir_suffix PATH_SUFFIXES x86_64)
  else()
    set(__opencl_libdir_suffix PATH_SUFFIXES x86)
  endif()
else()
  set(__opencl_lib_dir_hint)
  set(__opencl_libdir_suffix)
endif()

find_library(OPENCL_LIBRARY NAMES OpenCL
  ${__opencl_lib_dir_hint}
  ${__opencl_libdir_suffix}
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenCL DEFAULT_MSG
  OPENCL_LIBRARY OPENCL_INCLUDE_DIR)

set(OPENCL_INCLUDE_DIRS ${OPENCL_INCLUDE_DIR})
set(OPENCL_LIBRARIES ${OPENCL_LIBRARY})

mark_as_advanced(OPENCL_INCLUDE_DIR OPENCL_LIBRARY)
