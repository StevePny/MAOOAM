# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

################################################################################
# Modular Arbitrary-Order Ocean-Atmosphere Model (MAOOAM)
# https://github.com/StevePny/MAOOAM
# branch forked from
# https://github.com/Climdyn/MAOOAM
################################################################################

cmake_minimum_required( VERSION 3.3.2 FATAL_ERROR )

project( maooam Fortran )

set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH} )

set( ECBUILD_DEFAULT_BUILD_TYPE Release )
set( ENABLE_OS_TESTS           OFF CACHE BOOL "Disable OS tests" FORCE )
set( ENABLE_LARGE_FILE_SUPPORT OFF CACHE BOOL "Disable testing of large file support" FORCE )

include( ecbuild_system NO_POLICY_SCOPE )

ecbuild_requires_macro_version( 2.7 )

################################################################################
# Project
################################################################################

ecbuild_declare_project()

ecbuild_enable_fortran( REQUIRED )

################################################################################
# Dependencies
################################################################################


################################################################################
# Definitions
################################################################################


################################################################################
# Export package info
################################################################################

list( APPEND MAOOAM_TPLS )

set( MAOOAM_INCLUDE_DIRS ${CMAKE_Fortran_MODULE_DIRECTORY})
set( MAOOAM_LIBRARIES maooam )

get_directory_property( MAOOAM_DEFINITIONS COMPILE_DEFINITIONS )

foreach( _tpl ${MAOOAM_TPLS} )
  string( TOUPPER ${_tpl} TPL )
  list( APPEND MAOOAM_EXTRA_DEFINITIONS   ${${TPL}_DEFINITIONS}  ${${TPL}_TPL_DEFINITIONS}  )
  list( APPEND MAOOAM_EXTRA_INCLUDE_DIRS  ${${TPL}_INCLUDE_DIRS} ${${TPL}_TPL_INCLUDE_DIRS} )
  list( APPEND MAOOAM_EXTRA_LIBRARIES     ${${TPL}_LIBRARIES}    ${${TPL}_TPL_LIBRARIES}    )
endforeach()

################################################################################
# Sources
################################################################################

set( MAOOAM_LINKER_LANGUAGE Fortran )

include( maooam_compiler_flags )
include_directories( ${MAOOAM_INCLUDE_DIRS} )

add_subdirectory(fortran_split)

list ( APPEND maooam_src_files
	${maooam_modules_files}
	${maooam_toolbox_files}
)

################################################################################
# Library
################################################################################

ecbuild_add_library( TARGET   maooam
                     SOURCES  ${maooam_src_files}
                     LINKER_LANGUAGE ${MAOOAM_LINKER_LANGUAGE}
                    )

################################################################################
# Finalise configuration
################################################################################

# prepares a tar.gz of the sources and/or binaries
ecbuild_install_project( NAME maooam )

# print the summary of the configuration
ecbuild_print_summary()
