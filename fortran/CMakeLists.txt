#set ( maooam_modules_files
#    ${CMAKE_CURRENT_SOURCE_DIR}/util.f90
#    ${CMAKE_CURRENT_SOURCE_DIR}/params.f90
#    ${CMAKE_CURRENT_SOURCE_DIR}/inprod_analytic.f90
#    ${CMAKE_CURRENT_SOURCE_DIR}/aotensor_def.f90
#    ${CMAKE_CURRENT_SOURCE_DIR}/ic_def.f90
#    ${CMAKE_CURRENT_SOURCE_DIR}/stat.f90
#    ${CMAKE_CURRENT_SOURCE_DIR}/tensor.f90
#    ${CMAKE_CURRENT_SOURCE_DIR}/tl_ad_tensor.f90
#    ${CMAKE_CURRENT_SOURCE_DIR}/rk2_integrator.f90
#    ${CMAKE_CURRENT_SOURCE_DIR}/rk2_tl_ad_integrator.f90
#    ${CMAKE_CURRENT_SOURCE_DIR}/ifport.f90
#    PARENT_SCOPE
#)

cmake_minimum_required(VERSION 3.12)

project(maooam Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-O2 -Wall -cpp")
    set(debug "-g -O0 -fbounds-check -Wall -Wextra -Wconversion -pedantic -ffpe-trap=zero,overflow,underflow -fsanitize=address -cpp")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-assume byterecl -O2 -fpp")
    set(debug "-O0 -check all -traceback -fpe0 -check bounds -debug all -check uninit -ftrapuv -assume byterecl -fpp")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${debug}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")
find_package(LAPACK REQUIRED)

# Copy config files
file(GLOB nml_files  ./*.nml)
file(COPY ${nml_files} DESTINATION .)

#
# Compile.
#
file(GLOB sources  ./*.f90)
list(FILTER sources EXCLUDE REGEX test)
add_executable(maooam ${sources})
target_link_libraries(maooam ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

