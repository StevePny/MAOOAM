set ( maooam_modules_files
    ${CMAKE_CURRENT_SOURCE_DIR}/util.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/params.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/inprod_analytic.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/aotensor_def.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ic_def.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/stat.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/tensor.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/tl_ad_tensor.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/rk2_integrator.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/rk2_tl_ad_integrator.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ifport.f90
    PARENT_SCOPE
)
