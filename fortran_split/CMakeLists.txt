set ( maooam_modules_files
    ${CMAKE_CURRENT_SOURCE_DIR}/coupled/maooam_coupled_wrapper.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/coupled/util.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/coupled/params.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/coupled/inprod_analytic.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/coupled/aotensor_def.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/coupled/ic_def.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/coupled/stat.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/coupled/tensor.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/coupled/tl_ad_tensor.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/coupled/rk2_integrator.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/coupled/rk2_tl_ad_integrator.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/coupled/ifport.f90
    PARENT_SCOPE
)
