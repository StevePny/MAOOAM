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

    ${CMAKE_CURRENT_SOURCE_DIR}/atmos/maooam_atmos_wrapper.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/atmos/util.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/atmos/params.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/atmos/inprod_analytic.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/atmos/aotensor_def.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/atmos/ic_def.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/atmos/stat.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/atmos/tensor.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/atmos/tl_ad_tensor.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/atmos/rk2_integrator_decoupled.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/atmos/rk2_tl_ad_integrator_decoupled.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/atmos/ifport.f90

    ${CMAKE_CURRENT_SOURCE_DIR}/ocean/maooam_ocean_wrapper.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ocean/util.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ocean/params.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ocean/inprod_analytic.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ocean/aotensor_def.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ocean/ic_def.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ocean/stat.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ocean/tensor.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ocean/tl_ad_tensor.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ocean/rk2_integrator_decoupled.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ocean/rk2_tl_ad_integrator_decoupled.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/ocean/ifport.f90

    PARENT_SCOPE
)
