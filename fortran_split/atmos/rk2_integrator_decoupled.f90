! integrator.f90
!
!>  Module with the integration routines.
!
!> @copyright                                                               
!> 2015 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!> Edited by Eviatar Bach (UMD) to allow forced atmos/ocean integration (2018)
!
!---------------------------------------------------------------------------!
!                                                                           
!>  @remark                                                                 
!>  This module actually contains the Heun algorithm routines.
!>  The user can modify it according to its preferred integration scheme.
!>  For higher-order schemes, additional buffers will probably have to be defined.
!                                                                           
!---------------------------------------------------------------------------

MODULE atmos_integrator
  USE atmos_params, only: ndim, natm, noc
  USE atmos_tensor, only:sparse_mul3
  USE atmos_aotensor_def, only: aotensor

  IMPLICIT NONE

  PUBLIC :: init_integrator, step

! PRIVATE

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_y1 !< Buffer to hold the intermediate position (Heun algorithm)
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_f0 !< Buffer to hold tendencies at the initial position
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: buf_f1 !< Buffer to hold tendencies at the intermediate position


CONTAINS
  
  !> Routine to initialise the integration buffers.
  SUBROUTINE init_integrator
    INTEGER :: AllocStat
    
    if (.not. allocated(buf_y1) .or. .not. allocated(buf_f0) .or. .not. allocated(buf_f1)) then !STEVE: fixing init clash
      ALLOCATE(buf_y1(0:ndim),buf_f0(0:ndim),buf_f1(0:ndim) ,STAT=AllocStat)
    endif
    IF (AllocStat /= 0) STOP "*** Not enough memory ! ***"
  END SUBROUTINE init_integrator
  
  !> Routine computing the tendencies of the model
  !> @param t Time at which the tendencies have to be computed. Actually not needed for autonomous systems.
  !> @param y Point at which the tendencies have to be computed.
  !> @param res vector to store the result.
  !> @remark Note that it is NOT safe to pass `y` as a result buffer, 
  !> as this operation does multiple passes.
  SUBROUTINE tendencies(t,y,res)
    REAL(KIND=8), INTENT(IN) :: t
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res
    CALL sparse_mul3(aotensor, y, y, res)
  END SUBROUTINE tendencies

  !> Routine to perform an integration step (Heun algorithm). The incremented time is returned.
  !> @param y Initial point.
  !> @param t Actual integration time
  !> @param dt Integration timestep.
  !> @param res Final point after the step.
  SUBROUTINE step(y,t,dt,res,component)
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: y
    REAL(KIND=8), INTENT(INOUT) :: t
    REAL(KIND=8), INTENT(IN) :: dt
    REAL(KIND=8), DIMENSION(0:ndim), INTENT(OUT) :: res
    CHARACTER*3, OPTIONAL, INTENT(IN) :: component

    CALL tendencies(t,y,buf_f0)

    if (present(component)) then
      if (component == 'atm') then
        ! Zero out ocean tendencies
        buf_f0(1+2*natm:ndim) = 0.0d0
      endif
    endif

    buf_y1 = y+dt*buf_f0

    CALL tendencies(t+dt,buf_y1,buf_f1)

    if (present(component)) then
      if (component == 'atm') then
        ! Zero out ocean tendencies
        buf_f1(1+2*natm:ndim) = 0.0d0
      endif
    endif

    res=y+0.5*(buf_f0+buf_f1)*dt

    t=t+dt

  END SUBROUTINE step

END MODULE atmos_integrator
