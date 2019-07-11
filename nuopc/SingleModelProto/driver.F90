!==============================================================================
! Earth System Modeling Framework
! Copyright 2002-2019, University Corporation for Atmospheric Research, 
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics 
! Laboratory, University of Michigan, National Centers for Environmental 
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory, 
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

module driver

  !-----------------------------------------------------------------------------
  ! Code that specializes generic NUOPC_Driver
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Driver, &
    driver_routine_SS             => SetServices, &
    driver_label_SetModelServices => label_SetModelServices
  
  use MODEL, only: modelSS => SetServices
  
  implicit none
  
  private
  
  ! private module data --> ONLY PARAMETERS
  integer, parameter            :: stepCount = 1
  real(ESMF_KIND_R8), parameter :: stepTime  = 3600.D0  ! step time [s]
                                                      ! should be parent step

  public SetServices
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

  subroutine SetServices(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc
    logical :: local_verbose = .true.
    
    rc = ESMF_SUCCESS
    
    ! NUOPC_Driver registers the generic methods
    if (local_verbose) print *, "driver::SetServices:: calling NUOPC_CompDerive..."
    call NUOPC_CompDerive(driver, driver_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (local_verbose) print *, "driver::SetServices:: finished NUOPC_CompDerive."
      
    ! attach specializing method(s)
    if (local_verbose) print *, "driver::SetServices:: calling NUOPC_CompSpecialize..."
    call NUOPC_CompSpecialize(driver, specLabel=driver_label_SetModelServices, &
      specRoutine=SetModelServices, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (local_verbose) print *, "driver::SetServices:: finished NUOPC_CompSpecialize."

    ! set driver verbosity
!   if (local_verbose) print *, "driver::SetServices:: calling NUOPC_CompAttributeSet..."
!   call NUOPC_CompAttributeSet(driver, name="Verbosity", value="low", rc=rc)
!   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!     line=__LINE__, &
!     file=__FILE__)) &
!     return  ! bail out
!   if (local_verbose) print *, "driver::SetServices:: finished NUOPC_CompAttributeSet."

  end subroutine SetServices

  !-----------------------------------------------------------------------------

  subroutine SetModelServices(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_GridComp)           :: child
    type(ESMF_CplComp)            :: connector
    type(ESMF_Time)               :: startTime
    type(ESMF_Time)               :: stopTime
    type(ESMF_TimeInterval)       :: timeStep
    type(ESMF_Clock)              :: internalClock

    logical :: local_verbose = .true.

    rc = ESMF_SUCCESS
    
    ! SetServices for MODEL component
    if (local_verbose) print *, "driver::SetModelServices:: calling NUOPC_DriverAddComp..."
    call NUOPC_DriverAddComp(driver, "MODEL", modelSS, comp=child, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (local_verbose) print *, "driver::SetModelServices:: finished NUOPC_DriverAddComp."

    if (local_verbose) print *, "driver::SetModelServices:: calling NUOPC_CompAttributeSet..."
!   call NUOPC_CompAttributeSet(child, name="Verbosity", value="low", rc=rc)
!   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!     line=__LINE__, &
!     file=__FILE__)) &
!     return  ! bail out
    if (local_verbose) print *, "driver::SetModelServices:: finished NUOPC_CompAttributeSet."
      
    ! set the driver clock
    if (local_verbose) print *, "driver::SetModelServices:: calling ESMF_TimeSet..."
    call ESMF_TimeSet(startTime, s = 0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    if (local_verbose) print *, "driver::SetModelServices:: finished ESMF_TimeSet."

    if (local_verbose) print *, "driver::SetModelServices:: calling ESMF_TimeSet..."
    call ESMF_TimeSet(stopTime, s_r8 = stepTime * stepCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    if (local_verbose) print *, "driver::SetModelServices:: finished ESMF_TimeSet."

    if (local_verbose) print *, "driver::SetModelServices:: calling ESMF_TimeIntervalSet..."
    call ESMF_TimeIntervalSet(timeStep, s_r8 = stepTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    if (local_verbose) print *, "driver::SetModelServices:: finished ESMF_TimeIntervalSet."

    if (local_verbose) print *, "driver::SetModelServices:: calling ESMF_ClockCreate..."
    internalClock = ESMF_ClockCreate(name="Driver Clock", &
      timeStep=timeStep, startTime=startTime, stopTime=stopTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
    if (local_verbose) print *, "driver::SetModelServices:: finished ESMF_ClockCreate."
      
    if (local_verbose) print *, "driver::SetModelServices:: calling ESMF_GridCompSet..."
    call ESMF_GridCompSet(driver, clock=internalClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (local_verbose) print *, "driver::SetModelServices:: finished ESMF_GridCompSet."
      
  end subroutine SetModelServices

end module
