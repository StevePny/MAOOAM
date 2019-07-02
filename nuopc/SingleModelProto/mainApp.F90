!==============================================================================
! Earth System Modeling Framework
! Copyright 2002-2019, University Corporation for Atmospheric Research, 
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics 
! Laboratory, University of Michigan, National Centers for Environmental 
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory, 
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

program mainApp

  !-----------------------------------------------------------------------------
  ! Generic ESMF Main
  !-----------------------------------------------------------------------------

  use ESMF

  use driver, only: &
    driver_SS => SetServices

  implicit none
  
  integer                       :: rc, userRc
  type(ESMF_GridComp)           :: drvComp

  logical :: local_verbose = .true.

  ! Initialize ESMF
  if (local_verbose) print *, "mainApp:: calling ESMF_Initialize..."
  call ESMF_Initialize(defaultCalkind=ESMF_CALKIND_GREGORIAN, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (local_verbose) print *, "mainApp:: finished ESMF_Initialize."
    
  if (local_verbose) print *, "mainApp:: calling ESMF_LogWrite..."
  call ESMF_LogWrite("mainApp STARTING", ESMF_LOGMSG_INFO, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (local_verbose) print *, "mainApp:: finished ESMF_LogWrite."
    
  !-----------------------------------------------------------------------------
  
  ! -> CREATE THE DRIVER
  if (local_verbose) print *, "mainApp:: calling ESMF_GridCompCreate..."
  drvComp = ESMF_GridCompCreate(name="driver", rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (local_verbose) print *, "mainApp:: finished ESMF_GridCompCreate."
    
  ! -> SET DRIVER SERVICES
  if (local_verbose) print *, "mainApp:: calling ESMF_GridCompSetServices..."
  call ESMF_GridCompSetServices(drvComp, driver_SS, userRc=userRc, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (ESMF_LogFoundError(rcToCheck=userRc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (local_verbose) print *, "mainApp:: finished ESMF_GridCompSetServices."

  ! INITIALIZE THE DRIVER
  if (local_verbose) print *, "mainApp:: calling ESMF_GridCompInitialize..."
  call ESMF_GridCompInitialize(drvComp, userRc=userRc, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (ESMF_LogFoundError(rcToCheck=userRc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (local_verbose) print *, "mainApp:: finished ESMF_GridCompInitialize."
      
  ! RUN THE DRIVER
  if (local_verbose) print *, "mainApp:: calling ESMF_GridCompRun..."
  call ESMF_GridCompRun(drvComp, userRc=userRc, rc=rc)
  if (local_verbose) print *, "mainApp:: calling rc check..."
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (local_verbose) print *, "mainApp:: calling userRc check..."
  if (ESMF_LogFoundError(rcToCheck=userRc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (local_verbose) print *, "mainApp:: finished ESMF_GridCompRun."
  
  ! FINALIZE THE DRIVER
  if (local_verbose) print *, "mainApp:: calling ESMF_GridCompFinalize..."
  call ESMF_GridCompFinalize(drvComp, userRc=userRc, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (ESMF_LogFoundError(rcToCheck=userRc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (local_verbose) print *, "mainApp:: finished ESMF_GridCompFinalize."

  !-----------------------------------------------------------------------------
  
  if (local_verbose) print *, "mainApp:: calling ESMF_LogWrite..."
  call ESMF_LogWrite("mainApp FINISHED", ESMF_LOGMSG_INFO, rc=rc)
  if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
    line=__LINE__, &
    file=__FILE__)) &
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if (local_verbose) print *, "mainApp:: finished ESMF_LogWrite."

  ! Finalize ESMF
  if (local_verbose) print *, "mainApp:: calling ESMF_Finalize..."
  call ESMF_Finalize()
  if (local_verbose) print *, "mainApp:: finished ESMF_Finalize."

  print *, "End Program."
  
end program mainApp 
