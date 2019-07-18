!==============================================================================
! Earth System Modeling Framework
! Copyright 2002-2019, University Corporation for Atmospheric Research, 
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics 
! Laboratory, University of Michigan, National Centers for Environmental 
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory, 
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

! #define TESTAUTOADDCONNECTORS
#define FFLOG

module ESM

  !-----------------------------------------------------------------------------
  ! Code that specializes generic ESM Component code.
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Driver, inheritDriver => SetServices
  
  use ATM, only: atmSS => SetServices
  use OCN, only: ocnSS => SetServices
  use MED, only: medSS => SetServices
  
  use NUOPC_Connector, only: cplSS => SetServices
  
  implicit none
  
  private
  
  public SetServices
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

  subroutine SetServices(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc
    
    type(ESMF_Config)           :: config

    logical :: local_verbose = .true.

    rc = ESMF_SUCCESS
    
    ! NUOPC_Driver registers the generic methods
    if (local_verbose) print *, "ESM::SetServices:: Calling NUOPC_CompDerive..."
    call NUOPC_CompDerive(driver, inheritDriver, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! attach specializing method(s)
    if (local_verbose) print *, "ESM::SetServices:: Calling NUOPC_CompSpecialize..."
    call NUOPC_CompSpecialize(driver, specLabel=label_SetModelServices, specRoutine=SetModelServices, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    if (local_verbose) print *, "ESM::SetServices:: Calling NUOPC_CompSpecialize..."
    call NUOPC_CompSpecialize(driver, specLabel=label_SetRunSequence, specRoutine=SetRunSequence, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! create, open and set the config
    if (local_verbose) print *, "ESM::SetServices:: Calling ESMF_ConfigCreate..."
    config = ESMF_ConfigCreate(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (local_verbose) print *, "ESM::SetServices:: Calling ESMF_ConfigLoadFile..."
    call ESMF_ConfigLoadFile(config, "esmApp.runconfig", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (local_verbose) print *, "ESM::SetServices:: Calling ESMF_GridCompSet..."
    call ESMF_GridCompSet(driver, config=config, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine SetServices

  !-----------------------------------------------------------------------------

  subroutine SetModelServices(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_GridComp)           :: child
    type(ESMF_Time)               :: startTime
    type(ESMF_Time)               :: stopTime
    type(ESMF_TimeInterval)       :: timeStep
    type(ESMF_Clock)              :: internalClock
    type(ESMF_Config)             :: config
    type(NUOPC_FreeFormat)        :: attrFF

    logical :: local_verbose = .true.
    
    rc = ESMF_SUCCESS
    
    ! read free format driver attributes
    if (local_verbose) print *, "ESM::SetModelServices:: Calling ESMF_GridCompGet..."
    call ESMF_GridCompGet(driver, config=config, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (local_verbose) print *, "ESM::SetModelServices:: Calling NUOPC_FreeFormatCreate..."
    attrFF = NUOPC_FreeFormatCreate(config, label="driverAttributes::", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! ingest FreeFormat driver attributes
    if (local_verbose) print *, "ESM::SetModelServices:: Calling NUOPC_CompAttributeIngest..."
    call NUOPC_CompAttributeIngest(driver, attrFF, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! clean-up
    if (local_verbose) print *, "ESM::SetModelServices:: Calling NUOPC_FreeFormatDestroy..."
    call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
  
    ! SetServices for ATM
    if (local_verbose) print *, "ESM::SetModelServices:: Calling NUOPC_DriverAddComp..."
    call NUOPC_DriverAddComp(driver, "ATM", atmSS, comp=child, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! set default ATM attributes
    if (local_verbose) print *, "ESM::SetModelServices:: Calling NUOPC_CompAttributeSet..."
    call NUOPC_CompAttributeSet(child, name="Verbosity", value="1", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! read ATM attributes from config file into FreeFormat
    if (local_verbose) print *, "ESM::SetModelServices:: Calling NUOPC_FreeFormatCreate..."
    attrFF = NUOPC_FreeFormatCreate(config, label="atmAttributes::", relaxedflag=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! ingest FreeFormat atm attributes
    if (local_verbose) print *, "ESM::SetModelServices:: Calling NUOPC_CompAttributeIngest..."
    call NUOPC_CompAttributeIngest(child, attrFF, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! clean-up
    if (local_verbose) print *, "ESM::SetModelServices:: Calling NUOPC_FreeFormatDestroy..."
    call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! SetServices for OCN
    if (local_verbose) print *, "ESM::SetModelServices:: Calling NUOPC_DriverAddComp..."
    call NUOPC_DriverAddComp(driver, "OCN", ocnSS, comp=child, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! set default OCN attributes
    if (local_verbose) print *, "ESM::SetModelServices:: Calling NUOPC_CompAttributeSet..."
    call NUOPC_CompAttributeSet(child, name="Verbosity", value="1", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! read OCN attributes from config file into FreeFormat
    if (local_verbose) print *, "ESM::SetModelServices:: Calling NUOPC_FreeFormatCreate..."
    attrFF = NUOPC_FreeFormatCreate(config, label="ocnAttributes::", relaxedflag=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! ingest FreeFormat ocn attributes
    if (local_verbose) print *, "ESM::SetModelServices:: Calling NUOPC_CompAttributeIngest..."
    call NUOPC_CompAttributeIngest(child, attrFF, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! clean-up
    if (local_verbose) print *, "ESM::SetModelServices:: Calling NUOPC_FreeFormatDestroy..."
    call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! SetServices for MED
    if (local_verbose) print *, "ESM::SetModelServices:: Calling NUOPC_DriverAddComp..."
    call NUOPC_DriverAddComp(driver, "MED", medSS, comp=child, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! set default MED attributes
    if (local_verbose) print *, "ESM::SetModelServices:: Calling NUOPC_CompAttributeSet..."
    call NUOPC_CompAttributeSet(child, name="Verbosity", value="1", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! read MED attributes from config file into FreeFormat
    if (local_verbose) print *, "ESM::SetModelServices:: Calling NUOPC_FreeFormatCreate..."
    attrFF = NUOPC_FreeFormatCreate(config, label="medAttributes::", relaxedflag=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! ingest FreeFormat med attributes
    if (local_verbose) print *, "ESM::SetModelServices:: Calling NUOPC_CompAttributeIngest..."
    call NUOPC_CompAttributeIngest(child, attrFF, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    ! clean-up
    if (local_verbose) print *, "ESM::SetModelServices:: Calling NUOPC_FreeFormatDestroy..."
    call NUOPC_FreeFormatDestroy(attrFF, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

#ifndef TESTAUTOADDCONNECTORS
    ! SetServices for atm2med
    if (local_verbose) print *, "ESM::SetModelServices:: Calling NUOPC_DriverAddComp..."
    call NUOPC_DriverAddComp(driver, srcCompLabel="ATM", dstCompLabel="MED", &
      compSetServicesRoutine=cplSS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! SetServices for ocn2med
    if (local_verbose) print *, "ESM::SetModelServices:: Calling NUOPC_DriverAddComp..."
    call NUOPC_DriverAddComp(driver, srcCompLabel="OCN", dstCompLabel="MED", &
      compSetServicesRoutine=cplSS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! SetServices for med2atm
    if (local_verbose) print *, "ESM::SetModelServices:: Calling NUOPC_DriverAddComp..."
    call NUOPC_DriverAddComp(driver, srcCompLabel="MED", dstCompLabel="ATM", &
      compSetServicesRoutine=cplSS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! SetServices for med2ocn
    if (local_verbose) print *, "ESM::SetModelServices:: Calling NUOPC_DriverAddComp..."
    call NUOPC_DriverAddComp(driver, srcCompLabel="MED", dstCompLabel="OCN", &
      compSetServicesRoutine=cplSS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
#endif

    ! set the driver clock
    if (local_verbose) print *, "ESM::SetModelServices:: Calling ESMF_TimeIntervalSet..."
    call ESMF_TimeIntervalSet(timeStep, s=10, rc=rc) ! 10 second default step
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (local_verbose) print *, "ESM::SetModelServices:: Calling ESMF_TimeIntervalSet..."
    call ESMF_TimeSet(startTime, yy=2010, mm=6, dd=1, h=0, m=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (local_verbose) print *, "ESM::SetModelServices:: Calling ESMF_TimeSet..."
    call ESMF_TimeSet(stopTime, yy=2010, mm=6, dd=1, h=0, m=1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (local_verbose) print *, "ESM::SetModelServices:: Calling ESMF_ClockCreate..."
    internalClock = ESMF_ClockCreate(name="Application Clock", &
      timeStep=timeStep, startTime=startTime, stopTime=stopTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    if (local_verbose) print *, "ESM::SetModelServices:: Calling ESMF_GridCompSet..."
    call ESMF_GridCompSet(driver, clock=internalClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine SetModelServices

  !-----------------------------------------------------------------------------

  subroutine SetRunSequence(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc
    
    ! local variables
    character(ESMF_MAXSTR)        :: name, connName
    type(ESMF_Config)             :: config
    type(NUOPC_FreeFormat)        :: ff
    type(ESMF_CplComp), pointer   :: connectorList(:)
    integer                       :: i

    logical :: local_verbose = .true.

    rc = ESMF_SUCCESS
    
    ! query the Component for info
    if (local_verbose) print *, "ESM::SetRunSequence:: calling ESMF_GridCompGet..."
    call ESMF_GridCompGet(driver, name=name, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    
    ! read free format run sequence from config
    if (local_verbose) print *, "ESM::SetRunSequence:: calling ESMF_GridCompGet..."
    call ESMF_GridCompGet(driver, config=config, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    ff = NUOPC_FreeFormatCreate(config, label="runSeq::", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
      
#ifdef FFLOG
    if (local_verbose) print *, "ESM::SetRunSequence:: calling NUOPC_FreeFormatLog..."
    call NUOPC_FreeFormatLog(ff, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#endif

    ! ingest FreeFormat run sequence
    if (local_verbose) print *, "ESM::SetRunSequence:: calling NUOPC_DriverIngestRunSequence..."
    call NUOPC_DriverIngestRunSequence(driver, ff, &
#ifdef TESTAUTOADDCONNECTORS
      autoAddConnectors=.true., &
#endif
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

#ifdef FFLOG
    ! Diagnostic output
    if (local_verbose) print *, "ESM::SetRunSequence:: calling NUOPC_DriverPrint..."
    call NUOPC_DriverPrint(driver, orderflag=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#endif

    ! clean-up
    if (local_verbose) print *, "ESM::SetRunSequence:: calling NUOPC_FreeFormatDestroy..."
    call NUOPC_FreeFormatDestroy(ff, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
      
    ! set Verbosity on Connectors no matter how they were added
    ! also read in potential connector attributes from config
    nullify(connectorList)
    if (local_verbose) print *, "ESM::SetRunSequence:: calling NUOPC_DriverGetComp..."
    call NUOPC_DriverGetComp(driver, connectorList, rc=rc)    
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    do i=1, size(connectorList)
      if (local_verbose) print *, "ESM::SetRunSequence:: connectorList iteration i = ", i, " / ", size(connectorList)

      ! default Verbosity, may be overridden from config below
      if (local_verbose) print *, "ESM::SetRunSequence:: calling NUOPC_CompAttributeSet..."
      call NUOPC_CompAttributeSet(comp=connectorList(i), name="Verbosity", value="4097", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      ! read connector Attributes from config
!STEVE:TEST
      if (local_verbose) print *, "ESM::SetRunSequence:: calling NUOPC_CompGet..."
      call NUOPC_CompGet(comp=connectorList(i), name=connName, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

!     if (local_verbose) print *, "ESM::SetRunSequence:: calling NUOPC_CompAttributeGet..."
!     call NUOPC_CompAttributeGet(comp=connectorList(i), name=connName, rc=rc)
!     if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!       line=__LINE__, &
!       file=__FILE__)) &
!       return  ! bail out

      if (local_verbose) print *, "ESM::SetRunSequence:: calling ESMF_LogWrite..."
      call ESMF_LogWrite("Reading Attributes for Connector: "//trim(connName), &
        ESMF_LOGMSG_INFO, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      if (local_verbose) print *, "ESM::SetRunSequence:: calling NUOPC_FreeFormatCreate..."
      ff = NUOPC_FreeFormatCreate(config, relaxedflag=.true., &
        label=trim(connName)//"-Attributes::", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

#ifdef FFLOG 
      if (local_verbose) print *, "ESM::SetRunSequence:: calling NUOPC_FreeFormatLog..."
      call NUOPC_FreeFormatLog(ff, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
#endif

      ! ingest FreeFormat driver attributes
      if (local_verbose) print *, "ESM::SetRunSequence:: calling NUOPC_CompAttributeIngest..."
      call NUOPC_CompAttributeIngest(connectorList(i), ff, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      ! clean-up
      if (local_verbose) print *, "ESM::SetRunSequence:: calling NUOPC_FreeFormatDestroy..."
      call NUOPC_FreeFormatDestroy(ff, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

    enddo
    
  end subroutine SetRunSequence

  !-----------------------------------------------------------------------------

end module ESM
