!==============================================================================
! Earth System Modeling Framework
! Copyright 2002-2019, University Corporation for Atmospheric Research, 
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics 
! Laboratory, University of Michigan, National Centers for Environmental 
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory, 
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

module ESM

  !-----------------------------------------------------------------------------
  ! Code that specializes generic ESM Component code.
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Driver, &
    driver_routine_SS             => SetServices, &
    driver_label_SetModelServices => label_SetModelServices, &
    driver_label_ModifyCplLists   => label_ModifyCplLists
  
  use ATM, only: atmSS => SetServices
  use OCN, only: ocnSS => SetServices
  
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

    logical :: local_verbose = .true.
    
    rc = ESMF_SUCCESS
    
    ! NUOPC_Driver registers the generic methods
    if (local_verbose) print *, "ESM::SetServices calling NUOPC_CompDerive..."
    call NUOPC_CompDerive(driver, driver_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! attach specializing method(s)
    if (local_verbose) print *, "ESM::SetServices calling NUOPC_CompSpecialize..."
    call NUOPC_CompSpecialize(driver, specLabel=driver_label_SetModelServices, &
      specRoutine=SetModelServices, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(driver, specLabel=driver_label_ModifyCplLists, &
      specRoutine=ModifyCplLists, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! set driver verbosity
    call NUOPC_CompAttributeSet(driver, name="Verbosity", value="1", rc=rc)
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
!   type(ESMF_Grid)               :: grid
!   type(ESMF_Field)              :: field
    type(ESMF_Time)               :: startTime
    type(ESMF_Time)               :: stopTime
    type(ESMF_TimeInterval)       :: timeStep
    type(ESMF_Clock)              :: internalClock
    type(ESMF_GridComp)           :: child
    type(ESMF_CplComp)            :: connector

    logical :: local_verbose = .true.
    character(12) :: method
    character(ESMF_MAXSTR) :: test_cplList
    character(ESMF_MAXSTR), pointer :: ptr_cplList
    integer :: test_verbosity

    rc = ESMF_SUCCESS
    
    ! SetServices for ATM
    if (local_verbose) print *, "ESM::SetModelServices:: calling NUOPC_DriverAddComp for ATM..."
    call NUOPC_DriverAddComp(driver, "ATM", atmSS, comp=child, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompAttributeSet(child, name="Verbosity", value="1", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! SetServices for OCN
    if (local_verbose) print *, "ESM::SetModelServices:: calling NUOPC_DriverAddComp for OCN..."
    call NUOPC_DriverAddComp(driver, "OCN", ocnSS, comp=child, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompAttributeSet(child, name="Verbosity", value="1", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! Disabling the following macro, e.g. renaming to WITHCONNECTORS_disable,
    ! will result in a driver that does not call connectors between the model
    ! components. This mode can be used if all model components are driven 
    ! as independent models. However, even for independent models the
    ! connectors can be set here, but will turn into no-ops.
#define WITHCONNECTORS
#ifdef WITHCONNECTORS
    ! SetServices for atm2ocn
    if (local_verbose) print *, "ESM::SetModelServices:: calling NUOPC_DriverAddComp for ATM to OCN..."
    call NUOPC_DriverAddComp(driver, srcCompLabel="ATM", dstCompLabel="OCN", &
      compSetServicesRoutine=cplSS, comp=connector, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompAttributeSet(connector, name="Verbosity", value="1", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

!   call NUOPC_CompAttributeSet(comp=connector, name="ConnectionOptions", value=":remapmethod=redist", rc=rc)
!   call NUOPC_CompAttributeSet(comp=connector, name="cplList", value=":remapmethod=redist", rc=rc)
!   call NUOPC_CompAttributeGet(comp=connector, name="cplList", value=test_cplList, rc=rc)
!   print *, "cplList (NUOPC_CompAttributeGet) = "
!   print *, test_cplList
!   call NUOPC_CompAttributeGet(comp=connector, name="Verbosity", value=test_verbosity, rc=rc)
!   print *, "Verbosity = "
!   print *, test_verbosity
!   call NUOPC_ConnectorGet(connector=connector, srcFields, dstFields, rh, state, CplSet, cplSetList, srcVM, dstVM, rc)
!   call NUOPC_ConnectorGet(connector=connector, cplSetList=ptr_cplList, rc=rc)
!   call NUOPC_ConnectorSet(connector=connector, srcFields, dstFields, rh, state, CplSet, srcVM, dstVM, rc)
!   print *, "cplList (NUOPC_ConnectorGet) = "
!   print *, ptr_cplList

    ! SetServices for ocn2atm
    if (local_verbose) print *, "ESM::SetModelServices:: calling NUOPC_DriverAddComp for OCN to ATM..."
    call NUOPC_DriverAddComp(driver, srcCompLabel="OCN", dstCompLabel="ATM", &
      compSetServicesRoutine=cplSS, comp=connector, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompAttributeSet(connector, name="Verbosity", value="1", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

!   call NUOPC_CompAttributeSet(connector, name="cplList", value=":remapmethod=redist", rc=rc)
!   call NUOPC_CompAttributeSet(connector, name="ConnectionOptions", value=":remapmethod=redist", rc=rc)
   
      
#endif
      
    ! set the driver clock
!   call ESMF_TimeIntervalSet(timeStep, m=15, rc=rc) ! 15 minute steps
    if (local_verbose) print *, "ESM::SetModelServices:: calling ESMF_TimeIntervalSet for s=1..."
    call ESMF_TimeIntervalSet(timeStep, s=1, rc=rc) ! 15 minute steps
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (local_verbose) print *, "ESM::SetModelServices:: calling ESMF_TimeSet start yy=2010, mm=6, dd=1, h=0, m=0..."
    call ESMF_TimeSet(startTime, yy=2010, mm=6, dd=1, h=0, m=0, calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (local_verbose) print *, "ESM::SetModelServices:: calling ESMF_TimeSet stop yy=2010, mm=6, dd=1, h=1, m=0..."
    call ESMF_TimeSet(stopTime, yy=2010, mm=6, dd=1, h=1, m=0, calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (local_verbose) print *, "ESM::SetModelServices:: calling ESMF_ClockCreate..."
    internalClock = ESMF_ClockCreate(name="Application Clock", &
      timeStep=timeStep, startTime=startTime, stopTime=stopTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    if (local_verbose) print *, "ESM::SetModelServices:: calling ESMF_GridCompSet to set internal clock..."
    call ESMF_GridCompSet(driver, clock=internalClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
  end subroutine SetModelServices

  !-----------------------------------------------------------------------------

  subroutine ModifyCplLists(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc

    ! local variables
    character(len=160)              :: msg
    type(ESMF_CplComp), pointer     :: connectorList(:)
    integer                         :: i, j, cplListSize
    character(len=160), allocatable :: cplList(:)
    character(len=160)              :: tempString

    rc = ESMF_SUCCESS

    call ESMF_LogWrite("Driver is in ModifyCplLists()", ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    nullify(connectorList)
    call NUOPC_DriverGetComp(driver, compList=connectorList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    write (msg,*) "Found ", size(connectorList), " Connectors."// &
      " Modifying CplList Attribute...."
    call ESMF_LogWrite(trim(msg), ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    do i=1, size(connectorList)
      ! query the cplList for connector i
      call NUOPC_CompAttributeGet(connectorList(i), name="CplList", itemCount=cplListSize, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      if (cplListSize>0) then
        allocate(cplList(cplListSize))
        call NUOPC_CompAttributeGet(connectorList(i), name="CplList", valueList=cplList, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        ! go through all of the entries in the cplList
        do j=1, cplListSize
          if (trim(cplList(j))=="ocean_barotropic_streamfunction") then
            ! switch remapping to redist, b/c spectral model
            cplList(j) = trim(cplList(j))//":REMAPMETHOD=redist"
          elseif (trim(cplList(j))=="sea_water_temperature") then
            ! switch remapping to redist, b/c spectral model
            cplList(j) = trim(cplList(j))//":REMAPMETHOD=redist"
          elseif (trim(cplList(j))=="atmosphere_horizontal_streamfunction") then
            ! switch remapping to redist, b/c spectral model
            cplList(j) = trim(cplList(j))//":REMAPMETHOD=redist"
          elseif (trim(cplList(j))=="air_temperature") then
            ! switch remapping to redist, b/c spectral model
            cplList(j) = trim(cplList(j))//":REMAPMETHOD=redist"
          else
            print *, "ESM::ModifyCplLists:: Unsupported option: cplList(j) = ", cplList(j)
            print *, "ESM::ModifyCplLists:: changing to REMAPMETHOD=redist."
            ! switch remapping to redist, b/c spectral model
            cplList(j) = trim(cplList(j))//":REMAPMETHOD=redist"
          endif
        enddo
        ! store the modified cplList in CplList attribute of connector i
        call NUOPC_CompAttributeSet(connectorList(i), name="CplList", valueList=cplList, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        deallocate(cplList)
      endif
    enddo

    deallocate(connectorList)

  end subroutine ModifyCplLists

  !-----------------------------------------------------------------------------

end module ESM
