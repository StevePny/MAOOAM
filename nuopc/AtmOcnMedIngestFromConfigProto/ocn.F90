!==============================================================================
! Earth System Modeling Framework
! Copyright 2002-2019, University Corporation for Atmospheric Research, 
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics 
! Laboratory, University of Michigan, National Centers for Environmental 
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory, 
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

module OCN

  !-----------------------------------------------------------------------------
  ! OCN Component.
  !-----------------------------------------------------------------------------

  use maooam_ocean_wrapper, only: maooam_ocean_initialize, maooam_ocean_run, maooam_ocean_finalize
  use maooam_ocean_wrapper, only: maooam_nocn, maooam_natm

  use ESMF
  use NUOPC
  use NUOPC_Model, inheritModel => SetServices
  
  implicit none
  
  private

  real(kind=8), save :: f0=0
  integer, dimension(4) :: si, ei

  public SetServices
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  
  subroutine SetServices(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc
    
    rc = ESMF_SUCCESS
    
    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(model, inheritModel, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p1"/), userRoutine=InitializeP1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p2"/), userRoutine=InitializeP2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! attach specializing method(s)
!   call NUOPC_CompSpecialize(model, specLabel=label_SetClock, &
!     specRoutine=SetClock, rc=rc)
!   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!     line=__LINE__, &
!     file=__FILE__)) &
!     return  ! bail out

    call NUOPC_CompSpecialize(model, specLabel=label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! attach specializing method(s)
    ! -> NUOPC specializes by default --->>> first need to remove the default
    call ESMF_MethodRemove(model, label_CheckImport, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call NUOPC_CompSpecialize(model, specLabel=label_CheckImport, &
      specRoutine=CheckImport, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine InitializeP1(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    logical :: local_verbose = .true.
    
    rc = ESMF_SUCCESS

    if (local_verbose) print *, "OCN::InitializeP1:: Calling NUOPC_Advertise for atmosphere_horizontal_streamfunction..."
    call NUOPC_Advertise(importState, StandardName="atmosphere_horizontal_streamfunction", name="psi", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (local_verbose) print *, "OCN::InitializeP1:: Calling NUOPC_Advertise for air_temperature..."
    call NUOPC_Advertise(importState, StandardName="air_temperature", name="theta", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (local_verbose) print *, "OCN::InitializeP1:: Calling NUOPC_Advertise for ocean_barotropic_streamfunction..."
    call NUOPC_Advertise(exportState, StandardName="ocean_barotropic_streamfunction", name="A", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (local_verbose) print *, "OCN::InitializeP1:: Calling NUOPC_Advertise for sea_water_temperature..."
    call NUOPC_Advertise(exportState, StandardName="sea_water_temperature", name="T", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine InitializeP1
  
  !-----------------------------------------------------------------------------

  subroutine InitializeP2(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! local variables    
    type(ESMF_TimeInterval) :: stabilityTimeStep
    type(ESMF_Field)        :: field
    type(ESMF_Grid)         :: gridAtm
    type(ESMF_Grid)         :: gridOcn

    !STEVE: ESMF pointer to store grid
    type(ESMF_DistGrid)         :: distgridAtm       ! DistGrid object
    type(ESMF_DistGrid)         :: distgridOcn       ! DistGrid object

    real(ESMF_KIND_R8), pointer   :: dataPtr(:) => null()

    logical :: local_verbose = .true.

    rc = ESMF_SUCCESS

    !--------------------------------------------------------------------------
    ! Call model initialization routines, load initial conditions and assign to farrayP
    !--------------------------------------------------------------------------
    print *, "OCN::InitializeP2:: Calling maooam_ocean_initialize..."
    !STEVE: have to bypass this since it's identical to the atmosphere initialization
    !       and tries to reallocate the same arrays
    call maooam_ocean_initialize()
    print *, "OCN::InitializeP2:: Finished maooam_ocean_initialize."
    ! Set up grid indices:
    si(1) = 1
    ei(1) = maooam_natm
    si(2) = maooam_natm + 1
    ei(2) = maooam_natm + maooam_natm
    si(3) = maooam_natm*2 + 1
    ei(3) = maooam_natm*2 + maooam_nocn
    si(4) = maooam_natm*2 + maooam_nocn + 1
    ei(4) = maooam_natm*2 + maooam_nocn + maooam_nocn

    !--------------------------------------------------------------------------
    ! Set up ESMF grid objects
    !--------------------------------------------------------------------------

    ! create a distribution Grid object
    distgridAtm = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/maooam_natm/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Create the Grid objects
    gridAtm = ESMF_GridCreate(distgrid=distgridAtm, name="atmos_grid", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! create a distribution Grid object
    distgridOcn = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/maooam_nocn/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Create the Grid objects
    gridOcn = ESMF_GridCreate(distgrid=distgridOcn, name="ocean_grid", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !--------------------------------------------------------------------------
    ! Get the importable arrays
    !--------------------------------------------------------------------------

    ! importable array: Atmospheric streamfunction
    if (local_verbose) print *, "OCN::InitializeP2:: Calling ESMF_FieldCreate for psi..."
    field = ESMF_FieldCreate(grid=gridAtm, typekind=ESMF_TYPEKIND_R8, name="psi", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (local_verbose) print *, "OCN::InitializeP2:: Finished ESMF_FieldCreate for psi."
    call NUOPC_Realize(state=importState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! importable array: Atmospheric temperature
    if (local_verbose) print *, "OCN::InitializeP2:: Calling ESMF_FieldCreate for theta..."
    field = ESMF_FieldCreate(grid=gridAtm, typekind=ESMF_TYPEKIND_R8, name="theta", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (local_verbose) print *, "OCN::InitializeP2:: Finished ESMF_FieldCreate for theta."
    call NUOPC_Realize(importState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !--------------------------------------------------------------------------
    ! Get the exportable arrays
    !--------------------------------------------------------------------------

    ! exportable array: Ocean streamfunction
    if (local_verbose) print *, "OCN::InitializeP2:: Calling ESMF_FieldCreate for A..."
    field = ESMF_FieldCreate(grid=gridOcn, typekind=ESMF_TYPEKIND_R8, name="A", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (local_verbose) print *, "OCN::InitializeP2:: Finished ESMF_FieldCreate for A."
    call NUOPC_Realize(state=exportState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! exportable array: Ocean temperature
    if (local_verbose) print *, "OCN::InitializeP2:: Calling ESMF_FieldCreate for T..."
    field = ESMF_FieldCreate(grid=gridOcn, typekind=ESMF_TYPEKIND_R8, name="T", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (local_verbose) print *, "OCN::InitializeP2:: Finished ESMF_FieldCreate for T."
    call NUOPC_Realize(state=exportState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !--------------------------------------------------------------------------
    ! fill export with some data for reference
    !--------------------------------------------------------------------------
    call ESMF_StateGet(exportState, itemName="A", field=field, rc=rc)
    call ESMF_FieldGet(field, farrayPtr=dataPtr, rc=rc)
    dataPtr = 3.0d0

    call ESMF_StateGet(exportState, itemName="T", field=field, rc=rc)
    call ESMF_FieldGet(field, farrayPtr=dataPtr, rc=rc)
    dataPtr = 4.0d0

    if (local_verbose) print *, "OCN::InitializeP2:: finished."

  end subroutine InitializeP2
  
  !-----------------------------------------------------------------------------

  subroutine SetClock(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)              :: clock
    type(ESMF_TimeInterval)       :: stabilityTimeStep

    rc = ESMF_SUCCESS
    
    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(model, modelClock=clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! initialize internal clock
    ! here: parent Clock and stability timeStep determine actual model timeStep
    !TODO: stabilityTimeStep should be read in from configuation
    !TODO: or computed from internal Grid information
    call ESMF_TimeIntervalSet(stabilityTimeStep, s=10, rc=rc) ! 10 second steps
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_CompSetClock(model, clock, stabilityTimeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine SetClock

  !-----------------------------------------------------------------------------

  subroutine ModelAdvance(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(ESMF_Time)             :: currTime
    type(ESMF_TimeInterval)     :: timeStep
    character(len=160)          :: msgString

    type(ESMF_StateItem_Flag)   :: itemType
    type(ESMF_Field)            :: field
    real(ESMF_KIND_R8), pointer :: farrayPtr(:) => null()
    real(ESMF_KIND_R8), pointer :: dataPtr_psi(:)   => null()
    real(ESMF_KIND_R8), pointer :: dataPtr_theta(:)   => null()
    real(ESMF_KIND_R8), pointer :: dataPtr_A(:)   => null()
    real(ESMF_KIND_R8), pointer :: dataPtr_T(:)   => null()

    integer(kind=8) :: seconds
    real(kind=8) :: t,dt
    integer :: Nt
    integer :: ndim

    logical :: local_verbose = .true.

    if (local_verbose) print *, "OCN::ModelAdvance:: commencing..."

    rc = ESMF_SUCCESS
    
    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(model, modelClock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep
    
    ! Because of the way that the internal Clock was set in SetClock(),
    ! its timeStep is likely smaller than the parent timeStep. As a consequence
    ! the time interval covered by a single parent timeStep will result in 
    ! multiple calls to the ModelAdvance() routine. Every time the currTime
    ! will come in by one internal timeStep advanced. This goes until the
    ! stopTime of the internal Clock has been reached.
    
    call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing OCN from: ", unit=msgString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    call ESMF_TimePrint(currTime + timeStep, &
      preString="---------------------> to: ", unit=msgString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !--------------------------------------------------------------------------
    ! Setup model input pointer to state data
    !--------------------------------------------------------------------------
    ndim = 2*maooam_natm+2*maooam_nocn
    print *, "ndim = ", ndim
    allocate(farrayPtr(0:ndim))    ! user controlled allocation

    !--------------------------------------------------------------------------
    ! Set up data array farrayPtr
    !--------------------------------------------------------------------------
    call getDataPtr(importState,itemName='psi',dataPtr=dataPtr_psi)
    call getDataPtr(importState,itemName='theta',dataPtr=dataPtr_theta)
    call getDataPtr(exportState,itemName='A',dataPtr=dataPtr_A)
    call getDataPtr(exportState,itemName='T',dataPtr=dataPtr_T)

    farrayPtr(0) = f0
    farrayPtr(si(1):ei(1)) = dataPtr_psi
    farrayPtr(si(2):ei(2)) = dataPtr_theta
    farrayPtr(si(3):ei(3)) = dataPtr_A
    farrayPtr(si(4):ei(4)) = dataPtr_T

    !--------------------------------------------------------------------------
    ! Run model integration/step
    !--------------------------------------------------------------------------

    if (local_verbose) print *, "OCN::ModelAdvance:: calling ESMF_TimeGet..."
    call ESMF_TimeGet(currTime, s_i8=seconds, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return
    t = real(seconds)
    if (local_verbose) print *, "OCN::ModelAdvance:: calling ESMF_TimeIntervalGet..."
    call ESMF_TimeIntervalGet(timeStep, s_i8=seconds, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return
    dt = real(seconds)
    Nt = 1 !STEVE: just run one step of dt

    !--------------------------------------------------------------------------
    ! Step the model
    !--------------------------------------------------------------------------
    t = t/1000
    dt = dt/1000
    print *, "Using t = ", t
    print *, "Using dt = ", dt
    if (local_verbose) print *, "OCN::ModelAdvance:: calling maooam_ocean_run..."

    call maooam_ocean_run(X=farrayPtr,t=t,dt=dt,Nt=Nt) !,component)

    ! Only update the ocean field
    f0 = farrayPtr(0)
    dataPtr_A = farrayPtr(si(3):ei(3))
    dataPtr_T = farrayPtr(si(4):ei(4))
    deallocate(farrayPtr)

  end subroutine ModelAdvance

  subroutine getDataPtr(state,itemName,dataPtr)
    type(ESMF_State), intent(inout)         :: state
    character(*), intent(in)                :: itemName
    real(ESMF_KIND_R8), pointer, intent(in) :: dataPtr(:)

    integer                     :: rc
    type(ESMF_Field)            :: field
    type(ESMF_StateItem_Flag)   :: itemType
    character(len=160)          :: msgString

    logical :: local_verbose = .true.

    if (local_verbose) print *, "OCN::getDataPtr ..."
    
    call ESMF_StateGet(state, itemName=itemName, itemType=itemType, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (itemType /= ESMF_STATEITEM_NOTFOUND) then
      call ESMF_StateGet(state, itemName=itemName, field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_FieldGet(field, farrayPtr=dataPtr, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif

    if (local_verbose) print *, "OCN::dataPtr (",trim(itemName),")"," = "
    if (local_verbose) print *, dataPtr

  end subroutine getDataPtr

  !-----------------------------------------------------------------------------

  subroutine CheckImport(model, rc)
    type(ESMF_GridComp)   :: model
    integer, intent(out)  :: rc
    
    ! local variables
    type(ESMF_Clock)                :: clock
    type(ESMF_Time)                 :: currTime, invalidTime
    type(ESMF_State)                :: importState
    logical                         :: timeCheck
    type(ESMF_Field),       pointer :: fieldList(:)
    integer                         :: i
    
    rc = ESMF_SUCCESS

    ! query the Component for its clock and importState
    call ESMF_GridCompGet(model, clock=clock, importState=importState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! get the current time out of the clock
    call ESMF_ClockGet(clock, currTime=currTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! set up invalid time (by convention)
    call ESMF_TimeSet(invalidTime, yy=99999999, mm=01, dd=01, h=00, m=00, s=00, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! Loop through all the field in the importState, and test whether they
    ! are at invalidTime (ignore them for now), or at currTime. Any other
    ! time coming in would flag an incompatibility.
    
    nullify(fieldList)
    call NUOPC_GetStateMemberLists(importState, fieldList=fieldList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    do i=1, size(fieldList)

      timeCheck = NUOPC_IsAtTime(fieldList(i), invalidTime, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      if (timeCheck) then
        ! The field is at invalidTime
        ! -> In a real application mark the field with a flag as invalid 
        !    so the actual model code can act accordingly.
        ! Here for purpose of demonstration just log a message and continue on.
      
        call ESMF_LogWrite("OCN: detected import field at invalidTime", &
          ESMF_LOGMSG_INFO, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      else
        ! The field is NOT at invalidTime -> it must then be at currTime or it
        ! is incompatible!
      
        ! check that Fields in the importState show correct timestamp
        timeCheck = NUOPC_IsAtTime(fieldList(i), currTime, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out

        if (.not.timeCheck) then
          !TODO: introduce and use INCOMPATIBILITY return codes!!!!
          call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
            msg="NUOPC INCOMPATIBILITY DETECTED: "//&
            "Import Field not at current time", &
            line=__LINE__, file=__FILE__, &
            rcToReturn=rc)
          return  ! bail out
        endif
      
      endif
    enddo
    
  end subroutine
    
  !-----------------------------------------------------------------------------

end module OCN
