!==============================================================================
! Earth System Modeling Framework
! Copyright 2002-2019, University Corporation for Atmospheric Research, 
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics 
! Laboratory, University of Michigan, National Centers for Environmental 
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory, 
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

module MODEL

  !-----------------------------------------------------------------------------
  ! MODEL Component.
  !-----------------------------------------------------------------------------

  use maooam_coupled_wrapper, only: maooam_initialize, maooam_run, maooam_finalize
  use maooam_coupled_wrapper, only: maooam_natm, maooam_nocn

  use ESMF
  use NUOPC
  use NUOPC_Model, only: &
    model_routine_SS    => SetServices, &
    model_label_Advance => label_Advance
  
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

    logical :: local_verbose = .true.
    
    rc = ESMF_SUCCESS
    
    ! the NUOPC model component will register the generic methods
    if (local_verbose) print *, "SetServices:: Calling NUOPC_CompDerive..."
    call NUOPC_CompDerive(model, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (local_verbose) print *, "SetServices:: Finished NUOPC_CompDerive."
    
    ! set entry point for methods that require specific implementation
    if (local_verbose) print *, "SetServices:: Calling NUOPC_CompSetEntryPoint for p1..."
   call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, phaseLabelList=(/"IPDv00p1"/), userRoutine=InitializeP1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (local_verbose) print *, "SetServices:: Finished NUOPC_CompSetEntryPoint for p1."

    if (local_verbose) print *, "SetServices:: Calling NUOPC_CompSetEntryPoint for p2..."
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, phaseLabelList=(/"IPDv00p2"/), userRoutine=InitializeP2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (local_verbose) print *, "SetServices:: Finished NUOPC_CompSetEntryPoint for p2."
    
    ! attach specializing method(s)
    if (local_verbose) print *, "SetServices:: Calling NUOPC_CompSpecialize..."
    call NUOPC_CompSpecialize(model, specLabel=model_label_Advance, specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (local_verbose) print *, "SetServices:: Finished NUOPC_CompSpecialize."

    !STEVE: add a finalize routine 
!   call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_FINALIZE, phaseLabelList=(//), userRoutine=Finalize, rc=rc)
!   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!     line=__LINE__, &
!     file=__FILE__)) &
!     return  ! bail out

    print *, "SetServices:: Finished SetServices."
    
  end subroutine SetServices
  
  !-----------------------------------------------------------------------------

  subroutine InitializeP1(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    logical :: local_verbose = .true.
    
    rc = ESMF_SUCCESS

!   call NUOPC_FieldDictionarySetup(fileName='NUOPC_FieldDictionary.yaml', rc=rc)
!   STEVE: didn't work

    !--------------------------------------------------------------------------
    ! exportable field: atmosphere streamfunction
    !--------------------------------------------------------------------------
    call NUOPC_FieldDictionaryAddEntry(standardName='atmosphere_horizontal_streamfunction', canonicalUnits='m^2/s', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_Advertise(exportState, &
      StandardName="atmosphere_horizontal_streamfunction", name="psi", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !--------------------------------------------------------------------------
    ! exportable field: atmosphere temperature
    !--------------------------------------------------------------------------
    call NUOPC_FieldDictionaryAddEntry(standardName='air_temperature', canonicalUnits='K', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_Advertise(exportState, &
      StandardName="air_temperature", name="theta", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !--------------------------------------------------------------------------
    ! exportable field: ocean streamfunction
    !--------------------------------------------------------------------------
    call NUOPC_FieldDictionaryAddEntry(standardName='ocean_barotropic_streamfunction', canonicalUnits='m^2/s', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_Advertise(exportState, &
      StandardName="ocean_barotropic_streamfunction", name="A", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !--------------------------------------------------------------------------
    ! exportable field: ocean temperature
    !--------------------------------------------------------------------------
    call NUOPC_FieldDictionaryAddEntry(standardName='sea_water_temperature', canonicalUnits='m^2/s', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_Advertise(exportState, &
      StandardName="sea_water_temperature", name="T", rc=rc)
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
    ! Call model initialization routines, load initial conditions and assign default data
    !--------------------------------------------------------------------------
    print *, "InitializeP2:: Calling maooam_initialize..."
    call maooam_initialize()
    print *, "InitializeP2:: Finished maooam_initialize."
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
    ! create and open the config
    !--------------------------------------------------------------------------
!   config = ESMF_ConfigCreate(rc=rc)
!   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!     line=__LINE__, &
!     file=__FILE__)) &
!     return  ! bail out
!   call ESMF_ConfigLoadFile(config, "test.config", rc=rc)
!   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!     line=__LINE__, &
!     file=__FILE__)) &

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
    ! Create Field objects and 'realize' (i.e. assign to exportState)
    !--------------------------------------------------------------------------

    ! exportable array: Atmospheric streamfunction
    if (local_verbose) print *, "InitializeP2:: Calling ESMF_FieldCreate for psi..."
    field = ESMF_FieldCreate(grid=gridAtm, typekind=ESMF_TYPEKIND_R8, name="psi", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (local_verbose) print *, "InitializeP2:: Finished ESMF_FieldCreate for psi."
    call NUOPC_Realize(exportState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! exportable array: Atmospheric temperature
    if (local_verbose) print *, "InitializeP2:: Calling ESMF_FieldCreate for theta..."
    field = ESMF_FieldCreate(grid=gridAtm, typekind=ESMF_TYPEKIND_R8, name="theta", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (local_verbose) print *, "InitializeP2:: Finished ESMF_FieldCreate for theta."
    call NUOPC_Realize(state=exportState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! exportable array: Ocean streamfunction
    if (local_verbose) print *, "InitializeP2:: Calling ESMF_FieldCreate for A..."
    field = ESMF_FieldCreate(grid=gridOcn, typekind=ESMF_TYPEKIND_R8, name="A", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (local_verbose) print *, "InitializeP2:: Finished ESMF_FieldCreate for A."
    call NUOPC_Realize(state=exportState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! exportable array: Ocean temperature
    if (local_verbose) print *, "InitializeP2:: Calling ESMF_FieldCreate for T..."
    field = ESMF_FieldCreate(grid=gridOcn, typekind=ESMF_TYPEKIND_R8,  name="T", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (local_verbose) print *, "InitializeP2:: Finished ESMF_FieldCreate for T."
    call NUOPC_Realize(state=exportState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! fill export with some data for reference
    call ESMF_StateGet(exportState, itemName="psi", field=field, rc=rc)
    call ESMF_FieldGet(field, farrayPtr=dataPtr, rc=rc)
    dataPtr = 1.0d0

    call ESMF_StateGet(exportState, itemName="theta", field=field, rc=rc)
    call ESMF_FieldGet(field, farrayPtr=dataPtr, rc=rc)
    dataPtr = 2.0d0

    call ESMF_StateGet(exportState, itemName="A", field=field, rc=rc)
    call ESMF_FieldGet(field, farrayPtr=dataPtr, rc=rc)
    dataPtr = 3.0d0

    call ESMF_StateGet(exportState, itemName="T", field=field, rc=rc)
    call ESMF_FieldGet(field, farrayPtr=dataPtr, rc=rc)
    dataPtr = 4.0d0

    if (local_verbose) print *, "InitializeP2 :: finished."

  end subroutine InitializeP2
  
  !-----------------------------------------------------------------------------

  subroutine ModelAdvance(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    character(len=160)          :: msgString

    !STEVE: ESMF time information for the model
    type(ESMF_Time)             :: startTime
    type(ESMF_Time)             :: currTime
    type(ESMF_TimeInterval)     :: timeStep

    type(ESMF_StateItem_Flag)   :: itemType
    type(ESMF_Field)        :: field
    real(ESMF_KIND_R8), pointer :: farrayPtr(:) => null()
    real(ESMF_KIND_R8), pointer :: dataPtr(:)   => null()
    real(ESMF_KIND_R8), pointer :: dataPtr_psi(:)   => null()
    real(ESMF_KIND_R8), pointer :: dataPtr_theta(:)   => null()
    real(ESMF_KIND_R8), pointer :: dataPtr_A(:)   => null()
    real(ESMF_KIND_R8), pointer :: dataPtr_T(:)   => null()
    real(kind=8) :: X0,Xf
    integer(kind=8) :: seconds
    real(kind=8) :: t,dt
    integer :: Nt
    integer :: ndim
    integer :: j

    logical :: local_verbose = .true.

    if (local_verbose) print *, "ModelAdvance :: commencing..."

    rc = ESMF_SUCCESS
   
    ! 17.2.3 Implement a user-code Run routine
    ! During the execution loop, the run routine may be called many times. 
    ! Each time it should read data from the importState, use the clock to 
    ! determine what the current time is in the calling component, compute 
    ! new values or process the data, and produce any output and place it 
    ! in the exportState.
    
    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(model, clock=clock, importState=importState, exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep
    
    ! Because of the way that the internal Clock was set by default,
    ! its timeStep is equal to the parent timeStep. As a consequence the
    ! currTime + timeStep is equal to the stopTime of the internal Clock
    ! for this call of the ModelAdvance() routine.

    call ESMF_ClockPrint(clock, options="currTime", &
      preString="---->Advancing Model from: ", unit=msgString, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockPrint(clock, options="stopTime", &
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
    ! Get model input data
    !--------------------------------------------------------------------------

    ndim = 2*maooam_natm+2*maooam_nocn
    allocate(farrayPtr(0:ndim))    ! user controlled allocation
    farrayPtr(0) = f0

    ! call ESMF_StateGet(), etc to get fields, bundles, arrays from import state.
    call ESMF_StateGet(exportState, itemName="psi", itemType=itemType, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (itemType /= ESMF_STATEITEM_NOTFOUND) then
      call ESMF_StateGet(exportState, itemName="psi", field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_FieldGet(field, farrayPtr=dataPtr_psi, rc=rc) 
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (local_verbose) print *, "dataPtr (psi) = "
    if (local_verbose) print *, dataPtr_psi
!   farrayPtr(si:ei) => dataPtr
    farrayPtr(si(1):ei(1)) = dataPtr_psi
    if (local_verbose) print *, "farrayPtr(",si(1),":",ei(1),") = "
    if (local_verbose) print *, farrayPtr(si(1):ei(1))

    call ESMF_StateGet(exportState, itemName="theta", itemType=itemType, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (itemType /= ESMF_STATEITEM_NOTFOUND) then
      call ESMF_StateGet(exportState, itemName="theta", field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_FieldGet(field, farrayPtr=dataPtr_theta, rc=rc) 
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (local_verbose) print *, "dataPtr (theta) = "
    if (local_verbose) print *, dataPtr_theta
!   farrayPtr(si:ei) => dataPtr
    farrayPtr(si(2):ei(2)) = dataPtr_theta
    if (local_verbose) print *, "farrayPtr(",si(2),":",ei(2),") = "
    if (local_verbose) print *, farrayPtr(si(2):ei(2))

    call ESMF_StateGet(exportState, itemName="A", itemType=itemType, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (itemType /= ESMF_STATEITEM_NOTFOUND) then
      call ESMF_StateGet(exportState, itemName="A", field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_FieldGet(field, farrayPtr=dataPtr_A, rc=rc) 
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (local_verbose) print *, "dataPtr (A) = "
    if (local_verbose) print *, dataPtr_A
!   farrayPtr(si:ei) => dataPtr
    farrayPtr(si(3):ei(3)) = dataPtr_A
    if (local_verbose) print *, "farrayPtr(",si(3),":",ei(3),") = "
    if (local_verbose) print *, farrayPtr(si(4):ei(4))

    call ESMF_StateGet(exportState, itemName="T", itemType=itemType, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (itemType /= ESMF_STATEITEM_NOTFOUND) then
      call ESMF_StateGet(exportState, itemName="T", field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_FieldGet(field, farrayPtr=dataPtr_T, rc=rc) 
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    if (local_verbose) print *, "dataPtr (T) = "
    if (local_verbose) print *, dataPtr_T
!   farrayPtr(si:ei) => dataPtr
    farrayPtr(si(4):ei(4)) = dataPtr_T
    if (local_verbose) print *, "farrayPtr(",si(4),":",ei(4),") = "
    if (local_verbose) print *, farrayPtr(si(4):ei(4))

!   print *, "farrayPtr = "
    if (local_verbose) print *, "farrayPtr(",si(1),":",ei(4),") = "
    if (local_verbose) print *, farrayPtr(si(1):ei(4))

    !----------------------------
    ! Model timing info
    !----------------------------
    call ESMF_ClockGet(clock, startTime=startTime, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_TimeGet(currTime, s_i8=seconds, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return

    t = dble(seconds)
    if (local_verbose) print *, "t seconds = ", seconds

    call ESMF_TimeIntervalGet(timeStep, s_i8=seconds, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return

    dt = dble(seconds)
    if (local_verbose) print *, "dt seconds = ", seconds

    Nt = 1 !STEVE: just run one step of dt

!   print *, "ModelAdvance:: Pre- maooam model run:  farrayP = "
!   print *, farrayP            ! print PET-local farrayA directly
    ! Transform seconds to model non-dimensional time (approximate)
    t = t/1000
    dt = dt/1000
    if (local_verbose) print *, "Using t = ", t
    if (local_verbose) print *, "Using dt = ", dt

    call maooam_run(X=farrayPtr,t=t,dt=dt,Nt=Nt) !,component)

    ! Fill export state here using ESMF_StateAdd(), etc
    f0 = farrayPtr(0)
    dataPtr_psi = farrayPtr(si(1):ei(1))
    dataPtr_theta = farrayPtr(si(2):ei(2))
    dataPtr_A = farrayPtr(si(3):ei(3))
    dataPtr_T = farrayPtr(si(4):ei(4))

!   print *, "ModelAdvance:: Post- maooam model run: farrayPtr = "
!   print *, farrayPtr       ! print PET-local farrayA directly
    
    if (local_verbose) print *, "ModelAdvance:: finished."
    
  end subroutine ModelAdvance

  subroutine Finalize(model,importState,exportState,clock,rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    !STEVE: ESMF pointer to store array data (i.e. model state)
    type(ESMF_DistGrid)         :: distgrid       ! DistGrid object
    type(ESMF_Field)            :: field          ! Field object
    type(ESMF_Grid)             :: gridOut

    ! Clean up
    call ESMF_GridCompFinalize(model)
    
  end subroutine Finalize

end module MODEL
