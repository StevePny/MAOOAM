!==============================================================================
! Earth System Modeling Framework
! Copyright 2002-2019, University Corporation for Atmospheric Research, 
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics 
! Laboratory, University of Michigan, National Centers for Environmental 
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory, 
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================

module ATM

  !-----------------------------------------------------------------------------
  ! ATM Component.
  !-----------------------------------------------------------------------------

  use maooam_atmos_wrapper, only: maooam_atmos_initialize, maooam_atmos_run, maooam_atmos_finalize
  use maooam_atmos_wrapper, only: maooam_natm, maooam_nocn

  use ESMF
  use NUOPC
  use NUOPC_Model, &
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
    
    rc = ESMF_SUCCESS
    
    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(model, model_routine_SS, rc=rc)
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
    call NUOPC_CompSpecialize(model, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
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

    !-------------------------------------------------------------------------- 
    ! importable field: ocean_barotropic_streamfunction
    !-------------------------------------------------------------------------- 
    if (local_Verbose) print *, "ATM::InitializeP1:: calling NUOPC_FieldDictionaryAddEntry..."
    call NUOPC_FieldDictionaryAddEntry(standardName='ocean_barotropic_streamfunction', canonicalUnits='m^2/s', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (local_verbose) print *, "ATM::InitializeP1:: Calling NUOPC_Advertise for ocean_barotropic_streamfunction..."
    call NUOPC_Advertise(importState, StandardName="ocean_barotropic_streamfunction", name="A", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !-------------------------------------------------------------------------- 
    ! importable field: sea_water_temperature
    !-------------------------------------------------------------------------- 
    if (local_Verbose) print *, "ATM::InitializeP1:: calling NUOPC_FieldDictionaryAddEntry..."
    call NUOPC_FieldDictionaryAddEntry(standardName='sea_water_temperature', canonicalUnits='m^2/s', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (local_verbose) print *, "ATM::InitializeP1:: Calling NUOPC_Advertise for sea_water_temperature..."
    call NUOPC_Advertise(importState, StandardName="sea_water_temperature", name="T", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !-------------------------------------------------------------------------- 
    ! exportable field: atmosphere_horizontal_streamfunction
    !-------------------------------------------------------------------------- 
    if (local_Verbose) print *, "ATM::InitializeP1:: calling NUOPC_FieldDictionaryAddEntry..."
    call NUOPC_FieldDictionaryAddEntry(standardName='atmosphere_horizontal_streamfunction', canonicalUnits='m^2/s', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (local_verbose) print *, "ATM::InitializeP1:: Calling NUOPC_Advertise for atmosphere_horizontal_streamfunction..."
    call NUOPC_Advertise(exportState, StandardName="atmosphere_horizontal_streamfunction", name="psi", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !-------------------------------------------------------------------------- 
    ! exportable field: air_temperature
    !-------------------------------------------------------------------------- 
    if (local_Verbose) print *, "ATM::InitializeP1:: calling NUOPC_FieldDictionaryAddEntry..."
    call NUOPC_FieldDictionaryAddEntry(standardName='air_temperature', canonicalUnits='K', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    if (local_verbose) print *, "ATM::InitializeP1:: Calling NUOPC_Advertise for air_temperature..."
    call NUOPC_Advertise(exportState, StandardName="air_temperature", name="theta", rc=rc)
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
    ! Call model initialization routines, load initial conditions and assign to farrayP
    !--------------------------------------------------------------------------
    print *, "ATM::InitializeP2:: Calling maooam_atmos_initialize..."
    call maooam_atmos_initialize()
    print *, "ATM::InitializeP2:: Finished maooam_atmos_initialize."
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

    ! importable array: Ocean streamfunction
    if (local_verbose) print *, "ATM::InitializeP2:: Calling ESMF_FieldCreate for A..."
    field = ESMF_FieldCreate(grid=gridOcn, typekind=ESMF_TYPEKIND_R8, name="A", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (local_verbose) print *, "ATM::InitializeP2:: Finished ESMF_FieldCreate for A."
    call NUOPC_Realize(state=importState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! importable array: Ocean temperature
    if (local_verbose) print *, "ATM::InitializeP2:: Calling ESMF_FieldCreate for T..."
    field = ESMF_FieldCreate(grid=gridOcn, typekind=ESMF_TYPEKIND_R8, name="T", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (local_verbose) print *, "ATM::InitializeP2:: Finished ESMF_FieldCreate for T."
    call NUOPC_Realize(state=importState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !--------------------------------------------------------------------------
    ! Get the exportable arrays
    !--------------------------------------------------------------------------

    ! exportable array: Atmospheric streamfunction
    if (local_verbose) print *, "ATM::InitializeP2:: Calling ESMF_FieldCreate for psi..."
    field = ESMF_FieldCreate(grid=gridAtm, typekind=ESMF_TYPEKIND_R8, name="psi", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (local_verbose) print *, "ATM::InitializeP2:: Finished ESMF_FieldCreate for psi."
    call NUOPC_Realize(state=exportState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! exportable array: Atmospheric temperature
    if (local_verbose) print *, "ATM::InitializeP2:: Calling ESMF_FieldCreate for theta..."
    field = ESMF_FieldCreate(grid=gridAtm, typekind=ESMF_TYPEKIND_R8, name="theta", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (local_verbose) print *, "ATM::InitializeP2:: Finished ESMF_FieldCreate for theta."
    call NUOPC_Realize(exportState, field=field, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !--------------------------------------------------------------------------
    ! fill export with some data for reference
    !--------------------------------------------------------------------------
    call ESMF_StateGet(exportState, itemName="psi", field=field, rc=rc)
    call ESMF_FieldGet(field, farrayPtr=dataPtr, rc=rc)
    dataPtr = 1.0d0

    call ESMF_StateGet(exportState, itemName="theta", field=field, rc=rc)
    call ESMF_FieldGet(field, farrayPtr=dataPtr, rc=rc)
    dataPtr = 2.0d0

    if (local_verbose) print *, "ATM::InitializeP2:: finished."

  end subroutine InitializeP2
  
  !-----------------------------------------------------------------------------

  subroutine ModelAdvance(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(ESMF_Time)             :: currTime, startTime
    type(ESMF_TimeInterval)     :: timeStep
    character(len=160)          :: msgString

    type(ESMF_StateItem_Flag)   :: itemType
    type(ESMF_Field)        :: field
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
    logical :: local_writeout = .true.
    character(14) :: outfile = 'evol_atmos.dat'
    integer :: iunit
    integer :: fid = 10

    if (local_verbose) print *, "ATM::ModelAdvance :: commencing..."

    rc = ESMF_SUCCESS
    
    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(model, modelClock=clock, importState=importState, &
      exportState=exportState, rc=rc)
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
      preString="------>Advancing ATM from: ", unit=msgString, rc=rc)
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
    ! Set up data array farrayPtr
    !--------------------------------------------------------------------------
    ndim = 2*maooam_natm+2*maooam_nocn
    print *, "ndim = ", ndim
    allocate(farrayPtr(0:ndim))    ! user controlled allocation

    !--------------------------------------------------------------------------
    ! Set up data array farrayPtr
    !--------------------------------------------------------------------------
    call getDataPtr(exportState,itemName='psi',dataPtr=dataPtr_psi)
    call getDataPtr(exportState,itemName='theta',dataPtr=dataPtr_theta)
    call getDataPtr(importState,itemName='A',dataPtr=dataPtr_A)
    call getDataPtr(importState,itemName='T',dataPtr=dataPtr_T)

    farrayPtr(0) = f0
    farrayPtr(si(1):ei(1)) = dataPtr_psi
    farrayPtr(si(2):ei(2)) = dataPtr_theta
    farrayPtr(si(3):ei(3)) = dataPtr_A
    farrayPtr(si(4):ei(4)) = dataPtr_T
    if (local_verbose) print *, "farrayPtr = "
    if (local_verbose) print *, farrayPtr

    !--------------------------------------------------------------------------
    ! Get timing info
    !--------------------------------------------------------------------------
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
    call ESMF_TimeIntervalGet(timeStep, s_i8=seconds, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return
    dt = dble(seconds)
    Nt = 1 !STEVE: just run one step of dt

    t = t/1000
    dt = dt/1000
    if (local_verbose) print *, "Using t = ", t
    if (local_verbose) print *, "Using dt = ", dt

    !--------------------------------------------------------------------------
    ! Step the model
    !--------------------------------------------------------------------------
    if (local_verbose) print *, "ATM::ModelAdvance:: calling maooam_atmos_run..."
    call maooam_atmos_run(X=farrayPtr,t=t,dt=dt,Nt=Nt) !,component)

    ! Write to file
    if (local_writeout) then
      inquire(file=outfile, number=iunit)
      if (iunit < 0) open(fid,file=outfile)
      write(fid,*) t, farrayPtr(1:ndim)
    endif

    ! Only update the atmospheric field
    f0 = farrayPtr(0)
    dataPtr_psi = farrayPtr(si(1):ei(1))
    dataPtr_theta = farrayPtr(si(2):ei(2))

  end subroutine ModelAdvance

  subroutine getDataPtr(state,itemName,dataPtr)
    type(ESMF_State), intent(in)            :: state
    character(*), intent(in)                :: itemName
    real(ESMF_KIND_R8), pointer, intent(in) :: dataPtr(:)

    integer                     :: rc
    type(ESMF_Field)            :: field
    type(ESMF_StateItem_Flag)   :: itemType
    character(len=160)          :: msgString

    logical :: local_verbose = .true.

    if (local_verbose) print *, "ATM::getDataPtr ..."

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
    else
      print *, "ATM::getDataPtr:: importState ESMF_STATEITEM_NOTFOUND : ", trim(itemName)
      stop 8
    endif
    if (local_verbose) print *, "ATM::dataPtr (",trim(itemName),")"," = "
    if (local_verbose) print *, dataPtr

  end subroutine getDataPtr

end module ATM
