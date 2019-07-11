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

  real(ESMF_KIND_R8), pointer :: farrayP(:)   ! Fortran array pointer
  integer, parameter :: fdim2 = 2
  
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

    !STEVE: testing IPDv05p3
!   call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
!     phaseLabelList=(/"IPDv05p3"/), userRoutine=InitializeP3, rc=rc)
!   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!     line=__LINE__, &
!     file=__FILE__)) &
    
    ! attach specializing method(s)
    call NUOPC_CompSpecialize(model, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine SetServices
  
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
    integer :: ndim, si, ei   ! ndim = model dimension, si = start index, ei = end index

    logical :: local_verbose = .true.
    
    rc = ESMF_SUCCESS

    !--------------------------------------------------------------------------
    ! Call model initialization routines, load initial conditions and assign to farrayP
    !--------------------------------------------------------------------------
    print *, "ATM::InitializeP2:: Calling maooam_atmos_initialize..."
    call maooam_atmos_initialize()
    print *, "ATM::InitializeP2:: Finished maooam_atmos_initialize."

    ! Get model state dimension
    ndim = 2*maooam_natm+2*maooam_nocn

    if (ndim == 0 .or. ndim <= maooam_natm .or. ndim <= maooam_nocn) then
      print *, "ATM::InitializeP2::ERROR:: ndim = ", 0
      print *, "ATM::InitializeP2::EXITING..."
      stop 'ATM'
    endif

    ! Allocate pointer:
    print *, "ndim = ", ndim
    print *, "allocating farrayP(ndim)..."
    allocate(farrayP(ndim))    ! user controlled allocation
    farrayP(1:10)  = 1.0d0            ! initialize to some value
    farrayP(11:20) = 2.0d0            ! initialize to some value
    farrayP(21:28) = 3.0d0            ! initialize to some value
    farrayP(29:36) = 4.0d0            ! initialize to some value
    print *, "farrayP = ", farrayP

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
!   gridAtm = ESMF_GridCreateNoPeriDim(minIndex=(/1/), maxIndex=(/maooam_natm/), name="atmos_grid", rc=rc)
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
!   gridOcn = ESMF_GridCreateNoPeriDim(minIndex=(/1), maxIndex=(/maooam_nocn/), name="ocean_grid", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

!   call ESMF_GridAddCoord(gridAtm, staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)
!   call ESMF_GridAddCoord(gridOcn, staggerLoc=ESMF_STAGGERLOC_CENTER, rc=rc)

    !--------------------------------------------------------------------------
    ! Get the importable arrays
    !--------------------------------------------------------------------------

    ! importable array: Ocean temperature
    si = maooam_natm*2 + 1
    ei = maooam_natm*2 + maooam_nocn
    print *, "si, ei = ", si, ei
    if (local_verbose) print *, "ATM::InitializeP2:: Calling ESMF_FieldCreate for T..."
    field = ESMF_FieldCreate(grid=gridOcn, farray=farrayP(si:ei), indexflag=ESMF_INDEX_DELOCAL,  name="T", rc=rc)
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

    ! importable array: Ocean streamfunction
    si = maooam_natm*2 + maooam_nocn + 1
    ei = maooam_natm*2 + maooam_nocn + maooam_nocn
    print *, "si, ei = ", si, ei
    if (local_verbose) print *, "ATM::InitializeP2:: Calling ESMF_FieldCreate for A..."
    field = ESMF_FieldCreate(grid=gridOcn, farray=farrayP(si:ei), indexflag=ESMF_INDEX_DELOCAL, name="A", rc=rc)
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

    !--------------------------------------------------------------------------
    ! Get the exportable arrays
    !--------------------------------------------------------------------------

    ! exportable array: Atmospheric temperature
    si = 1
    ei = maooam_natm
    print *, "si, ei = ", si, ei
    if (local_verbose) print *, "ATM::InitializeP2:: Calling ESMF_FieldCreate for theta..."
    field = ESMF_FieldCreate(grid=gridAtm, farray=farrayP(si:ei), indexflag=ESMF_INDEX_DELOCAL, name="theta", rc=rc)
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

    ! exportable array: Atmospheric streamfunction
    si = maooam_natm + 1
    ei = maooam_natm + maooam_natm
    print *, "si, ei = ", si, ei
    if (local_verbose) print *, "ATM::InitializeP2:: Calling ESMF_FieldCreate for psi..."
    field = ESMF_FieldCreate(grid=gridAtm, farray=farrayP(si:ei), indexflag=ESMF_INDEX_DELOCAL, name="psi", rc=rc)
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

    real(kind=8), dimension(:), pointer :: X
    integer(kind=8) :: seconds
    real(kind=8) :: t,dt
    integer :: Nt

    logical :: local_verbose = .true.

    if (local_verbose) print *, "ATM::ModelAdvance :: commencing..."

#define NUOPC_TRACE__OFF
#ifdef NUOPC_TRACE
    call ESMF_TraceRegionEnter("ATM:ModelAdvance")
#endif

    rc = ESMF_SUCCESS
 
    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(model, modelClock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep
    
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
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    !STEVE: Assuming the model goes here:
    !--------------------------------------------------------------------------
    call ESMF_ClockGet(clock, startTime=startTime, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    !STEVE: To get time information from the ESMF_Clock:
    !STEVE: http://www.earthsystemmodeling.org/esmf_releases/non_public/ESMF_1_0_8/ESMF_refdoc/node5.html#SECTION050441000000000000000

    call ESMF_TimeGet(currTime, s_i8=seconds, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    t = real(seconds)
    call ESMF_TimeIntervalGet(timeStep, s_i8=seconds, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
    dt = real(seconds)
    Nt = 1 !STEVE: just run one step of dt

    !STEVE: I'm assuming all I need to do is update the data array referenced by the pointer that is registered with the state object
!   print *, "ModelAdvance:: Pre- maooam model run:  farrayP = "
!   print *, farrayP            ! print PET-local farrayA directly
!   allocate(X(36))
!   X = farrayP(:,1)
    call maooam_atmos_run(X=farrayP,t=t,dt=dt,Nt=Nt) !,component)
!   farrayP(:,1) = X
!   print *, "ModelAdvance:: Post- maooam model run: farrayP = "
!   print *, farrayP            ! print PET-local farrayA directly


#ifdef NUOPC_TRACE
    call ESMF_TraceRegionExit("ATM:ModelAdvance")
#endif
    
  end subroutine ModelAdvance

end module ATM
