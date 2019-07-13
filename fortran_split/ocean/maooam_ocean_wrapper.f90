module maooam_ocean_wrapper
!STEVE:
! Wrapper module to (i)nitialize, (r)un, and (f)inalize
! the MAOOAM model


public :: maooam_ocean_initialize, maooam_ocean_run, maooam_ocean_finalize
public :: maooam_nocn
public :: maooam_natm

integer :: maooam_nocn
integer :: maooam_natm

private

logical, save :: first_flag = .true.
character(3), parameter :: component = 'ocn'

contains


!------------------------------------------------------------------------------
subroutine maooam_ocean_initialize
!------------------------------------------------------------------------------

use ocean_aotensor_def, only: init_aotensor
use ocean_integrator, only: init_integrator
use ocean_stat, only: init_stat  ! Running statistics data

use ocean_params, only: noc
use ocean_params, only: natm
use ocean_aotensor_def, only: aotensor

!real(kind=8), dimension(:), pointer, intent(inout) :: X0

if (allocated(aotensor)) then
  print *, "maooam_ocean_wrapper::maooam_ocean_initialize:: aotensor already initialized. skipping initialization..."
else
  print *, "maooam_ocean_wrapper :: call init_aotensor... (Reads the namelist files: params.nml, modeselection.nml, int_params.nml)"
  call init_aotensor
  print *, "maooam_ocean_wrapper :: call init_integrator..."
  call init_integrator  ! Initialize the integrator
  print *, "maooam_ocean_wrapper :: call init_stat..."
  call init_stat        ! Initialize the statistics operations
endif

! Initialize needed parameters for ESMF cap
maooam_nocn = noc
maooam_natm = natm

end subroutine maooam_ocean_initialize



!------------------------------------------------------------------------------
subroutine maooam_ocean_run(X,t,dt,Nt)
!------------------------------------------------------------------------------

use ocean_integrator, only: step
use ocean_ic_def, only: load_IC

use ocean_ic_def, only: IC
use ocean_params, only: ndim

real(kind=8), dimension(:), pointer, intent(inout)  :: X
real(kind=8), intent(inout) :: t
real(kind=8), intent(in) :: dt
integer, intent(in) :: Nt
real(kind=8), allocatable, dimension(:) :: Xnew
integer :: n

logical :: local_verbose = .true.

! Load the initial conditions
if (first_flag) then
  call load_IC

  if (local_verbose) then
    print *, "ocean::maooam_run :: X = "
    print *, X
    print *, "ocean::maooam_run :: IC = "
    print *, IC
  endif

  X = IC
  first_flag=.false.
endif

! Setup array for iterating
n = size(X)
if (ndim+1 .ne. n) then
  print *, "ocean::maooam_run:: ERROR"
  print *, "ndim+1 = ", ndim+1
  print *, "size(X) = ", n
  stop 1
endif
allocate(Xnew(0:ndim))

if (local_verbose) then
  print *, "====================================================================="
  print *, "ocean::maooam_run :: init X = "
  print *, X
  print *, "---------------------------------------------------------------------"
endif

!Cycle MAOOAM through this time step of JEDI
do i = 1,Nt
  if (local_verbose) then
    print *, "i, t, dt = ", i, t, dt
    print *, "X0 = "
    print *, X
  endif
  call step(X,t,dt,Xnew,component)
  X=Xnew
  ! Debug:
  if (local_verbose) then
    print *, "Xnew = "
    print *, Xnew
    print *, "-------------------------------------------------------------------"
  endif
enddo

if (local_verbose) then
  print *, "---------------------------------------------------------------------"
  print *, "ocean::maooam_run:: final X = "
  print *, X
  print *, "====================================================================="
endif

end subroutine maooam_ocean_run



!------------------------------------------------------------------------------
subroutine maooam_ocean_finalize()
!------------------------------------------------------------------------------

use ocean_params, only: ams, oms    
use ocean_inprod_analytic, only: awavenum, owavenum
use ocean_aotensor_def, only: aotensor !, count_elems
use ocean_integrator, only: buf_y1,buf_f0,buf_f1
use ocean_stat, only: m,mprev,v,mtmp,stat_i

print *, "maooam_ocean_finalize :: start..."

!STEVE: these data structures need to be manually
!       removed because they are initiated inside
!       the MAOOAM model and cause errors if
!       the model is re-initiated externally

! From maooam file params.f90:
if (allocated(ams)) then
  deallocate(ams)
endif
if (allocated(oms)) then
  deallocate(oms)
endif
! From maooam file inprod_analytic.f90:
if (allocated(awavenum)) then
  deallocate(awavenum)
endif
if (allocated(owavenum)) then
  deallocate(owavenum)
endif
! From maooam file aotensor_def.f90:
if (allocated(aotensor)) then
  deallocate(aotensor) !,count_elems)
endif
! From maooam file rkX.integrator.f90:
if (allocated(buf_y1)) then
  deallocate(buf_y1)
endif
if (allocated(buf_f0)) then
  deallocate(buf_f0)
endif
if (allocated(buf_f1)) then
  deallocate(buf_f1)
endif
! From maooamfile stat.f90:
stat_i = 0
if (allocated(m)) then
  deallocate(m)
endif
if (allocated(mprev)) then
  deallocate(mprev)
endif
if (allocated(v)) then
  deallocate(v)
endif
if (allocated(mtmp)) then
  deallocate(mtmp)
endif

print *, "maooam_ocean_finalize :: finished."

end subroutine maooam_ocean_finalize


end module maooam_ocean_wrapper
