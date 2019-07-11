module maooam_atmos_wrapper
!STEVE:
! Wrapper module to (i)nitialize, (r)un, and (f)inalize
! the MAOOAM model


public :: maooam_atmos_initialize, maooam_atmos_run, maooam_atmos_finalize
public :: maooam_natm
public :: maooam_nocn

integer :: maooam_natm
integer :: maooam_nocn

private

logical, save :: first_flag = .true.
character(3), parameter :: component = 'atm'

contains


!------------------------------------------------------------------------------
subroutine maooam_atmos_initialize
!------------------------------------------------------------------------------

use atmos_aotensor_def, only: init_aotensor
use atmos_integrator, only: init_integrator
use atmos_stat, only: init_stat  ! Running statistics data

use atmos_params, only: natm
use atmos_params, only: noc

!real(kind=8), dimension(:), pointer, intent(inout) :: X0

print *, "maooam_atmos_wrapper :: call init_aotensor... (Reads the namelist files: params.nml, modeselection.nml, int_params.nml)"
call init_aotensor
print *, "maooam_atmos_wrapper :: call init_integrator..."
call init_integrator  ! Initialize the integrator
print *, "maooam_atmos_wrapper :: call init_stat..."
call init_stat        ! Initialize the statistics operations

! Initialize needed parameters for ESMF cap
maooam_natm = natm
maooam_nocn = noc

end subroutine maooam_atmos_initialize



!------------------------------------------------------------------------------
subroutine maooam_atmos_run(X,t,dt,Nt,component)
!------------------------------------------------------------------------------

use atmos_integrator, only: step
use atmos_ic_def, only: load_IC

use atmos_ic_def, only: IC

real(kind=8), dimension(:), pointer, intent(inout)  :: X
real(kind=8), intent(inout) :: t
real(kind=8), intent(in) :: dt
integer, intent(in) :: Nt
character(3), intent(in), optional :: component !empty, 'atm', or 'ocn'. If the latter two cases, uses the other part as fixed forcing
! https://stackoverflow.com/questions/43351338/passing-a-value-for-an-optional-fortran-parameter-that-will-return-false-for-pre
!character(3), allocatable, intent(in), optional :: component !empty, 'atm', or 'ocn'. If the latter two cases, uses the other part as fixed forcing
real(kind=8), allocatable, dimension(:) :: Xnew
integer :: n

! Load the initial conditions
if (first_flag) then
  call load_IC
  X = IC
  first_flag=.false.
endif

! Setup array for iterating
n = size(X)
allocate(Xnew(n))

print *, "maooam_atmos_run :: init X = "
print *, X

!Cycle MAOOAM through this time step of JEDI
do n = 1,Nt
! call step(X,t,dt,Xnew)
  ! Use the following to support forcing from atmosphere or ocean on the full state.
  call step(X,t,dt,Xnew,component)
  X=Xnew
  ! Debug:
  print *, X
enddo
print *, "maooam_atmos_run:: final X = "
print *, X

end subroutine maooam_atmos_run



!------------------------------------------------------------------------------
subroutine maooam_atmos_finalize()
!------------------------------------------------------------------------------

use atmos_params, only: ams, oms    
use atmos_inprod_analytic, only: awavenum, owavenum
use atmos_aotensor_def, only: aotensor !, count_elems
use atmos_integrator, only: buf_y1,buf_f0,buf_f1
use atmos_stat, only: m,mprev,v,mtmp,stat_i

print *, "maooam_atmos_finalize :: start..."

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

print *, "maooam_atmos_finalize :: finished."

end subroutine maooam_atmos_finalize


end module maooam_atmos_wrapper
