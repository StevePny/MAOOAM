program ESMF_ArrayFarrayEx
#include "ESMF.h"

  use ESMF
  use ESMF_TestMod
  
  implicit none

  ! local variables
  real(ESMF_KIND_R8)       :: farrayE(10,10)  ! explicit shape Fortran array
  real(ESMF_KIND_R8), allocatable :: farrayA(:,:) ! allocatable Fortran array
  real(ESMF_KIND_R8), pointer :: farrayP(:,:)   ! Fortran array pointer

  real(ESMF_KIND_R8), pointer :: farrayPtr(:,:) ! matching Fortran array ptr 
  type(ESMF_DistGrid)         :: distgrid       ! DistGrid object
  type(ESMF_Array)            :: array          ! Array object
  integer                     :: rc

  call ESMF_Initialize(defaultlogfilename="ArrayFarrayEx.Log", logkindflag=ESMF_LOGKIND_MULTI, rc=rc)
  if (rc /= ESMF_SUCCESS) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  farrayE = 12.45d0 ! initialize to some value

  ! Set up grid object for arrays
  distgrid = ESMF_DistGridCreate(minIndex=(/1,1/), maxIndex=(/40,10/), rc=rc)

  ! Explicit version:
  array = ESMF_ArrayCreate(farray=farrayE, distgrid=distgrid, indexflag=ESMF_INDEX_DELOCAL, rc=rc)
  print *, farrayE
  call ESMF_ArrayGet(array, farrayPtr=farrayPtr, rc=rc)
  print *, farrayPtr
  call ESMF_ArrayDestroy(array, rc=rc)

  ! Allocatable version:
  allocate(farrayA(10,10))    ! user controlled allocation
  farrayA = 23.67d0           ! initialize to some value
  array = ESMF_ArrayCreate(farray=farrayA, distgrid=distgrid, &
    indexflag=ESMF_INDEX_DELOCAL, rc=rc)
  print *, farrayA            ! print PET-local farrayA directly
  call ESMF_ArrayGet(array, farrayPtr=farrayPtr, rc=rc)! obtain array pointer
  print *, farrayPtr          ! print PET-local piece of Array through pointer
  call ESMF_ArrayDestroy(array, rc=rc) ! destroy the Array
  deallocate(farrayA)         ! user controlled de-allocation

  ! Pointer version:
  allocate(farrayP(10,10))    ! user controlled allocation
  farrayP = 56.81d0           ! initialize to some value
  array = ESMF_ArrayCreate(farray=farrayP, distgrid=distgrid, &
    indexflag=ESMF_INDEX_DELOCAL, rc=rc)
  print *, farrayP            ! print PET-local farrayA directly
  call ESMF_ArrayGet(array, farrayPtr=farrayPtr, rc=rc)! obtain array pointer
  print *, farrayPtr          ! print PET-local piece of Array through pointer
  call ESMF_ArrayDestroy(array, rc=rc) ! destroy the Array
  deallocate(farrayP)         ! user controlled de-allocation

  ! wrap up
  call ESMF_DistGridDestroy(distgrid, rc=rc) ! destroy the DistGrid
  call ESMF_Finalize(rc=rc)
end program
