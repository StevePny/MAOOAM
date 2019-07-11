
! stat.f90
!
!>  Statistics accumulators
!
!> @copyright                                                               
!> 2015 Lesley De Cruz & Jonathan Demaeyer.
!> See LICENSE.txt for license information.                                  
!
!---------------------------------------------------------------------------!



MODULE atmos_stat
  USE atmos_params, only: ndim

  IMPLICIT NONE

! PRIVATE
  PUBLIC !STEVE
  
  INTEGER :: stat_i=0 !< Number of stats accumulated  !STEVE: changed "i" to "stat_i" after making this public
                      ! (change to public is needed to implement 'model_finalize' to deallocate these arrays
  
  ! Vectors holding the stats
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: m       !< Vector storing the inline mean
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: mprev   !< Previous mean vector
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: v       !< Vector storing the inline variance
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: mtmp  


  PUBLIC :: acc,init_stat,mean,var,iter,reset

  CONTAINS

    !> Initialise the accumulators
    SUBROUTINE init_stat
      INTEGER :: AllocStat
      
      if (.not. allocated(m) .or. .not. allocated(mprev) .or. .not. allocated(v) .or. .not. allocated(mtmp)) then !STEVE: fixing init clash
        ALLOCATE(m(0:ndim),mprev(0:ndim),v(0:ndim),mtmp(0:ndim), STAT=AllocStat)
      endif
      IF (AllocStat /= 0) then
        print *, "init_stat:: AllocStat = ", AllocStat
        STOP '*** Not enough memory ***'
      ENDIF
      m=0.D0
      mprev=0.D0
      v=0.D0
      mtmp=0.D0
      
    END SUBROUTINE init_stat

    !> Accumulate one state
    SUBROUTINE acc(x)
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(0:ndim), INTENT(IN) :: x
      stat_i=stat_i+1
      mprev=m+(x-m)/stat_i
      mtmp=mprev
      mprev=m
      m=mtmp
      v=v+(x-mprev)*(x-m)
    END SUBROUTINE acc

    !> Function returning the mean
    FUNCTION mean()
      REAL(KIND=8), DIMENSION(0:ndim) :: mean
      mean=m
    END FUNCTION mean

    !> Function returning the variance
    FUNCTION var()
      REAL(KIND=8), DIMENSION(0:ndim) :: var
      var=v/(stat_i-1)
    END FUNCTION var

    !> Function returning the number of data accumulated
    FUNCTION iter()
      INTEGER :: iter
      iter=stat_i
    END FUNCTION iter

    !> Routine resetting the accumulators
    SUBROUTINE reset
      m=0.D0
      mprev=0.D0
      v=0.D0
      stat_i=0
    END SUBROUTINE reset
      

  END MODULE atmos_stat
