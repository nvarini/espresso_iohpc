!
! Copyright (C) 2015 M. Govoni 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors to this file: 
! Marco Govoni
!
!-----------------------------------------------------------------------
PROGRAM wstat
  !-----------------------------------------------------------------------
  ! 
  ! This is the main program that calculates the static screening.
  !
  USE check_stop,           ONLY : check_stop_init
  USE mp_global,            ONLY : mp_startup, mp_global_end
  USE west_environment,     ONLY : west_environment_start, west_environment_end
  USE mp,                   ONLY : mp_sum,mp_barrier
  USE control_flags,        ONLY : gamma_only
  ! 
  IMPLICIT NONE
  !
  CHARACTER(LEN=9) :: code = 'WSTAT'
  !
  ! *** START *** 
  !
  CALL check_stop_init ()
  !
  ! Initialize MPI, clocks, print initial messages
  !
#ifdef __MPI
  CALL mp_startup ( start_images = .TRUE. )
#endif
  !
  CALL west_environment_start ( code )
  !
  CALL wstat_readin ( )
  !
  CALL wstat_setup ( )
  !
  CALL davidson_diago ( )
  !
  CALL exx_ungo ( )
  !
  CALL clean_scratchfiles( )
  !
  CALL print_clock(' ')
  !
  CALL west_environment_end( code )
  !
  CALL mp_global_end()
  !
9000 FORMAT (/5x,'Program ',a12,' starts ...',/5x,'Today is ',a9,' at ',a9)
  !
END PROGRAM
