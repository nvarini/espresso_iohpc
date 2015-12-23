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
PROGRAM wfreq
  !-----------------------------------------------------------------------
  ! 
  ! This is the main program that calculates the GW.
  !
  USE check_stop,           ONLY : check_stop_init
  USE mp_global,            ONLY : mp_startup, mp_global_end
  USE west_environment,     ONLY : west_environment_start, west_environment_end 
  USE mp,                   ONLY : mp_sum,mp_barrier
  USE control_flags,        ONLY : gamma_only
  USE westcom,              ONLY : wfreq_calculation
  ! 
  IMPLICIT NONE
  !
  CHARACTER(LEN=9) :: code = 'WFREQ'
  LOGICAL :: lgate(6)
  INTEGER :: i
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
  CALL wfreq_readin( )
  !
  CALL wfreq_setup( )
  !
  lgate=.FALSE.
  DO i = 1, 8 
     IF( wfreq_calculation(i:i) == 'X' ) lgate(1) = .TRUE.
     IF( wfreq_calculation(i:i) == 'W' ) lgate(2) = .TRUE.
     IF( wfreq_calculation(i:i) == 'w' ) lgate(3) = .TRUE.
     IF( wfreq_calculation(i:i) == 'G' ) lgate(4) = .TRUE.
     IF( wfreq_calculation(i:i) == 'g' ) lgate(5) = .TRUE.
     IF( wfreq_calculation(i:i) == 'Q' ) lgate(6) = .TRUE.
  ENDDO
  !
  IF(lgate(1)) THEN
     CALL solve_hf( )
  ENDIF
  !
  IF(lgate(2)) THEN
     CALL solve_wfreq(.FALSE.)
  ENDIF
  !
  IF(lgate(3)) THEN
     CALL solve_wfreq(.TRUE.)
  ENDIF
  !
  IF(lgate(4)) THEN
     CALL solve_gfreq(.FALSE.)
  ENDIF
  !
  IF(lgate(5)) THEN
     CALL solve_gfreq(.TRUE.)
  ENDIF
  !
  IF(lgate(6)) THEN
     CALL solve_qp( )
  ENDIF
  !
  CALL clean_scratchfiles( )
  !
  CALL print_clock(' ')
  !
  CALL west_environment_end( code )
  !
  CALL mp_global_end()
  !
  STOP
  !
9000 FORMAT (/5x,'Program ',a12,' starts ...',/5x,'Today is ',a9,' at ',a9)
  !
END PROGRAM
