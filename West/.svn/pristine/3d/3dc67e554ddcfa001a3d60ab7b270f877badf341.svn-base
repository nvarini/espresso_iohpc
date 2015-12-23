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
SUBROUTINE exx_ungo()
  !-----------------------------------------------------------------------
  !
  USE funct,                  ONLY : dft_is_hybrid,init_dft_exxrpa,stop_exx
  USE exx,                    ONLY : x_gamma_extrapolation,exxdiv_treatment,exx_grid_init,exx_div_check,&
                                     &deallocate_exx,exxinit,vexx
  !
  IMPLICIT NONE
  !
  IF( dft_is_hybrid() ) THEN
     CALL stop_exx
     CALL deallocate_exx
  ENDIF
  !
END SUBROUTINE 
