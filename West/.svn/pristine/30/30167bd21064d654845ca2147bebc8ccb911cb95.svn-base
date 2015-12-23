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
SUBROUTINE wstat_setup
  !-----------------------------------------------------------------------
  !
  USE westcom,                ONLY : alphapv_dfpt,npwq0,npwq0_g,q0ig_l2g,sqvc,wstat_dirname,west_prefix,&
                                   & n_pdep_basis,n_pdep_eigen,n_pdep_times,isz
  USE mp,                     ONLY : mp_max
  USE mp_global,              ONLY : intra_bgrp_comm
  USE pwcom,                  ONLY : npw
  USE kinds,                  ONLY : DP
  USE gvect,                  ONLY : gstart,g,ig_l2g
  USE constants,              ONLY : e2,fpi
  USE cell_base,              ONLY : tpiba2
  USE io_files,               ONLY : tmp_dir
  !
  IMPLICIT NONE
  !
  REAL(DP) :: q(3)
  REAL(DP) :: qq
  COMPLEX(DP),EXTERNAL :: get_alpha_pv
  INTEGER :: ig
  !
  CALL do_setup ( ) 
  !
  ! Calculate ALPHA_PV
  !
  alphapv_dfpt = get_alpha_pv()
  npwq0 = npw
  ALLOCATE(q0ig_l2g(npwq0))
  q0ig_l2g(1:npwq0) = ig_l2g(1:npwq0)
  npwq0_g=MAXVAL( q0ig_l2g(1:npwq0))
  CALL mp_max(npwq0_g,intra_bgrp_comm)
  !
  ALLOCATE(sqvc(npwq0))
  !
  CALL store_sqvc(sqvc,npwq0,1,isz)
  !
  CALL set_nbndocc()
  !
  wstat_dirname = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.wstat.save'
  CALL my_mkdir( wstat_dirname )
  !
  n_pdep_basis = n_pdep_eigen * n_pdep_times
  !
END SUBROUTINE 
