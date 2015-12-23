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
SUBROUTINE wfreq_setup
  !-----------------------------------------------------------------------
  !
  USE westcom,                ONLY : alphapv_dfpt,npwq0,npwq0_g,q0ig_l2g,sqvc,west_prefix,wfreq_dirname,&
                                   & n_pdep_eigen_to_use,n_imfreq,nbnd_occ,l_macropol,macropol_calculation,io_comm,&
                                   & n_refreq,isz,qp_bandrange
  USE mp,                     ONLY : mp_max
  USE mp_global,              ONLY : intra_bgrp_comm
  USE pwcom,                  ONLY : npw,nbnd
  USE kinds,                  ONLY : DP
  USE gvect,                  ONLY : gstart,g,ig_l2g
  USE constants,              ONLY : e2,fpi
  USE cell_base,              ONLY : tpiba2
  USE io_files,               ONLY : tmp_dir
  USE distribution_center,    ONLY : index_distr_init,pert,macropert,ifr,rfr,aband
  USE wavefunctions_module,   ONLY : evc
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
  IF(qp_bandrange(1)>nbnd) CALL errore('wfreq_setup','Err: qp_bandrange(1)>nbnd', 1) 
  IF(qp_bandrange(2)>nbnd) CALL errore('wfreq_setup','Err: qp_bandrange(2)>nbnd', 1) 
  !
  CALL set_nbndocc()
  !
  wfreq_dirname = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.wfreq.save'
  CALL my_mkdir( wfreq_dirname )
  !
  CALL index_distr_init(pert,n_pdep_eigen_to_use,'i','npdep',.TRUE.)
  CALL index_distr_init(macropert,n_pdep_eigen_to_use+3,'i','npdep+macro',.TRUE.)
  CALL index_distr_init(ifr,n_imfreq,'z','n_imfreq',.TRUE.)
  CALL index_distr_init(rfr,n_refreq,'z','n_refreq',.TRUE.)
  CALL index_distr_init(aband,nbnd,'i','nbnd',.TRUE.)
  !
  CALL set_freqlists( )
  !
  SELECT CASE(macropol_calculation)
  CASE('c','C')
     l_macropol = .TRUE.
  END SELECT
  !
  CALL define_the_iocomm( io_comm ) ! this defines communicator between heads of each image (me_bgrp==0) 
  !
END SUBROUTINE 
