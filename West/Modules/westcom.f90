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
MODULE scratch_area
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY :  DP
  !
  SAVE
  !
  ! Coulomb
  REAL(DP),ALLOCATABLE :: sqvc(:)
  INTEGER :: npwq0, npwq0_g
  INTEGER,ALLOCATABLE :: q0ig_l2g(:)
  INTEGER,ALLOCATABLE :: iks_l2g(:)
  !
  ! DBS
  REAL(DP),ALLOCATABLE :: ev(:)
  REAL(DP),ALLOCATABLE :: ev_distr(:)
  COMPLEX(DP),ALLOCATABLE :: dng(:,:)
  COMPLEX(DP),ALLOCATABLE :: dvg(:,:)
  LOGICAL,ALLOCATABLE :: conv(:)
  !
  INTEGER,ALLOCATABLE :: nbnd_occ(:) 
  !
  ! EPSILON
  REAL(DP),ALLOCATABLE :: dmat_inv(:,:,:)
  COMPLEX(DP),ALLOCATABLE :: zmat_inv(:,:,:)
  !
  ! CORR
  REAL(DP),ALLOCATABLE :: ihead(:)
  REAL(DP),ALLOCATABLE :: ibody1(:,:,:,:)
  REAL(DP),ALLOCATABLE :: ibody2(:,:,:,:,:)
  REAL(DP),ALLOCATABLE :: idiago(:,:,:,:)
  COMPLEX(DP),ALLOCATABLE :: rhead(:)
  COMPLEX(DP),ALLOCATABLE :: rbody(:,:,:,:) 
  !
  ! I/O 
  INTEGER :: io_comm ! communicator for head of images (me_bgrp==0)
  !
  REAL(DP) :: isz
  !
END MODULE 
!
!
MODULE wstat_center
  !
  USE kinds, ONLY :  DP
  !
  SAVE
  !
  ! INPUT FOR wstat_control
  !
  CHARACTER(LEN=1) :: wstat_calculation 
  INTEGER :: n_pdep_basis
  INTEGER :: n_pdep_times
  INTEGER :: n_pdep_eigen
  INTEGER :: n_pdep_maxiter
  INTEGER :: n_dfpt_maxiter
  INTEGER :: n_steps_write_restart
  INTEGER :: n_pdep_restart_from_itr
  INTEGER :: n_pdep_read_from_file
  REAL(DP) :: trev_pdep
  REAL(DP) :: tr2_dfpt
  LOGICAL :: l_deflate
  LOGICAL :: l_kinetic_only
  !
  ! Common workspace
  !
  CHARACTER(LEN=256) :: west_prefix
  COMPLEX(DP) :: alphapv_dfpt
  CHARACTER(LEN=256) :: wstat_dirname
  LOGICAL :: l_minimize_exx_if_active
  !
END MODULE 
!
!
MODULE wfreq_center
  !
  USE kinds, ONLY : DP
  !
  SAVE
  !
  ! INPUT FOR wfreq_control
  !
  CHARACTER(LEN=8) :: wfreq_calculation
  INTEGER :: n_lanczos
  INTEGER :: n_pdep_eigen_to_use
  INTEGER :: n_imfreq
  INTEGER :: n_refreq
  INTEGER :: qp_bandrange(2)
  REAL(DP) :: ecut_imfreq
  REAL(DP) :: ecut_refreq
  REAL(DP) :: wfreq_eta
  INTEGER :: n_secant_maxiter
  REAL(DP) :: trev_secant
  LOGICAL :: l_enable_lanczos
  CHARACTER(LEN=1) :: macropol_calculation
  LOGICAL :: l_enable_gwetot
  REAL(DP) :: exx_etot
  REAL(DP) :: o_restart_time
  !
  ! Common workspace
  !
  CHARACTER(LEN=256) :: wfreq_dirname
  LOGICAL,PARAMETER :: l_skip_nl_part_of_hcomr=.FALSE.
  LOGICAL :: l_macropol
  !
  ! re freq 
  !
  REAL(DP),ALLOCATABLE :: refreq_list(:)
  !
  ! im freq
  !
  REAL(DP),ALLOCATABLE :: imfreq_list(:)
  REAL(DP),ALLOCATABLE :: imfreq_list_integrate(:,:)
  REAL(DP),PARAMETER :: frequency_list_power = 2._DP
  INTEGER :: div_kind_hf ! 1=spherical region, 2=GB, 3=cut_ws
  !
  ! gw_etot 
  !
  REAL(DP) :: dft_etot
  REAL(DP) :: dft_exc
  REAL(DP) :: gw_ecorr
  REAL(DP) :: gw_exx
  REAL(DP) :: gw_exc
  ! 
  !
END MODULE
!
!
MODULE wan_center
  !
  USE kinds, ONLY : DP
  !
  SAVE
  !
  REAL(DP),ALLOCATABLE :: wanc(:,:)
  REAL(DP),ALLOCATABLE :: wanu(:,:)
  INTEGER :: wantot
  !
END MODULE
!
!
MODULE io_unit_numbers
  !
  SAVE
  !
  INTEGER,PARAMETER :: iuwfc=20
  INTEGER :: lrwfc 
  !
END MODULE 
!
!
MODULE westcom
  !
  USE scratch_area
  USE wstat_center
  USE wfreq_center
  USE wan_center
  USE io_unit_numbers
  !
END MODULE
