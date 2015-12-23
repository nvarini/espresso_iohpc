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
SUBROUTINE preallocate_solvewfreq( ik, nb_1, nb_2, mypara )
   !-----------------------------------------------------------------------
   !
   USE kinds,   ONLY : DP
   USE distribution_center, ONLY : index_distribution
   USE mp_images,   ONLY : nimage,my_image_id,inter_image_comm
   USE mp,          ONLY : mp_sum
   USE westcom,     ONLY : n_pdep_eigen_to_use,iks_l2g,n_lanczos,wfreq_dirname
   !
   IMPLICIT NONE
   !
   ! I/O
   !
   INTEGER,INTENT(IN) :: ik, nb_1, nb_2
   TYPE(index_distribution), INTENT(IN) :: mypara
   !
   ! Workspace
   !
   INTEGER :: ib, ii
   CHARACTER(LEN=5) :: c_glob_iks
   CHARACTER(LEN=5) :: c_glob_ib
   CHARACTER(LEN=256) :: fname
   !
   WRITE(c_glob_iks,'(i5.5)') ik
   DO ib = nb_1, nb_2
      !
      WRITE(c_glob_ib, '(i5.5)') ib
      !
      fname = TRIM( wfreq_dirname )//"/w_diag_K"//c_glob_iks//"B"//c_glob_ib//".dat"
      CALL mp_master_creates_and_preallocates( fname, n_lanczos*mypara%nglob )   
      fname = TRIM( wfreq_dirname )//"/w_brak_K"//c_glob_iks//"B"//c_glob_ib//".dat"
      CALL mp_master_creates_and_preallocates( fname, n_lanczos*mypara%nglob*mypara%nglob )   
      !
   ENDDO
   !
END SUBROUTINE
!
!
!
SUBROUTINE preallocate_solvegfreq( ik, nb_1, nb_2, mypara )
   !
   USE kinds,   ONLY : DP
   USE distribution_center, ONLY : index_distribution
   USE mp_images,   ONLY : nimage,my_image_id,inter_image_comm
   USE mp,          ONLY : mp_sum
   USE westcom,     ONLY : n_pdep_eigen_to_use,iks_l2g,n_lanczos,wfreq_dirname
   !
   IMPLICIT NONE
   !
   ! I/O
   !
   INTEGER,INTENT(IN) :: ik, nb_1, nb_2
   TYPE(index_distribution), INTENT(IN) :: mypara
   !
   ! Workspace
   !
   INTEGER :: ib, ii
   CHARACTER(LEN=5) :: c_glob_iks
   CHARACTER(LEN=5) :: c_glob_ib
   CHARACTER(LEN=256) :: fname
   !
   WRITE(c_glob_iks,'(i5.5)') ik
   DO ib = nb_1, nb_2
      !
      WRITE(c_glob_ib, '(i5.5)') ib
      !
      fname = TRIM( wfreq_dirname )//"/g_diag_K"//c_glob_iks//"B"//c_glob_ib//".dat"
      CALL mp_master_creates_and_preallocates( fname, n_lanczos*mypara%nglob )   
      fname = TRIM( wfreq_dirname )//"/g_brak_K"//c_glob_iks//"B"//c_glob_ib//".dat"
      CALL mp_master_creates_and_preallocates( fname, n_lanczos*mypara%nglob*mypara%nglob )   
      !
   ENDDO
   !
END SUBROUTINE
!
!
!
SUBROUTINE writeout_solvewfreq( glob_iks, glob_ib, diago, braket, io_comm, nloc, nglob, myoffset)
  !
  ! WHO WRITES? ... root_bgrp
  ! WHAT?       ... diago(        n_lanczos ) 
  ! WHAT?       ... braket( ndim, n_lanczos ) 
  !
  USE kinds,          ONLY : DP
  USE westcom,        ONLY : n_pdep_eigen_to_use,iks_l2g,n_lanczos,wfreq_dirname
  USE mp_global,      ONLY : my_image_id,intra_bgrp_comm,me_bgrp,root_bgrp
  USE mp,             ONLY : mp_bcast
  USE distribution_center, ONLY : index_distribution
  !
  IMPLICIT NONE
  !
  ! I/O 
  !
  INTEGER,INTENT(IN) :: glob_iks
  INTEGER,INTENT(IN) :: glob_ib
  REAL(DP),INTENT(IN) :: diago(n_lanczos,nloc)
  REAL(DP),INTENT(IN) :: braket(nglob,n_lanczos,nloc)
  INTEGER,INTENT(IN) :: io_comm
  INTEGER,INTENT(IN) :: nloc, nglob, myoffset 
  !
  ! Workspace
  !
  CHARACTER(LEN=5) :: c_glob_iks
  CHARACTER(LEN=5) :: c_glob_ib
  CHARACTER(LEN=256) :: fname
  !
  ! Generate the filename
  ! 
  WRITE(c_glob_iks,'(i5.5)') glob_iks
  WRITE(c_glob_ib, '(i5.5)') glob_ib
  !
  fname = TRIM( wfreq_dirname )//"/w_diag_K"//c_glob_iks//"B"//c_glob_ib//".dat"
  CALL mp_write_d2_at( fname, io_comm, diago, n_lanczos, nloc, myoffset ) 
  !
  fname = TRIM( wfreq_dirname )//"/w_brak_K"//c_glob_iks//"B"//c_glob_ib//".dat"
  CALL mp_write_d3_at( fname, io_comm, braket, nglob, n_lanczos, nloc, myoffset ) 
  !
END SUBROUTINE
!
!
!
SUBROUTINE readin_solvewfreq( glob_iks, glob_ib, diago, braket, io_comm, nloc, nglob, myoffset )
  !
  ! WHO READS?  ... root_bgrp
  ! WHAT?       ... diago(        n_lanczos ) 
  ! WHAT?       ... braket( ndim, n_lanczos ) 
  !
  USE kinds,          ONLY : DP
  USE westcom,        ONLY : n_pdep_eigen_to_use,iks_l2g,n_lanczos,wfreq_dirname
  USE mp_global,      ONLY : my_image_id,intra_bgrp_comm,me_bgrp,root_bgrp
  USE mp,             ONLY : mp_bcast
  USE iotk_module
  !
  IMPLICIT NONE
  !
  ! I/O 
  !
  INTEGER,INTENT(IN) :: glob_iks
  INTEGER,INTENT(IN) :: glob_ib
  REAL(DP),INTENT(OUT) :: diago(n_lanczos,nloc)
  REAL(DP),INTENT(OUT) :: braket(nglob,n_lanczos,nloc)
  INTEGER,INTENT(IN) :: io_comm
  INTEGER,INTENT(IN) :: nloc,nglob,myoffset
  !
  ! Workspace
  !
  CHARACTER(LEN=5) :: c_glob_iks
  CHARACTER(LEN=5) :: c_glob_ib
  CHARACTER(LEN=256) :: fname
  !
  ! Generate the filename
  ! 
  WRITE(c_glob_iks,'(i5.5)') glob_iks
  WRITE(c_glob_ib, '(i5.5)') glob_ib
  !
  fname = TRIM( wfreq_dirname )//"/w_diag_K"//c_glob_iks//"B"//c_glob_ib//".dat"
  CALL mp_read_d2_at( fname, io_comm, diago, n_lanczos, nloc, myoffset )
  !
  fname = TRIM( wfreq_dirname )//"/w_brak_K"//c_glob_iks//"B"//c_glob_ib//".dat"
  CALL mp_read_d3_at( fname, io_comm, braket, nglob, n_lanczos, nloc, myoffset )
  !
END SUBROUTINE
!
!
SUBROUTINE writeout_imfreq( glob_ifreq,head,lambda )
  !
  ! WHO WRITES? ... me_bgrp ( IF my_image_id == 0 ) 
  ! WHAT?       ... head 
  ! WHAT?       ... lambda( n_pdep_eigen_to_use, n_pdep_eigen_to_use ) 
  !
  USE kinds,          ONLY : DP
  USE westcom,        ONLY : n_pdep_eigen_to_use,wfreq_dirname
  USE mp_global,      ONLY : my_image_id,inter_image_comm,me_image,root_image
  USE mp,             ONLY : mp_bcast
  USE iotk_module
  !
  IMPLICIT NONE
  !
  ! I/O 
  !
  INTEGER,INTENT(IN) :: glob_ifreq
  REAL(DP),INTENT(IN) :: head
  REAL(DP),INTENT(IN) :: lambda(n_pdep_eigen_to_use,n_pdep_eigen_to_use)
  !
  ! Workspace
  !
  CHARACTER(LEN=6) :: c_glob_ifreq
  INTEGER :: iunout, ierr
  CHARACTER(LEN=256) :: fname
  ! 
  !
  WRITE(c_glob_ifreq,'(i6.6)') glob_ifreq
  fname = TRIM( wfreq_dirname )//"/w_imF"//c_glob_ifreq//".dat"
  !
  IF ( my_image_id == 0 ) THEN
     !
     ! ... open XML descriptor
     !
     CALL iotk_free_unit( iunout, ierr )
     CALL iotk_open_write( iunout, FILE = fname , BINARY = .FALSE., IERR = ierr )
     !
  END IF
  !
  CALL mp_bcast( ierr, 0, inter_image_comm )
  CALL errore( 'freq', 'cannot open freq file for writing', ierr )
  !
  IF ( my_image_id == 0 ) THEN  
     !
     CALL iotk_write_begin( iunout, "imCHI" )
     CALL iotk_write_dat( iunout, "head", head )
     CALL iotk_write_dat( iunout, "lambda", lambda )
     CALL iotk_write_end( iunout, "imCHI" )
     !
     ! ... close XML descriptor
     !
     CALL iotk_close_write( iunout )
     !
  END IF
  !
END SUBROUTINE
!
!
SUBROUTINE readin_imfreq( glob_ifreq,head,lambda )
  !
  ! WHO READS?  ... me_bgrp .AND. my_image_id==0; bcasted via inter_image_comm
  ! WHAT?       ... head 
  ! WHAT?       ... lambda( n_pdep_eigen_to_use, n_pdep_eigen_to_use ) 
  !
  USE kinds,          ONLY : DP
  USE westcom,        ONLY : n_pdep_eigen_to_use,wfreq_dirname
  USE mp_global,      ONLY : my_image_id,inter_image_comm
  USE mp,             ONLY : mp_bcast
  USE iotk_module
  !
  IMPLICIT NONE
  !
  ! I/O 
  !
  INTEGER,INTENT(IN) :: glob_ifreq
  REAL(DP),INTENT(OUT) :: head
  REAL(DP),INTENT(OUT) :: lambda(n_pdep_eigen_to_use,n_pdep_eigen_to_use)
  !
  ! Workspace
  !
  CHARACTER(LEN=6) :: c_glob_ifreq
  INTEGER :: iunout, ierr
  CHARACTER(LEN=256) :: fname
  ! 
  WRITE(c_glob_ifreq,'(i6.6)') glob_ifreq
  fname = TRIM( wfreq_dirname )//"/w_imF"//c_glob_ifreq//".dat"
  !
  !
  IF ( my_image_id == 0 ) THEN
     !
     ! ... open XML descriptor
     !
     CALL iotk_free_unit( iunout, ierr )
     CALL iotk_open_read( iunout, FILE = fname , BINARY = .FALSE., IERR = ierr )
     !
  END IF
  !
  CALL mp_bcast( ierr, 0, inter_image_comm )
  CALL errore( 'freq', 'cannot open freq file for reading', ierr )
  !
  IF ( my_image_id == 0 ) THEN  
     !
     CALL iotk_scan_begin( iunout, "imCHI" )
     CALL iotk_scan_dat( iunout, "head", head )
     CALL iotk_scan_dat( iunout, "lambda", lambda )
     CALL iotk_scan_end( iunout, "imCHI" )
     !
     ! ... close XML descriptor
     !
     CALL iotk_close_read( iunout )
     !
  END IF
  !
  CALL mp_bcast( head  , 0, inter_image_comm)
  CALL mp_bcast( lambda, 0, inter_image_comm)
  !
END SUBROUTINE
!
!
!
SUBROUTINE writeout_solvegfreq( glob_iks, glob_ib, diago, braket, io_comm, nloc, nglob, myoffset)
  !
  ! WHO WRITES? ... root_bgrp
  ! WHAT?       ... diago(        n_lanczos ) 
  ! WHAT?       ... braket( ndim, n_lanczos ) 
  !
  USE kinds,          ONLY : DP
  USE westcom,        ONLY : n_pdep_eigen_to_use,iks_l2g,n_lanczos,wfreq_dirname
  USE mp_global,      ONLY : my_image_id,intra_bgrp_comm,me_bgrp,root_bgrp
  USE mp,             ONLY : mp_bcast
  USE distribution_center, ONLY : index_distribution
  !
  IMPLICIT NONE
  !
  ! I/O 
  !
  INTEGER,INTENT(IN) :: glob_iks
  INTEGER,INTENT(IN) :: glob_ib
  REAL(DP),INTENT(IN) :: diago(n_lanczos,nloc)
  REAL(DP),INTENT(IN) :: braket(nglob,n_lanczos,nloc)
  INTEGER,INTENT(IN) :: io_comm
  INTEGER,INTENT(IN) :: nloc, nglob, myoffset 
  !
  ! Workspace
  !
  CHARACTER(LEN=5) :: c_glob_iks
  CHARACTER(LEN=5) :: c_glob_ib
  CHARACTER(LEN=256) :: fname
  !
  ! Generate the filename
  ! 
  WRITE(c_glob_iks,'(i5.5)') glob_iks
  WRITE(c_glob_ib, '(i5.5)') glob_ib
  !
  fname = TRIM( wfreq_dirname )//"/g_diag_K"//c_glob_iks//"B"//c_glob_ib//".dat"
  CALL mp_write_d2_at( fname, io_comm, diago, n_lanczos, nloc, myoffset ) 
  !
  fname = TRIM( wfreq_dirname )//"/g_brak_K"//c_glob_iks//"B"//c_glob_ib//".dat"
  CALL mp_write_d3_at( fname, io_comm, braket, nglob, n_lanczos, nloc, myoffset ) 
  !
END SUBROUTINE
!
!
!
SUBROUTINE readin_solvegfreq( glob_iks, glob_ib, diago, braket, io_comm, nloc, nglob, myoffset )
  !
  ! WHO READS?  ... root_bgrp
  ! WHAT?       ... diago(        n_lanczos ) 
  ! WHAT?       ... braket( ndim, n_lanczos ) 
  !
  USE kinds,          ONLY : DP
  USE westcom,        ONLY : n_pdep_eigen_to_use,iks_l2g,n_lanczos,wfreq_dirname
  USE mp_global,      ONLY : my_image_id,intra_bgrp_comm,me_bgrp,root_bgrp
  USE mp,             ONLY : mp_bcast
  USE iotk_module
  !
  IMPLICIT NONE
  !
  ! I/O 
  !
  INTEGER,INTENT(IN) :: glob_iks
  INTEGER,INTENT(IN) :: glob_ib
  REAL(DP),INTENT(OUT) :: diago(n_lanczos,nloc)
  REAL(DP),INTENT(OUT) :: braket(nglob,n_lanczos,nloc)
  INTEGER,INTENT(IN) :: io_comm
  INTEGER,INTENT(IN) :: nloc,nglob,myoffset
  !
  ! Workspace
  !
  CHARACTER(LEN=5) :: c_glob_iks
  CHARACTER(LEN=5) :: c_glob_ib
  CHARACTER(LEN=256) :: fname
  !
  ! Generate the filename
  ! 
  WRITE(c_glob_iks,'(i5.5)') glob_iks
  WRITE(c_glob_ib, '(i5.5)') glob_ib
  !
  fname = TRIM( wfreq_dirname )//"/g_diag_K"//c_glob_iks//"B"//c_glob_ib//".dat"
  CALL mp_read_d2_at( fname, io_comm, diago, n_lanczos, nloc, myoffset )
  !
  fname = TRIM( wfreq_dirname )//"/g_brak_K"//c_glob_iks//"B"//c_glob_ib//".dat"
  CALL mp_read_d3_at( fname, io_comm, braket, nglob, n_lanczos, nloc, myoffset )
  !
END SUBROUTINE
!
!
!
SUBROUTINE writeout_overlap( labellina, glob_iks, glob_ib, overlap, no1, no2 )
  !
  ! WHO WRITES? ... root
  ! WHAT?       ... overlap( no1, no2 ) 
  !
  USE kinds,          ONLY : DP
  USE westcom,        ONLY : iks_l2g,wfreq_dirname
  USE mp_world,       ONLY : mpime,root,world_comm 
  USE mp,             ONLY : mp_bcast
  USE iotk_module
  !
  IMPLICIT NONE
  !
  ! I/O 
  !
  CHARACTER(LEN=1),INTENT(IN) :: labellina
  INTEGER,INTENT(IN) :: glob_iks
  INTEGER,INTENT(IN) :: glob_ib
  REAL(DP),INTENT(IN) :: overlap(no1,no2)
  INTEGER,INTENT(IN) :: no1,no2
  !
  ! Workspace
  !
  CHARACTER(LEN=5) :: c_glob_iks
  CHARACTER(LEN=5) :: c_glob_ib
  INTEGER :: iunout, ierr
  CHARACTER(LEN=256) :: fname
  !
  ! Generate the filename
  ! 
  WRITE(c_glob_iks,'(i5.5)') glob_iks
  WRITE(c_glob_ib, '(i5.5)') glob_ib
  !
  fname = TRIM( wfreq_dirname )//"/over_"//labellina//"_K"//c_glob_iks//"B"//c_glob_ib//".dat"
  !
  IF ( mpime == root ) THEN
     !
     ! ... open XML descriptor
     !
     CALL iotk_free_unit( iunout, ierr )
     CALL iotk_open_write( iunout, FILE = fname , BINARY = .TRUE., IERR = ierr )
     !
  END IF
  !
  CALL mp_bcast( ierr, root, world_comm )
  CALL errore( 'freq', 'cannot open freq file for writing', ierr )
  !
  IF ( mpime == root ) THEN  
     !
     CALL iotk_write_begin( iunout, "OVER" )
     CALL iotk_write_dat( iunout, "overlap", overlap )
     CALL iotk_write_end( iunout, "OVER" )
     !
     ! ... close XML descriptor
     !
     CALL iotk_close_write( iunout )
     !
  END IF
  !
END SUBROUTINE
!
!
!
!
SUBROUTINE readin_overlap( labellina, glob_iks, glob_ib, overlap, no1, no2 )
  !
  ! WHO WRITES? ... root
  ! WHAT?       ... overlap( no1, no2 ) 
  !
  USE kinds,          ONLY : DP
  USE westcom,        ONLY : iks_l2g,wfreq_dirname
  USE mp_world,       ONLY : mpime,root,world_comm 
  USE mp,             ONLY : mp_bcast
  USE iotk_module
  !
  IMPLICIT NONE
  !
  ! I/O 
  !
  CHARACTER(LEN=1),INTENT(IN) :: labellina
  INTEGER,INTENT(IN) :: glob_iks
  INTEGER,INTENT(IN) :: glob_ib
  REAL(DP),INTENT(OUT) :: overlap(no1,no2)
  INTEGER,INTENT(IN) :: no1,no2
  !
  ! Workspace
  !
  CHARACTER(LEN=5) :: c_glob_iks
  CHARACTER(LEN=5) :: c_glob_ib
  INTEGER :: iunout, ierr
  CHARACTER(LEN=256) :: fname
  !
  ! Generate the filename
  ! 
  WRITE(c_glob_iks,'(i5.5)') glob_iks
  WRITE(c_glob_ib, '(i5.5)') glob_ib
  !
  fname = TRIM( wfreq_dirname )//"/over_"//labellina//"_K"//c_glob_iks//"B"//c_glob_ib//".dat"
  !
  IF ( mpime == root ) THEN
     !
     ! ... open XML descriptor
     !
     CALL iotk_free_unit( iunout, ierr )
     CALL iotk_open_read( iunout, FILE = fname , BINARY = .TRUE., IERR = ierr )
     !
  END IF
  !
  CALL mp_bcast( ierr, root, world_comm )
  CALL errore( 'freq', 'cannot open freq file for reading', ierr )
  !
  IF ( mpime == root ) THEN  
     !
     CALL iotk_scan_begin( iunout, "OVER" )
     CALL iotk_scan_dat( iunout, "overlap", overlap )
     CALL iotk_scan_end( iunout, "OVER" )
     !
     ! ... close XML descriptor
     !
     CALL iotk_close_read( iunout )
     !
  END IF
  !
  CALL mp_bcast( overlap, root, world_comm )
  !
END SUBROUTINE
!
!
!
!
SUBROUTINE writeout_solvehf( sigma_hf, nb, nk )
  !
  ! WHO WRITES? ... root
  ! WHAT?       ... sigma_hf(nb, nk) 
  !
  USE kinds,          ONLY : DP
  USE westcom,        ONLY : wfreq_dirname
  USE mp_world,       ONLY : mpime,root,world_comm 
  USE mp,             ONLY : mp_bcast
  USE iotk_module
  !
  IMPLICIT NONE
  !
  ! I/O 
  !
  REAL(DP),INTENT(IN) :: sigma_hf(nb,nk)
  INTEGER,INTENT(IN) :: nb, nk
  !
  ! Workspace
  !
  INTEGER :: iunout, ierr
  CHARACTER(LEN=256) :: fname
  !
  ! Generate the filename
  !
  fname = TRIM( wfreq_dirname )//"/hf.dat"
  !
  IF ( mpime == root ) THEN
     !
     ! ... open XML descriptor
     !
     CALL iotk_free_unit( iunout, ierr )
     CALL iotk_open_write( iunout, FILE = fname , BINARY = .TRUE., IERR = ierr )
     !
  END IF
  !
  CALL mp_bcast( ierr, root, world_comm )
  CALL errore( 'freq', 'cannot open freq file for writing', ierr )
  !
  IF ( mpime == root ) THEN  
     !
     CALL iotk_write_begin( iunout, "HF_DB" )
     CALL iotk_write_dat( iunout, "sigma_hf", sigma_hf )
     CALL iotk_write_end( iunout, "HF_DB" )
     !
     ! ... close XML descriptor
     !
     CALL iotk_close_write( iunout )
     !
  END IF
  !
END SUBROUTINE
!
!
!
!
SUBROUTINE readin_solvehf( sigma_hf, nb, nk )
  !
  ! WHO WRITES? ... root
  ! WHAT?       ... sigma_hf(nb, nk) 
  !
  USE kinds,          ONLY : DP
  USE westcom,        ONLY : iks_l2g,wfreq_dirname
  USE mp_world,       ONLY : mpime,root,world_comm 
  USE mp,             ONLY : mp_bcast
  USE iotk_module
  !
  IMPLICIT NONE
  !
  ! I/O 
  !
  REAL(DP),INTENT(OUT) :: sigma_hf(nb,nk)
  INTEGER,INTENT(IN) :: nb, nk
  !
  ! Workspace
  !
  INTEGER :: iunout, ierr
  CHARACTER(LEN=256) :: fname
  !
  ! Generate the filename
  ! 
  fname = TRIM( wfreq_dirname )//"/hf.dat"
  !
  IF ( mpime == root ) THEN
     !
     ! ... open XML descriptor
     !
     CALL iotk_free_unit( iunout, ierr )
     CALL iotk_open_read( iunout, FILE = fname , BINARY = .TRUE., IERR = ierr )
     !
  END IF
  !
  CALL mp_bcast( ierr, root, world_comm )
  CALL errore( 'freq', 'cannot open freq file for reading', ierr )
  !
  IF ( mpime == root ) THEN  
     !
     CALL iotk_scan_begin( iunout, "HF_DB" )
     CALL iotk_scan_dat( iunout, "sigma_hf", sigma_hf )
     CALL iotk_scan_end( iunout, "HF_DB" )
     !
     ! ... close XML descriptor
     !
     CALL iotk_close_read( iunout )
     !
  END IF
  !
  CALL mp_bcast( sigma_hf, root, world_comm )
  !
END SUBROUTINE
