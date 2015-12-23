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
SUBROUTINE solve_gfreq(l_read_restart)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP 
  USE westcom,              ONLY : sqvc,west_prefix,n_pdep_eigen_to_use,n_lanczos,npwq0,qp_bandrange,iks_l2g,&
                                 & wfreq_dirname,l_enable_lanczos,nbnd_occ,iuwfc,lrwfc,io_comm,o_restart_time
  USE mp_global,            ONLY : my_image_id,nimage,inter_image_comm,intra_bgrp_comm,inter_pool_comm
  USE mp,                   ONLY : mp_bcast,mp_barrier,mp_sum
  USE io_global,            ONLY : stdout, ionode
  USE gvect,                ONLY : g,ngm,gstart
  USE cell_base,            ONLY : tpiba2,bg
  USE fft_base,             ONLY : dffts
  USE constants,            ONLY : tpi,fpi,e2
  USE pwcom,                ONLY : npw,npwx,et,nks,current_spin,isk,xk,nbnd,lsda,igk,g2kin,ecutwfc,nkstot,current_k
  USE wavefunctions_module, ONLY : evc,psic,psic_nc
  USE io_files,             ONLY : tmp_dir,nwordwfc,iunwfc
  USE fft_at_gamma,         ONLY : DOUBLEBAND_INVFFT,SINGLEBAND_INVFFT,DOUBLEBAND_FWFFT,SINGLEBAND_FWFFT
  USE fft_at_k,             ONLY : SINGLEBAND_INVFFT_k,SINGLEBAND_FWFFT_k
  USE becmod,               ONLY : becp,allocate_bec_type,deallocate_bec_type
  USE uspp,                 ONLY : vkb,nkb
  USE pdep_io,              ONLY : pdep_read_G_and_distribute 
  USE io_push,              ONLY : io_push_title
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin,npol 
  USE buffers,              ONLY : get_buffer
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE distribution_center,  ONLY : pert
  USE wfreq_restart,        ONLY : solvegfreq_restart_write,solvegfreq_restart_read,bks_type
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  LOGICAL,INTENT(IN) :: l_read_restart
  !
  ! Workspace
  !
  INTEGER :: i1,i2,i3,ip,ig,glob_ip,ir,ib,iks,m,im
  CHARACTER(LEN=256)    :: wstat_dirname, fname
  CHARACTER(LEN=6)      :: my_label_b
  COMPLEX(DP),ALLOCATABLE :: auxr(:)
  !INTEGER :: iuwfc, lrwfc
  INTEGER :: nbndval
  REAL(DP),ALLOCATABLE :: diago( :, : ), subdiago( :, :), bnorm(:), braket(:, :, :)
  COMPLEX(DP),ALLOCATABLE :: q_s( :, :, : )
  INTEGER :: info
  COMPLEX(DP),ALLOCATABLE :: dvpsi(:,:)
  REAL(DP) :: anorm
  LOGICAL :: conv_root
  COMPLEX(DP),ALLOCATABLE :: pertg(:),pertr(:)
  COMPLEX(DP) :: ps_c
  REAL(DP) :: ps_r
  REAL(DP),EXTERNAL :: DDOT 
  COMPLEX(DP),EXTERNAL :: ZDOTC 
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  REAL(DP),ALLOCATABLE :: overlap(:,:)
  LOGICAL :: l_iks_skip, l_ib_skip
  REAL(DP) :: time_spent(2)
  REAL(DP),EXTERNAL :: get_clock
  TYPE(bks_type) :: bks
  !
  CALL io_push_title("(G)-Lanczos")
  !
  ! This is to reduce memory
  !
  CALL deallocate_bec_type( becp ) 
  CALL allocate_bec_type ( nkb, pert%nloc, becp ) ! I just need 2 becp at a time
  !
  ! Remember the directory names
  !
  wstat_dirname = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.wstat.save'
  !
  IF(l_read_restart) THEN
     CALL solvegfreq_restart_read( bks )
  ELSE
     bks%lastdone_ks   = 0 
     bks%lastdone_band = 0 
     bks%old_ks        = 0 
     bks%old_band      = 0 
     bks%max_ks        = nks 
     bks%min_ks        = 1 
  ENDIF
  !
  barra_load = 0
  DO iks = 1, nks
     IF(iks<bks%lastdone_ks) CYCLE
     DO ib = qp_bandrange(1), qp_bandrange(2)
        IF(iks==bks%lastdone_ks .AND. ib <= bks%lastdone_band ) CYCLE
        barra_load = barra_load + 1
     ENDDO
  ENDDO
  ! 
  IF( barra_load == 0 ) THEN 
     CALL start_bar_type( barra, 'glanczos', 1 )
     CALL update_bar_type( barra, 'glanczos', 1 ) 
  ELSE
     CALL start_bar_type( barra, 'glanczos', barra_load )
  ENDIF
  !
  ! LOOP 
  !
  DO iks = 1, nks   ! KPOINT-SPIN
     IF(iks<bks%lastdone_ks) CYCLE
     !
     CALL gk_sort(xk(1,iks),ngm,g,ecutwfc/tpiba2,npw,igk,g2kin)
     g2kin=g2kin*tpiba2
     !
     ! reads unperturbed wavefuctions psi_k in G_space, for all bands
     !
     IF(nks.NE.1) THEN
        !CALL get_buffer( evc, nwordwfc, iunwfc, iks )
        !iuwfc = 20
        !lrwfc = nbnd * npwx * npol
        IF(my_image_id==0) CALL get_buffer(evc, lrwfc, iuwfc, iks)
        CALL mp_bcast(evc,0,inter_image_comm)
     ENDIF
     !
     CALL init_us_2 (npw, igk, xk (1, iks), vkb)
     !
     current_k = iks
     current_spin = isk(iks)
     !
     nbndval = nbnd_occ(iks)
     !
     bks%max_band=nbndval
     bks%min_band=1
     !
     ALLOCATE(dvpsi(npwx*npol,pert%nlocx))   
     CALL preallocate_solvegfreq( iks_l2g(iks), qp_bandrange(1), qp_bandrange(2), pert )
     !
     time_spent(1) = get_clock( 'glanczos' ) 
     !
     ! LOOP over band states 
     !
     DO ib = qp_bandrange(1), qp_bandrange(2)
        IF(iks==bks%lastdone_ks .AND. ib <= bks%lastdone_band ) CYCLE
        !
        ! PSIC
        !
        IF(gamma_only) THEN
           CALL SINGLEBAND_invfft(npw,evc(1,ib),npwx,psic,dffts%nnr) 
        ELSEIF(noncolin) THEN
           CALL SINGLEBAND_invfft_k(npw,evc(1     ,ib),npwx,psic_nc(1,1),dffts%nnr,.TRUE.)
           CALL SINGLEBAND_invfft_k(npw,evc(1+npwx,ib),npwx,psic_nc(1,2),dffts%nnr,.TRUE.)
        ELSE
           CALL SINGLEBAND_invfft_k(npw,evc(1,ib),npwx,psic,dffts%nnr,.TRUE.)
        ENDIF
        !
        ! ZEROS
        !
        dvpsi = 0._DP   
        !
        ! Read PDEP
        !
        ALLOCATE( pertg(npwq0) )
        ALLOCATE( pertr( dffts%nnr ) )
        !
        DO ip=1,pert%nloc
           glob_ip = pert%l2g(ip)
           !
           ! Exhume dbs eigenvalue
           !
           WRITE(my_label_b,'(i6.6)') glob_ip
           fname = TRIM( wstat_dirname ) // "/E"//TRIM(ADJUSTL(my_label_b))//".dat"
           CALL pdep_read_G_and_distribute(fname,pertg)
           !
           ! Multiply by sqvc
           pertg(:) = sqvc(:) * pertg(:) ! / SQRT(fpi*e2)     ! CONTROLLARE QUESTO
           !
           ! Bring it to R-space
           IF(gamma_only) THEN
              CALL SINGLEBAND_invfft(npwq0,pertg(1),npwx,pertr,dffts%nnr)
              DO ir=1,dffts%nnr 
                 pertr(ir)=psic(ir)*pertr(ir)
              ENDDO
              CALL SINGLEBAND_fwfft(npw,pertr,dffts%nnr,dvpsi(1,ip),npwx)
           ELSEIF(noncolin) THEN
              CALL SINGLEBAND_invfft_k(npwq0,pertg(1),npwx,pertr,dffts%nnr,.FALSE.)
              DO ir=1,dffts%nnr 
                 pertr(ir)=psic_nc(ir,1)*DCONJG(pertr(ir))
              ENDDO
              CALL SINGLEBAND_fwfft_k(npw,pertr,dffts%nnr,dvpsi(1,ip),npwx,.TRUE.)
              CALL SINGLEBAND_invfft_k(npwq0,pertg(1),npwx,pertr,dffts%nnr,.FALSE.)
              DO ir=1,dffts%nnr 
                 pertr(ir)=psic_nc(ir,2)*DCONJG(pertr(ir))
              ENDDO
              CALL SINGLEBAND_fwfft_k(npw,pertr,dffts%nnr,dvpsi(1+npwx,ip),npwx,.TRUE.)
           ELSE
              CALL SINGLEBAND_invfft_k(npwq0,pertg(1),npwx,pertr,dffts%nnr,.FALSE.)
              DO ir=1,dffts%nnr 
                 pertr(ir)=psic(ir)*DCONJG(pertr(ir))
              ENDDO
              CALL SINGLEBAND_fwfft_k(npw,pertr,dffts%nnr,dvpsi(1,ip),npwx,.TRUE.)
           ENDIF 
           !
           !
        ENDDO ! pert
        ! 
        DEALLOCATE(pertr)
        DEALLOCATE(pertg)
        !
        ! OVERLAP( glob_ip, im=1:n_hstates ) = < psi_im iks | dvpsi_glob_ip >
        !
        IF(ALLOCATED(overlap)) DEALLOCATE(overlap) 
        ALLOCATE(overlap(pert%nglob, nbnd ) )
        !
        overlap = 0._DP
        DO im = 1, nbnd
           DO ip = 1, pert%nloc
              glob_ip = pert%l2g(ip)
              IF( gamma_only ) THEN
                 !
                 ps_r = 2._DP * DDOT( 2*npw, evc(1,im), 1, dvpsi(1,ip), 1)
                 IF (gstart==2) ps_r = ps_r - REAL(evc(1,im),KIND=DP) * REAL(dvpsi(1,ip),KIND=DP)
                 overlap(glob_ip,im) = ps_r
                 !
              ELSE
                 !
                 ps_c = ZDOTC( npw, evc(1,im), 1, dvpsi(1,ip), 1)
                 IF(noncolin) THEN
                    ps_c = ps_c + ZDOTC( npw, evc(1+npwx,im), 1, dvpsi(1+npwx,ip), 1)
                 ENDIF
                 overlap(glob_ip,im) = REAL(ps_c,KIND=DP)
                 !
              ENDIF
           ENDDO
        ENDDO
        !
        CALL mp_sum(overlap,intra_bgrp_comm)
        CALL mp_sum(overlap,inter_image_comm)
        CALL writeout_overlap( 'g', iks_l2g(iks), ib, overlap, pert%nglob, nbnd )
        DEALLOCATE(overlap)
        !
        CALL apply_alpha_pc_to_m_wfcs(nbnd,pert%nloc,dvpsi,(1._DP,0._DP))
        !
        ! Now dvpsi is distributed according to eigen_distr (image), I need to use it for lanczos
        ! In the gamma_only case I need to process 2 dvpsi at a time (+ the odd last one, eventually), otherwise 1 at a time.
        !
        IF( l_enable_lanczos ) THEN  
           !
           ALLOCATE( bnorm    (                             pert%nloc ) )
           ALLOCATE( diago    (               n_lanczos   , pert%nloc ) )
           ALLOCATE( subdiago (               n_lanczos-1 , pert%nloc ) )
           ALLOCATE( q_s      ( npwx*npol   , pert%nloc   , n_lanczos   ) )  ! WARNING ORDER INVERTED TO SMOOTHEN LANCZOS ALGORITHM 
           !
           CALL solve_deflated_lanczos_w_full_ortho ( nbnd, pert%nloc, n_lanczos, dvpsi, diago, subdiago, q_s, bnorm)
           !
           ALLOCATE( braket   ( pert%nglob, n_lanczos   , pert%nloc ) )
           CALL get_brak_hyper_parallel(dvpsi,pert%nloc,n_lanczos,q_s,braket,pert%nloc,pert%nlocx,pert%nglob)
           DEALLOCATE( q_s )
           !
           DO ip = 1, pert%nloc
              CALL diago_lanczos( bnorm(ip), diago( :, ip), subdiago( :, ip), braket(:,:,ip), pert%nglob )
           ENDDO
           !
           DEALLOCATE( bnorm )
           DEALLOCATE( subdiago )
           !
           ! MPI-IO 
           !
           CALL writeout_solvegfreq( iks_l2g(iks), ib, diago, braket, io_comm, pert%nloc, pert%nglob, pert%myoffset )
           !
           DEALLOCATE( diago ) 
           DEALLOCATE( braket )
           !
        ENDIF ! l_enable_lanczos
        !
        time_spent(2) = get_clock( 'glanczos' ) 
        !
        IF( o_restart_time >= 0._DP ) THEN 
           IF( (time_spent(2)-time_spent(1)) > o_restart_time*60._DP .OR. ib == qp_bandrange(2) ) THEN 
              bks%lastdone_ks=iks
              bks%lastdone_band=ib 
              CALL solvegfreq_restart_write( bks )
              bks%old_ks=iks
              bks%old_band=ib
              time_spent(1) = get_clock( 'glanczos' )
           ENDIF
        ENDIF
        !
        CALL update_bar_type( barra, 'glanczos', 1 )
        !
     ENDDO ! BANDS
     !
     DEALLOCATE(dvpsi)
     !
  ENDDO ! KPOINT-SPIN 
  ! 
  CALL stop_bar_type( barra, 'glanczos' )
  !
END SUBROUTINE 
