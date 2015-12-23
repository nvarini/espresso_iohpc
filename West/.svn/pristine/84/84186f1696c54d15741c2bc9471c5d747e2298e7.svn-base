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
SUBROUTINE solve_wfreq(l_read_restart)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP 
  USE westcom,              ONLY : sqvc,west_prefix,n_pdep_eigen_to_use,n_lanczos,npwq0,l_macropol,iks_l2g,dmat_inv,zmat_inv,&
                                 & l_enable_lanczos,nbnd_occ,iuwfc,lrwfc,io_comm,wfreq_eta,imfreq_list,refreq_list,&
                                 & rhead,ihead,o_restart_time,l_skip_nl_part_of_hcomr
  USE mp_global,            ONLY : my_image_id,nimage,inter_image_comm,intra_bgrp_comm
  USE mp_world,             ONLY : mpime
  USE mp,                   ONLY : mp_bcast,mp_barrier,mp_sum
  USE io_global,            ONLY : stdout,ionode
  USE gvect,                ONLY : g,ngm,gstart,ig_l2g
  USE cell_base,            ONLY : tpiba2,bg,omega
  USE fft_base,             ONLY : dffts
  USE constants,            ONLY : tpi,fpi,e2
  USE pwcom,                ONLY : npw,npwx,et,nks,current_spin,isk,xk,nbnd,lsda,igk,g2kin,ecutwfc,current_k,wk
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
  USE distribution_center,  ONLY : pert,macropert,index_distribution,copy_container,ifr,rfr
  USE wfreq_restart,        ONLY : solvewfreq_restart_write,solvewfreq_restart_read,bks_type
  USE linear_algebra_kernel, ONLY : matdiago_dsy
  USE io_eigenfreq,         ONLY : eigenfreq_write
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  LOGICAL,INTENT(IN) :: l_read_restart
  !
  ! Workspace
  !
  INTEGER :: i1,i2,i3,im,ip,ig,glob_ip,ir,iv,iks,ipol,m
  CHARACTER(LEN=256)    :: wstat_dirname, fname
  CHARACTER(LEN=6)      :: my_label_b
  COMPLEX(DP),ALLOCATABLE :: auxr(:)
  INTEGER :: nbndval
  REAL(DP),ALLOCATABLE :: diago( :, : ), subdiago( :, :), bnorm(:), braket(:, :, :)
  COMPLEX(DP),ALLOCATABLE :: q_s( :, :, : )
  INTEGER :: info
  COMPLEX(DP),ALLOCATABLE :: dvpsi(:,:)
  COMPLEX(DP),ALLOCATABLE :: phi(:,:)
  COMPLEX(DP),ALLOCATABLE :: phi_tmp(:,:)
  REAL(DP) :: anorm
  LOGICAL :: conv_root
  COMPLEX(DP),ALLOCATABLE :: pertg(:),pertr(:)
  COMPLEX(DP) :: zkonstant
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  COMPLEX(DP) :: ps_c
  REAL(DP) :: ps_r
  REAL(DP),EXTERNAL :: DDOT 
  COMPLEX(DP),EXTERNAL :: ZDOTC 
  TYPE(index_distribution) :: mypara
  REAL(DP),ALLOCATABLE :: overlap(:,:)
  REAL(DP) :: mwo, ecv, dfactor, frequency, dhead
  COMPLEX(DP) :: zmwo, zfactor,zm,zp, zhead
  INTEGER :: glob_jp,ic,ifreq,il
  REAL(DP),ALLOCATABLE :: dmatilda(:,:), dlambda(:,:)
  COMPLEX(DP),ALLOCATABLE :: zmatilda(:,:), zlambda(:,:)
  REAL(DP),ALLOCATABLE :: dmat(:,:,:)
  COMPLEX(DP),ALLOCATABLE :: zmat(:,:,:)
  LOGICAL :: l_iks_skip, l_iv_skip
  REAL(DP),ALLOCATABLE :: eigenfreq(:)
  REAL(DP) :: time_spent(2)
  REAL(DP),EXTERNAL :: get_clock
  TYPE(bks_type) :: bks
  !
  CALL io_push_title("(W)-Lanczos")
  ! 
  ! DISTRIBUTING...
  !
  IF(l_macropol) THEN
     CALL copy_container( macropert, mypara )
  ELSE
     CALL copy_container( pert,      mypara )
  ENDIF
  !
  ! This is to reduce memory
  !
  CALL deallocate_bec_type( becp ) 
  CALL allocate_bec_type ( nkb, mypara%nloc, becp ) ! I just need 2 becp at a time
  !
  ! ALLOCATE dmat, zmat, where chi0 is stored
  !
  ALLOCATE( dmat( mypara%nglob, mypara%nloc, ifr%nloc) )
  ALLOCATE( zmat( mypara%nglob, mypara%nloc, rfr%nloc) )
  dmat = 0._DP
  zmat = 0._DP
  !
  ! Remember the directory name
  !
  wstat_dirname = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.wstat.save'
  !
  IF(l_read_restart) THEN
     CALL solvewfreq_restart_read( bks, dmat, zmat, mypara%nglob, mypara%nloc )
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
     DO iv = 1, nbnd_occ(iks)
        IF(iks==bks%lastdone_ks .AND. iv <= bks%lastdone_band ) CYCLE
        barra_load = barra_load + 1
     ENDDO
  ENDDO
  ! 
  IF( barra_load == 0 ) THEN 
     CALL start_bar_type ( barra, 'wlanczos', 1 )
     CALL update_bar_type( barra, 'wlanczos', 1 ) 
  ELSE
     CALL start_bar_type ( barra, 'wlanczos', barra_load )
  ENDIF
  !
  ! LOOP 
  !
  DO iks = 1, nks   ! KPOINT-SPIN
     IF(iks<bks%lastdone_ks) CYCLE
     !
     mwo = - wk(iks) / omega
     zmwo = CMPLX( - wk(iks) / omega, 0._DP, KIND=DP)
     !
     CALL gk_sort(xk(1,iks),ngm,g,ecutwfc/tpiba2,npw,igk,g2kin)
     g2kin=g2kin*tpiba2
     !
     ! reads unperturbed wavefuctions psi_k in G_space, for all bands
     !
     IF(nks.NE.1) THEN
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
     bks%max_band = nbndval
     bks%min_band = 1
     !
     ALLOCATE(dvpsi(npwx*npol,mypara%nlocx)) 
     !
     time_spent(1) = get_clock( 'wlanczos' ) 
     !
     ! LOOP over band states 
     !
     DO iv = 1, nbndval
        IF(iks==bks%lastdone_ks .AND. iv <= bks%lastdone_band ) CYCLE
        !
        ! PSIC
        !
        IF(gamma_only) THEN
           CALL SINGLEBAND_invfft(npw,evc(1,iv),npwx,psic,dffts%nnr) 
        ELSEIF(noncolin) THEN
           CALL SINGLEBAND_invfft_k(npw,evc(1     ,iv),npwx,psic_nc(1,1),dffts%nnr,.TRUE.)
           CALL SINGLEBAND_invfft_k(npw,evc(1+npwx,iv),npwx,psic_nc(1,2),dffts%nnr,.TRUE.)
        ELSE
           CALL SINGLEBAND_invfft_k(npw,evc(1,iv),npwx,psic,dffts%nnr,.TRUE.)
        ENDIF
        !
        IF(l_macropol) THEN
           !
           ! PHI 
           !
           ALLOCATE(phi(npwx*npol,3))
           ALLOCATE(phi_tmp(npwx*npol,3))
           CALL commutator_Hx_psi (iks, 1, 1, evc(1,iv), phi_tmp(1,1), l_skip_nl_part_of_hcomr)
           CALL commutator_Hx_psi (iks, 1, 2, evc(1,iv), phi_tmp(1,2), l_skip_nl_part_of_hcomr)
           CALL commutator_Hx_psi (iks, 1, 3, evc(1,iv), phi_tmp(1,3), l_skip_nl_part_of_hcomr)
           phi = 0._DP
           DO i1 = 1, 3
              DO i2 = 1, 3 
                 zkonstant = CMPLX( bg(i1,i2) , KIND=DP )  
                 CALL ZAXPY(npwx*npol,zkonstant,phi_tmp(1,i2),1,phi(1,i1),1) 
              ENDDO
           ENDDO
           DEALLOCATE(phi_tmp)
           !
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
        DO ip=1,mypara%nloc
           !
           glob_ip = mypara%l2g(ip)
           !
           ! Decide whether read dbs E or dhpi 
           !
           IF(glob_ip<=n_pdep_eigen_to_use) THEN
              !
              ! Exhume dbs eigenvalue
              !
              WRITE(my_label_b,'(i6.6)') glob_ip
              fname = TRIM( wstat_dirname ) // "/E"//TRIM(ADJUSTL(my_label_b))//".dat"
              CALL pdep_read_G_and_distribute(fname,pertg)
              !
              ! Multiply by sqvc
              pertg(:) = sqvc(:) * pertg(:)
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
                    pertr(ir)=psic_nc(ir,1)*pertr(ir)
                 ENDDO
                 CALL SINGLEBAND_fwfft_k(npw,pertr,dffts%nnr,dvpsi(1,ip),npwx,.TRUE.)
                 CALL SINGLEBAND_invfft_k(npwq0,pertg(1),npwx,pertr,dffts%nnr,.FALSE.)
                 DO ir=1,dffts%nnr 
                    pertr(ir)=psic_nc(ir,2)*pertr(ir)
                 ENDDO
                 CALL SINGLEBAND_fwfft_k(npw,pertr,dffts%nnr,dvpsi(1+npwx,ip),npwx,.TRUE.)
              ELSE
                 CALL SINGLEBAND_invfft_k(npwq0,pertg(1),npwx,pertr,dffts%nnr,.FALSE.)
                 DO ir=1,dffts%nnr 
                    pertr(ir)=psic(ir)*pertr(ir)
                 ENDDO
                 CALL SINGLEBAND_fwfft_k(npw,pertr,dffts%nnr,dvpsi(1,ip),npwx,.TRUE.)
              ENDIF 
              !
           ELSE
              !
              ipol = glob_ip-n_pdep_eigen_to_use
              !
              dvpsi(:,ip) = phi(:,ipol) * DSQRT(fpi * e2)
              !
           ENDIF
           !
        ENDDO ! pert
        ! 
        DEALLOCATE(pertr)
        DEALLOCATE(pertg)
        IF(l_macropol) THEN
           DEALLOCATE(phi)
        ENDIF
        !
        CALL apply_alpha_pc_to_m_wfcs(nbndval,mypara%nloc,dvpsi,(1._DP,0._DP))
        !
        ! OVERLAP( glob_ip, im=1:nbnd ) = < psi_im iks | dvpsi_glob_ip >
        !
        IF(ALLOCATED(overlap)) DEALLOCATE(overlap) 
        ALLOCATE(overlap(mypara%nglob, nbndval+1:nbnd ) )
        !
        overlap = 0._DP
        DO ic = nbndval+1, nbnd
           DO ip = 1, mypara%nloc
              !
              glob_ip = mypara%l2g(ip)
              !
              IF( gamma_only ) THEN
                 !
                 ps_r = 2._DP * DDOT( 2*npw, evc(1,ic), 1, dvpsi(1,ip), 1)
                 IF (gstart==2) ps_r = ps_r - REAL(evc(1,ic),KIND=DP) * REAL(dvpsi(1,ip),KIND=DP)
                 overlap(glob_ip,ic) = ps_r
                 !
              ELSE
                 !
                 ps_c = ZDOTC( npw, evc(1,ic), 1, dvpsi(1,ip), 1)
                 IF(noncolin) THEN
                    ps_c = ps_c + ZDOTC( npw, evc(1+npwx,ic), 1, dvpsi(1+npwx,ip), 1)
                 ENDIF
                 overlap(glob_ip,ic) = REAL(ps_c,KIND=DP)
                 !
              ENDIF
              IF( glob_ip > n_pdep_eigen_to_use ) overlap(glob_ip,ic) = overlap(glob_ip,ic) / ( et(ic,iks)-et(iv,iks) )
           ENDDO
        ENDDO
        !
        CALL mp_sum(overlap,intra_bgrp_comm)
        CALL mp_sum(overlap,inter_image_comm)
        !
        ! Update dmat with cond
        !
        DO ifreq = 1, ifr%nloc
           !
           frequency = imfreq_list( ifreq ) 
           !
           DO ic = nbndval+1, nbnd 
              !
              ecv = et(ic,iks)-et(iv,iks)
              dfactor = mwo * 2._DP * ecv / ( ecv**2 + frequency**2 )
              !
              DO ip = 1, mypara%nloc
                 glob_ip = mypara%l2g(ip)
                 DO glob_jp = 1, mypara%nglob
                    !
                    dmat(glob_jp,ip,ifreq) = dmat(glob_jp,ip,ifreq) + overlap( glob_jp, ic) * overlap( glob_ip, ic) * dfactor
                    !
                 ENDDO
              ENDDO
              !
           ENDDO ! ic
        ENDDO ! ifreq
        !
        ! Update zmat with cond
        !
        DO ifreq = 1, rfr%nloc
           !
           frequency = refreq_list( ifreq ) 
           !
           DO ic = nbndval+1, nbnd 
              !
              ecv = et(ic,iks)-et(iv,iks)
              zp = CMPLX( ecv + frequency, - wfreq_eta, KIND=DP )
              zm = CMPLX( ecv - frequency, - wfreq_eta, KIND=DP )
              zfactor = zmwo / zp + zmwo / zm 
              !
              DO ip = 1, mypara%nloc
                 glob_ip = mypara%l2g(ip)
                 DO glob_jp = 1, mypara%nglob
                    !
                    zmat(glob_jp,ip,ifreq) = zmat(glob_jp,ip,ifreq) + CMPLX( overlap( glob_jp, ic) * &
                    & overlap( glob_ip, ic), 0._DP, KIND=DP ) * zfactor
                    !
                 ENDDO
              ENDDO
              !
           ENDDO ! ic
        ENDDO ! ifreq 
        !
        DEALLOCATE(overlap)
        !
        ! Apply Pc, to be sure
        !
        CALL apply_alpha_pc_to_m_wfcs(nbnd,mypara%nloc,dvpsi,(1._DP,0._DP))
        !
        ! Now dvpsi is distributed according to eigen_distr (image), I need to use it for lanczos
        ! In the gamma_only case I need to process 2 dvpsi at a time (+ the odd last one, eventually), otherwise 1 at a time.
        !
        IF( l_enable_lanczos ) THEN  
           !
           ALLOCATE( bnorm    (                             mypara%nloc ) )
           ALLOCATE( diago    (               n_lanczos   , mypara%nloc ) )
           ALLOCATE( subdiago (               n_lanczos-1 , mypara%nloc ) )
           ALLOCATE( q_s      ( npwx*npol   , mypara%nloc , n_lanczos   ) )  ! WARNING ORDER INVERTED TO SMOOTHEN LANCZOS ALGORITHM 
           !
           CALL solve_deflated_lanczos_w_full_ortho ( nbnd, mypara%nloc, n_lanczos, dvpsi, diago, subdiago, q_s, bnorm)
           !
           ALLOCATE( braket   ( mypara%nglob, n_lanczos   , mypara%nloc ) )
           CALL get_brak_hyper_parallel(dvpsi,mypara%nloc,n_lanczos,q_s,braket,mypara%nloc,mypara%nlocx,mypara%nglob)
           DEALLOCATE( q_s )
           !
           DO ip = 1, mypara%nloc
              CALL diago_lanczos( bnorm(ip), diago( :, ip), subdiago( :, ip), braket(:,:,ip), mypara%nglob )
           ENDDO
           !
           DEALLOCATE( bnorm )
           DEALLOCATE( subdiago )
           !
           ! Update dmat with lanczos
           !
           DO ifreq = 1, ifr%nloc
              !
              frequency = imfreq_list( ifreq ) 
              !
              DO il = 1, n_lanczos 
                 !
                 DO ip = 1, mypara%nloc
                    glob_ip = mypara%l2g(ip)
                    ecv = diago( il, ip ) - et(iv,iks) 
                    dfactor = mwo * 2._DP * ecv / ( ecv**2 + frequency**2 )
                    IF( glob_ip > n_pdep_eigen_to_use ) dfactor = dfactor / ecv
                    DO glob_jp = 1, mypara%nglob
                       !
                       IF( glob_jp > n_pdep_eigen_to_use ) dfactor = dfactor / ecv
                       dmat(glob_jp,ip,ifreq) = dmat(glob_jp,ip,ifreq) + braket( glob_jp, il, ip) * dfactor
                       !
                    ENDDO
                 ENDDO
                 !
              ENDDO ! il
           ENDDO ! ifreq
           !
           ! Update zmat with lanczos
           !
           DO ifreq = 1, rfr%nloc
              !
              frequency = refreq_list( ifreq ) 
              !
              DO il = 1, n_lanczos 
                 !
                 DO ip = 1, mypara%nloc
                    glob_ip = mypara%l2g(ip)
                    ecv = diago( il, ip ) - et(iv,iks) 
                    zp = CMPLX( ecv + frequency, - wfreq_eta, KIND=DP )
                    zm = CMPLX( ecv - frequency, - wfreq_eta, KIND=DP )
                    zfactor = zmwo / zp + zmwo / zm
                    IF( glob_ip > n_pdep_eigen_to_use ) zfactor = zfactor / CMPLX( ecv, 0._DP, KIND=DP)
                    DO glob_jp = 1, mypara%nglob
                       !
                       IF( glob_jp > n_pdep_eigen_to_use ) zfactor = zfactor / CMPLX( ecv, 0._DP, KIND=DP)
                       zmat(glob_jp,ip,ifreq) = zmat(glob_jp,ip,ifreq) + CMPLX( braket( glob_jp, il, ip), 0._DP, KIND=DP ) * zfactor
                       !
                    ENDDO
                 ENDDO
                 !
              ENDDO ! il
           ENDDO ! ifreq
           !
           ! MPI-IO 
           !
           !CALL writeout_solvewfreq( iks_l2g(iks), iv, diago, braket, io_comm, mypara%nloc, mypara%nglob, mypara%myoffset )
           !
           DEALLOCATE( diago ) 
           DEALLOCATE( braket )
           !
        ENDIF ! l_enable_lanczos
        !
        time_spent(2) = get_clock( 'wlanczos' ) 
        !
        IF( o_restart_time >= 0._DP ) THEN 
           IF( (time_spent(2)-time_spent(1)) > o_restart_time*60._DP .OR. iv == nbndval ) THEN 
              bks%lastdone_ks=iks
              bks%lastdone_band=iv
              CALL solvewfreq_restart_write(bks,dmat,zmat,mypara%nglob,mypara%nloc)
              bks%old_ks = iks 
              bks%old_band = iv
              time_spent(1) = get_clock( 'wlanczos' )
           ENDIF
        ENDIF
        !
        CALL update_bar_type( barra, 'wlanczos', 1 )
        !
     ENDDO ! BANDS
     !
     DEALLOCATE(dvpsi)
     !
  ENDDO ! KPOINT-SPIN 
  !
  CALL stop_bar_type( barra, 'wlanczos' )
  !
  ! EPS-1 imfreq
  !
  ALLOCATE( dmatilda( mypara%nglob, mypara%nglob ) )
  ALLOCATE( dlambda( n_pdep_eigen_to_use, n_pdep_eigen_to_use) ) 
  ALLOCATE( dmat_inv( pert%nglob, pert%nloc, ifr%nloc) )
  dmat_inv = 0._DP
  IF(l_macropol) ALLOCATE( ihead( ifr%nloc) )
  !
  CALL start_clock('chi0_diag')
  !
  DO ifreq = 1, ifr%nloc
     !
     dmatilda = 0._DP
     DO ip = 1, mypara%nloc
        glob_ip = mypara%l2g(ip) 
        dmatilda( :, glob_ip) = dmat( :, ip, ifreq )
     ENDDO
     !
     CALL mp_sum( dmatilda, inter_image_comm )
     ! 
     CALL chi_invert_real( dmatilda, dhead, dlambda, mypara%nglob)
     !
     DO ip = 1, pert%nloc
        glob_ip = pert%l2g(ip)
        dmat_inv(1:n_pdep_eigen_to_use,ip,ifreq) = dlambda( 1:n_pdep_eigen_to_use, glob_ip)
     ENDDO 
     IF( l_macropol ) ihead( ifreq) = dhead
     !
  ENDDO
  
  CALL stop_clock('chi0_diag')
  
  !
  DEALLOCATE( dlambda )
  DEALLOCATE( dmatilda )
  DEALLOCATE(dmat)
  !
  ! EPS-1 refreq
  !
  ALLOCATE( zmatilda( mypara%nglob, mypara%nglob ) )
  ALLOCATE( zlambda( n_pdep_eigen_to_use, n_pdep_eigen_to_use ) ) 
  ALLOCATE( zmat_inv( pert%nglob, pert%nloc, rfr%nloc) )
  zmat_inv = 0._DP
  IF(l_macropol) ALLOCATE( rhead( rfr%nloc) )
  !
  DO ifreq = 1, rfr%nloc
     !
     zmatilda = 0._DP
     DO ip = 1, mypara%nloc
        glob_ip = mypara%l2g(ip) 
        zmatilda( :, glob_ip) = zmat( :, ip, ifreq )
     ENDDO
     !
     CALL mp_sum( zmatilda, inter_image_comm ) 
     CALL chi_invert_complex( zmatilda, zhead, zlambda, mypara%nglob)
     !
     DO ip = 1, pert%nloc
        glob_ip = pert%l2g(ip)
        zmat_inv(1:n_pdep_eigen_to_use,ip,ifreq) = zlambda( 1:n_pdep_eigen_to_use, glob_ip)
     ENDDO 
     IF( l_macropol ) rhead( ifreq) = zhead
     !
  ENDDO
  !
  DEALLOCATE( zlambda )
  DEALLOCATE( zmatilda )
  DEALLOCATE( zmat ) 
  !
END SUBROUTINE 
