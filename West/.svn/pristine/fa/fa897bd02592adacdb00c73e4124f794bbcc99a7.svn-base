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
SUBROUTINE dfpt (m,dvg,dng)
  !-----------------------------------------------------------------------
  !
  USE kinds,                 ONLY : DP
  USE io_global,             ONLY : stdout
  USE wvfct,                 ONLY : nbnd,g2kin,et
  USE fft_base,              ONLY : dfftp,dffts
  USE gvect,                 ONLY : nl,gstart,ig_l2g
  USE wavefunctions_module,  ONLY : evc,psic
  USE gvecs,                 ONLY : ngms
  USE mp,                    ONLY : mp_sum,mp_barrier,mp_bcast
  USE mp_global,             ONLY : inter_image_comm,inter_pool_comm,my_image_id
  USE fft_at_gamma,          ONLY : DOUBLEBAND_INVFFT,SINGLEBAND_INVFFT,DOUBLEBAND_FWFFT,SINGLEBAND_FWFFT
  USE fft_at_k,              ONLY : SINGLEBAND_INVFFT_k,SINGLEBAND_FWFFT_k
  USE buffers,               ONLY : get_buffer
  USE noncollin_module,      ONLY : noncolin,npol
  USE bar,                   ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE pwcom,                 ONLY : current_spin,wk,nks,nelup,neldw,isk,ecutwfc,g,igk,ngm,tpiba2,xk,omega,npw,npwx,lsda,nkstot,&
                                  & current_k
  USE control_flags,         ONLY : gamma_only, io_level
  USE io_files,              ONLY : tmp_dir, nwordwfc, iunwfc, diropn
  USE uspp,                  ONLY : nkb, vkb, okvan
  USE constants,             ONLY : e2,fpi
  USE westcom,               ONLY : npwq0,sqvc,nbnd_occ,iuwfc,lrwfc
  USE io_push,               ONLY : io_push_title
  USE mp_world,              ONLY : mpime,world_comm
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER,INTENT(IN) :: m
  COMPLEX(DP),INTENT(IN) :: dvg(npwx,m)
  COMPLEX(DP),INTENT(OUT) :: dng(npwx,m)
  !
  ! Workspace
  !
  INTEGER :: ipert, ig, ir, ibnd, iks
  INTEGER :: nbndval
  !
  REAL(DP) :: anorm, prod 
  REAL(DP),ALLOCATABLE :: h_diag(:,:)
  !
  COMPLEX(DP),ALLOCATABLE :: dvpsi(:,:),dpsi(:,:)
  COMPLEX(DP),ALLOCATABLE :: aux_r(:),aux_g(:)
  COMPLEX(DP),ALLOCATABLE :: dpsic(:)
  !
  TYPE(bar_type) :: barra
  !
  LOGICAL :: conv_dfpt
  LOGICAL :: exst,exst_mem
  !
  CALL mp_barrier( world_comm )
  !
  CALL report_dynamical_memory()
  !
  CALL io_push_title("Sternheimer eq. solver...")
  !
  dng=0.0_DP
  !
  CALL start_bar_type( barra, 'dfpt', MAX(m,1) * nks ) 
  !
  !IF(nks>1) CALL diropn(iuwfc,'wfc',lrwfc,exst) 
  !
  DO iks = 1, nks  ! KPOINT-SPIN LOOP
     !
     CALL gk_sort(xk(1,iks),ngm,g,ecutwfc/tpiba2,npw,igk,g2kin)
     g2kin=g2kin*tpiba2
     !
     ! reads unperturbed wavefuctions psi_k in G_space, for all bands
     !
     IF(nks>1) THEN
        !iuwfc = 20
        !lrwfc = nbnd * npwx * npol 
        !!CALL get_buffer( evc, nwordwfc, iunwfc, iks )
        IF(my_image_id==0) CALL get_buffer( evc, lrwfc, iuwfc, iks )
        !CALL mp_bcast(evc,0,inter_image_comm)
        !CALL davcio(evc,lrwfc,iuwfc,iks,-1)
        CALL mp_bcast(evc,0,inter_image_comm)
     ENDIF
     !
     CALL init_us_2 (npw, igk, xk (1, iks), vkb)
     !
     current_k = iks
     current_spin = isk(iks)
     !
     nbndval = nbnd_occ(iks)
     IF(nbndval==0) THEN
        CALL update_bar_type( barra,'dfpt', MAX(1,m) )
        CYCLE
     ENDIF
     !
     ALLOCATE(h_diag(npwx*npol,nbndval))
     CALL set_preconditioning(h_diag,nbndval)
     !
     ! set h_diag
     !
     ALLOCATE(dvpsi(npwx*npol,nbndval))   
     ALLOCATE(dpsi(npwx*npol,nbndval))
     !
     DO ipert = 1, m
        !
        ALLOCATE(aux_g(npwx))
        ALLOCATE(aux_r(dffts%nnr))
        !
        aux_g = 0._DP
        aux_r = 0._DP
        !
        DO ig = 1, npwq0  ! perturbation acts only on body
           aux_g(ig) = dvg(ig,ipert) * sqvc(ig) 
        ENDDO
        !
        IF(gamma_only) THEN
          CALL SINGLEBAND_invfft(npwq0,aux_g,npwx,aux_r,dffts%nnr)
        ELSE
          CALL SINGLEBAND_invfft_k(npwq0,aux_g,npwx,aux_r,dffts%nnr,.FALSE.)
        ENDIF
        !
        ! The perturbation is in aux_r
        !
        dvpsi=0._DP
        dpsi =0._DP
        !
        IF(gamma_only) THEN
           !
           ! double bands @ gamma
           DO ibnd=1,nbndval-MOD(nbndval,2),2
              !
              CALL DOUBLEBAND_invfft(npw,evc(1,ibnd),evc(1,ibnd+1),npwx,psic,dffts%nnr)
              DO ir=1,dffts%nnr
                 psic(ir) = psic(ir) * REAL(aux_r(ir),KIND=DP)
              ENDDO
              CALL DOUBLEBAND_fwfft(npw,psic,dffts%nnr,dvpsi(1,ibnd),dvpsi(1,ibnd+1),npwx)
              !
           ENDDO
           ! 
           ! single band @ gamma
           IF( MOD(nbndval,2) == 1 ) THEN
              ibnd=nbndval
              !
              CALL SINGLEBAND_invfft(npw,evc(1,ibnd),npwx,psic,dffts%nnr)
              DO ir=1,dffts%nnr
                 psic(ir) = CMPLX( REAL(psic(ir),KIND=DP) * REAL(aux_r(ir),KIND=DP), 0._DP, KIND=DP)
              ENDDO
              CALL SINGLEBAND_fwfft(npw,psic,dffts%nnr,dvpsi(1,ibnd),npwx) 
              !
           ENDIF
           !
        ELSE
           !
           ! only single bands
           DO ibnd=1,nbndval
              !
              CALL SINGLEBAND_invfft_k(npw,evc(1,ibnd),npwx,psic,dffts%nnr,.TRUE.)
              DO ir=1,dffts%nnr
                 psic(ir) = psic(ir) * aux_r(ir)
              ENDDO
              CALL SINGLEBAND_fwfft_k(npw,psic,dffts%nnr,dvpsi(1,ibnd),npwx,.TRUE.) 
              !
           ENDDO
           !
           IF(npol==2) THEN
              DO ibnd=1,nbndval
                 !
                 CALL SINGLEBAND_invfft_k(npw,evc(npwx+1,ibnd),npwx,psic,dffts%nnr,.TRUE.)
                 DO ir=1,dffts%nnr
                    psic(ir) = psic(ir) * aux_r(ir)
                 ENDDO
                 CALL SINGLEBAND_fwfft_k(npw,psic,dffts%nnr,dvpsi(npwx+1,ibnd),npwx,.TRUE.) 
                 !
              ENDDO
           ENDIF
           !
        ENDIF   
        !
        DEALLOCATE(aux_g)
        DEALLOCATE(aux_r)
        !
        CALL apply_alpha_pc_to_m_wfcs(nbndval,nbndval,dvpsi,(-1._DP,0._DP))
        !
        CALL linsolve_sternheimer_m_wfcts (nbndval, nbndval, dvpsi, dpsi, h_diag, et(1,iks), conv_dfpt, anorm) 
        !
        ! eliminare lter
        IF(.NOT.conv_dfpt) THEN 
           WRITE(stdout, '(7X,"** WARNING : PERT ",i8," not converged, res = ",E10.3)') ipert,anorm
        ENDIF
        !
        ALLOCATE(aux_r(dffts%nnr))
        !
        aux_r=0._DP
        ! 
        IF(gamma_only) THEN
           !
           ! double band @ gamma
           DO ibnd=1,nbndval
              !
              CALL DOUBLEBAND_invfft(npw, evc(1,ibnd), dpsi(1,ibnd),npwx, psic,dffts%nnr)
              DO ir=1,dffts%nnr
                 prod =  REAL( psic(ir),KIND=DP) * DIMAG( psic(ir)) 
                 aux_r(ir) = aux_r(ir) + CMPLX( prod, 0.0_DP, KIND=DP)
              ENDDO
              !
           ENDDO
           !
        ELSE
           !
           ALLOCATE(dpsic(dffts%nnr))
           !
           ! only single bands
           DO ibnd=1,nbndval 
              !
              CALL SINGLEBAND_invfft_k(npw, evc(1,ibnd),npwx, psic,dffts%nnr,.TRUE.)
              CALL SINGLEBAND_invfft_k(npw,dpsi(1,ibnd),npwx,dpsic,dffts%nnr,.TRUE.)
              DO ir=1,dffts%nnr 
                 aux_r(ir) = aux_r(ir) + DCONJG(psic(ir))*dpsic(ir)
              ENDDO
              !
           ENDDO
           !
           IF(npol==2) THEN
              DO ibnd=1,nbndval 
                 !
                 CALL SINGLEBAND_invfft_k(npw, evc(npwx+1,ibnd),npwx, psic,dffts%nnr,.TRUE.)
                 CALL SINGLEBAND_invfft_k(npw,dpsi(npwx+1,ibnd),npwx,dpsic,dffts%nnr,.TRUE.)
                 DO ir=1,dffts%nnr 
                    aux_r(ir) = aux_r(ir) + DCONJG(psic(ir))*dpsic(ir)
                 ENDDO
                 !
              ENDDO
           ENDIF
           !
           DEALLOCATE(dpsic)
           !
        ENDIF
        !
        ! The perturbation is in aux_r
        !
        ALLOCATE(aux_g(npwx))
        IF(gamma_only) THEN
           CALL SINGLEBAND_fwfft(npwq0,aux_r,dffts%nnr,aux_g,npwx) 
        ELSE
           CALL SINGLEBAND_fwfft_k(npwq0,aux_r,dffts%nnr,aux_g,npwx,.FALSE.) 
        ENDIF
        !
        DO ig=1,npwq0 ! pert acts only on body
           dng(ig,ipert) = dng(ig,ipert) + 2._DP * wk(iks) * aux_g(ig) * sqvc(ig) / omega
        ENDDO
        !
        DEALLOCATE(aux_g)
        DEALLOCATE(aux_r)
        !
        CALL update_bar_type( barra,'dfpt', 1 ) 
        !
     ENDDO ! ipert
     !
     IF( m==0 ) CALL update_bar_type( barra,'dfpt', 1 )
     !
     DEALLOCATE(h_diag)
     DEALLOCATE(dpsi)
     DEALLOCATE(dvpsi)
     !
  ENDDO ! K-POINT and SPIN
  !
  IF( gstart==2 ) dng(1,1:m) = CMPLX( 0._DP, 0._DP, KIND=DP )
  !
  CALL mp_sum(dng,inter_pool_comm)
  !
  CALL mp_barrier( world_comm )
  !
  CALL stop_bar_type( barra, 'dfpt' )
  !
!  CALL close_buffer(iuwfc,'delete')
  !IF( nks > 1 ) THEN
  !   IF ( exst ) THEN
  !      CLOSE(unit=iuwfc,status='keep')
  !   ELSE
  !      CLOSE(unit=iuwfc,status='delete')
  !   ENDIF
  !ENDIF
  !
END SUBROUTINE
!
!
!
SUBROUTINE set_preconditioning(h_diag,m)
  !
  ! Set preconditioning, first eprec, then h_diag
  !
  USE kinds,                 ONLY : DP
  USE wvfct,                 ONLY : g2kin
  USE wavefunctions_module,  ONLY : evc
  USE noncollin_module,      ONLY : noncolin,npol
  USE pwcom,                 ONLY : npw,npwx
  USE mp,                    ONLY : mp_sum
  USE mp_global,             ONLY : intra_bgrp_comm
  !
  ! I/O
  !
  INTEGER,INTENT(IN) :: m
  REAL(DP),INTENT(OUT) :: h_diag(npwx*npol,m)
  !
  ! Workspace
  !
  INTEGER :: ibnd
  COMPLEX(DP),ALLOCATABLE :: aux(:)
  REAL(DP),ALLOCATABLE :: eprec(:)
  COMPLEX(DP),EXTERNAL :: ZDOTC 
  !
  ALLOCATE( aux(npwx*npol) )
  ALLOCATE( eprec(m) )
  !
  h_diag=0._DP
  !
  DO ibnd = 1, m
     !
     ! define eprec
     aux=(0._DP,0._DP)
     DO ig = 1, npw
        aux (ig) = g2kin (ig) * evc (ig, ibnd)
     ENDDO
     IF (noncolin) THEN
        DO ig = 1, npw
           aux (ig+npwx) = g2kin (ig) * evc (ig+npwx, ibnd)
        ENDDO
     ENDIF
     eprec(ibnd) = 1.35_DP * ZDOTC(npw,evc(1,ibnd),1,aux(1),1)
     IF(noncolin) eprec(ibnd) = eprec(ibnd) + 1.35_DP * ZDOTC(npw,evc(1+npwx,ibnd),1,aux(1+npwx),1)
     !
  ENDDO
  !
  CALL mp_sum(eprec,intra_bgrp_comm)
  !
  DO ibnd = 1, m
     !
     ! set h_diag
     DO ig = 1, npw
        h_diag(ig,ibnd)=1._DP/MAX(1._DP,g2kin(ig)/eprec(ibnd))
     ENDDO
     IF (noncolin) THEN
       DO ig = 1, npw
          h_diag(ig+npwx,ibnd)=1._DP/MAX(1._DP,g2kin(ig)/eprec(ibnd))
       ENDDO
     END IF
     !
  ENDDO
  !
  DEALLOCATE( aux )
  DEALLOCATE( eprec )
  !
END SUBROUTINE
