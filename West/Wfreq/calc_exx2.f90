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
SUBROUTINE calc_exx2( sigma_exx, nb1, nb2 )
  !-----------------------------------------------------------------------
  !
  ! store in sigma_exx(n,iks) = < n,iks | V_exx | n,iks >     n = nb1, nb2
  !
  USE kinds,                ONLY : DP 
  USE mp_global,            ONLY : inter_image_comm,intra_bgrp_comm,my_image_id
  USE mp,                   ONLY : mp_sum,mp_bcast
  USE io_global,            ONLY : stdout, ionode
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE scf,                  ONLY : rho, rho_core, rhog_core
  USE gvect,                ONLY : g,nl,gstart,ngm_g,ig_l2g,ngm
  USE gvecs,                ONLY : ngms
  USE cell_base,            ONLY : tpiba2,omega,tpiba,at,alat
  USE fft_base,             ONLY : dfftp,dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE constants,            ONLY : tpi,fpi,rytoev,e2
  USE pwcom,                ONLY : npw,npwx,et,nks,current_spin,isk,xk,nbnd,lsda,igk,g2kin,ecutwfc,nspin,current_k
  USE fft_at_gamma,         ONLY : DOUBLEBAND_INVFFT,SINGLEBAND_INVFFT,DOUBLEBAND_FWFFT,SINGLEBAND_FWFFT,smooth_fwfft
  USE fft_at_k,             ONLY : SINGLEBAND_INVFFT_k,SINGLEBAND_FWFFT_k
  USE wavefunctions_module, ONLY : evc,psic,psic_nc
  USE westcom,              ONLY : iuwfc,lrwfc,npwq0,nbnd_occ,div_kind_hf
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin,npol 
  USE buffers,              ONLY : get_buffer
  USE funct,                ONLY : dft_is_hybrid,init_dft_exxrpa,stop_exx
  USE exx,                  ONLY : x_gamma_extrapolation,exxdiv_treatment,exx_grid_init,exx_div_check,&
                                   &deallocate_exx,exxinit,vexx,exx_grid_initialized
  USE uspp,                 ONLY : vkb,nkb
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE io_push,              ONLY : io_push_bar
  USE distribution_center,  ONLY : index_distribution,destroy_container,index_distr_init
  USE coulomb_vcut_module,  ONLY : vcut_init, vcut_type, vcut_info, &
                                   vcut_get,  vcut_spheric_get, vcut_destroy
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: nb1, nb2
  REAL(DP),INTENT(OUT) :: sigma_exx( nb1:nb2, nks) 
  !
  ! Workspace
  !
  COMPLEX(DP),ALLOCATABLE :: pertg(:),pertr(:)
  REAL(DP), EXTERNAL :: DDOT
  INTEGER :: ib,iv,i1,ir,iks,ig,iv_glob
  INTEGER :: nbndval
  TYPE(index_distribution) :: vband
  REAL(DP) :: peso
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  INTEGER :: nq1, nq2, nq3         ! integers defining the X integration mesh
  REAL(DP) :: atws(3,3)
  REAL(DP),ALLOCATABLE :: mysqvc(:)
  REAL(DP) :: q(3)
  REAL(DP) :: ecutvcut
  TYPE(vcut_type)   :: vcut
  REAL(DP) :: mydiv
  !
  WRITE(stdout,'(5x,a)') ' '
  CALL io_push_bar()
  WRITE(stdout,'(5x,a)') 'Sigma_X'
  CALL io_push_bar()
  !
  ALLOCATE( pertg( ngms ) )
  ALLOCATE( pertr( dffts%nnr ) )
  ALLOCATE( mysqvc(ngms) )
  !
  CALL store_sqvc(mysqvc,ngms,div_kind_hf,mydiv)
  !
  ! Set to zero
  !
  sigma_exx = 0._DP
  !
  barra_load = nks * ( nb2-nb1 + 1 )
  CALL start_bar_type( barra, 'sigmax', barra_load )
  !
  ! LOOP 
  !
  DO iks = 1, nks   ! KPOINT-SPIN
     !
     IF( gamma_only ) THEN 
        peso = 2._DP  
     ELSE
        peso = 1._DP
     ENDIF
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
     CALL index_distr_init(vband,nbndval,'i','nbndval',.FALSE.)
     !
     DO ib = nb1, nb2
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
        DO iv = 1, vband%nloc
           !
           iv_glob = vband%l2g(iv)
           !
           ! Bring it to R-space
!           IF(gamma_only) THEN
              CALL SINGLEBAND_invfft(npw,evc(1,iv_glob),npwx,pertr,dffts%nnr)
              DO ir=1,dffts%nnr 
                 pertr(ir)=psic(ir)*pertr(ir)
              ENDDO
              !CALL SINGLEBAND_fwfft(npw,pertr,dffts%nnr,pertg,npwx)
              CALL smooth_fwfft(ngms,pertr,dffts%nnr,pertg,ngms)
!           ELSEIF(noncolin) THEN
!              CALL SINGLEBAND_invfft_k(npwq0,evc(1,iv_glob),npwx,pertr,dffts%nnr,.FALSE.)
!              DO ir=1,dffts%nnr 
!                 pertr(ir)=psic_nc(ir,1)*pertr(ir)
!              ENDDO
!              CALL SINGLEBAND_fwfft_k(npw,pertr,dffts%nnr,pertg(1),npwx,.TRUE.)
!              CALL SINGLEBAND_invfft_k(npwq0,evc(1+npwx,iv_glob),npwx,pertr,dffts%nnr,.FALSE.)
!              DO ir=1,dffts%nnr 
!                 pertr(ir)=psic_nc(ir,2)*pertr(ir)
!              ENDDO
!              CALL SINGLEBAND_fwfft_k(npw,pertr,dffts%nnr,pertg(1+npwx),npwx,.TRUE.)
!           ELSE
!              CALL SINGLEBAND_invfft_k(npwq0,evc(1,iv_glob),npwx,pertr,dffts%nnr,.FALSE.)
!              DO ir=1,dffts%nnr 
!                 pertr(ir)=psic(ir)*pertr(ir)
!              ENDDO
!              CALL SINGLEBAND_fwfft_k(npw,pertr,dffts%nnr,pertg,npwx,.TRUE.)
!           ENDIF 
           !
           DO ig = 1,ngms
              pertg(ig) = pertg(ig) * mysqvc(ig) 
           ENDDO
           sigma_exx( ib, iks ) = sigma_exx( ib, iks ) - peso * DDOT( 2*ngms, pertg(1), 1, pertg(1), 1) / omega
           !IF(gstart==2) sigma_exx( ib, iks ) = sigma_exx( ib, iks ) + REAL( pertg(1), KIND = DP )**2 / omega
           IF( ib == iv_glob .AND. gstart == 2 ) sigma_exx( ib, iks ) = sigma_exx( ib, iks ) - mydiv
           !
!           IF( noncolin ) THEN
!              DO ig = 1,npwq0
!                 pertg(ig+npwx) = pertg(ig+npwx) * mysqvc(ig) 
!              ENDDO
!              sigma_exx( ib, iks ) = sigma_exx( ib, iks ) - DDOT( 2*npwq0, pertg(1+npwx), 1, pertg(1+npwx), 1) / omega
!              IF( ib == iv_glob .AND. gstart == 2 ) sigma_exx( ib, iks ) = sigma_exx( ib, iks ) - mydiv 
!           ENDIF 
           !
        ENDDO
        !
        CALL update_bar_type( barra, 'sigmax', 1 )
        !
     ENDDO ! ib
     !
     CALL destroy_container(vband)
     !
  ENDDO ! iks
  !
  CALL stop_bar_type( barra, 'sigmax' )
  !
  CALL mp_sum( sigma_exx, intra_bgrp_comm )
  CALL mp_sum( sigma_exx, inter_image_comm ) 
  !
  DEALLOCATE( pertr, pertg, mysqvc )
  !
END SUBROUTINE 
