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
SUBROUTINE solve_qp()
  !-----------------------------------------------------------------------
  !
  ! ... This subroutine solves the DBS problem for GAMMA, at non-zero freqeuncies. 
  ! ... Perturbations are distributed according to the POT mpi_communicator
  !
  USE kinds,                ONLY : DP 
  USE westcom,              ONLY : n_pdep_eigen_to_use,n_lanczos,qp_bandrange,iks_l2g,imfreq_list_integrate,&
                                 & n_secant_maxiter,nbnd_occ,trev_secant,io_comm,l_enable_lanczos,imfreq_list,n_imfreq,&
                                 & ihead,ibody1,ibody2,idiago,rhead,rbody,dmat_inv,zmat_inv,l_macropol
  USE mp_global,            ONLY : inter_image_comm,intra_bgrp_comm
  USE mp_world,             ONLY : mpime,root
  USE mp,                   ONLY : mp_barrier,mp_sum
  USE io_global,            ONLY : stdout, ionode
  USE pwcom,                ONLY : npw,npwx,et,nks,current_spin,isk,xk,nbnd,lsda,igk,g2kin,ecutwfc,current_k
  USE io_push,              ONLY : io_push_title,io_push_bar
  USE constants,            ONLY : rytoev
  USE west_io,              ONLY : serial_table_output
  USE distribution_center,  ONLY : pert,ifr,rfr,aband
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  REAL(DP),ALLOCATABLE :: sigma_hf(:,:)
  COMPLEX(DP),ALLOCATABLE :: sigma_cor_in (:,:)
  COMPLEX(DP),ALLOCATABLE :: sigma_cor_out(:,:)
  REAL(DP),ALLOCATABLE :: z_in(:,:)
  REAL(DP),ALLOCATABLE :: qp_energy(:,:)
  COMPLEX(DP),ALLOCATABLE :: sc(:,:,:) 
  REAL(DP),ALLOCATABLE :: en(:,:,:) 
  LOGICAL,ALLOCATABLE :: l_conv(:,:)
  REAL(DP),PARAMETER   :: eshift = 0.007349862_DP ! = 0.1 eV
  INTEGER :: ib, iks, ifixed, ip, glob_ip, ifreq, il, im, glob_im, glob_jp, glob_ifreq
  REAL(DP),ALLOCATABLE :: out_tab(:,:)
  CHARACTER(LEN=5) :: myglobk 
  CHARACTER(LEN=6) :: cifixed
  INTEGER :: notconv
  INTEGER,EXTERNAL :: get_nbndval
  INTEGER :: nbndval 
  REAL(DP),ALLOCATABLE :: dtemp(:)
  COMPLEX(DP),ALLOCATABLE :: ztemp(:)
  REAL(DP),ALLOCATABLE :: diago(:,:)
  REAL(DP),ALLOCATABLE :: braket(:,:,:)
  REAL(DP),ALLOCATABLE :: overlap(:,:)
  TYPE(bar_type) :: barra
  INTEGER :: barra_load
  ! 
  CALL start_clock('solve_qp')
  !
  CALL io_push_title("(Q)uasiparticle corrections")
  !
  ALLOCATE( imfreq_list_integrate( 2, ifr%nloc ) )
  ALLOCATE( dtemp( n_imfreq ) ) 
  !
  ALLOCATE( sigma_hf (qp_bandrange(1):qp_bandrange(2),nks) )
  ALLOCATE( sigma_cor_in (qp_bandrange(1):qp_bandrange(2),nks) )
  ALLOCATE( sigma_cor_out(qp_bandrange(1):qp_bandrange(2),nks) )
  ALLOCATE( z_in (qp_bandrange(1):qp_bandrange(2),nks) )
  ALLOCATE( qp_energy (qp_bandrange(1):qp_bandrange(2),nks) )
  !
  dtemp = 0._DP
  DO ifreq = 1, ifr%nloc 
     glob_ifreq = ifr%l2g(ifreq)
     dtemp( glob_ifreq ) = imfreq_list( ifreq )
  ENDDO
  CALL mp_sum( dtemp, intra_bgrp_comm ) 
  !
  DO ifreq = 1, ifr%nloc 
     glob_ifreq = ifr%l2g(ifreq)
     IF( glob_ifreq == 1 ) THEN
        imfreq_list_integrate(1,ifreq) = 0._DP
     ELSE
        imfreq_list_integrate(1,ifreq) = ( dtemp(glob_ifreq) + dtemp(glob_ifreq-1) ) * 0.5_DP
     ENDIF
     !
     IF( glob_ifreq == n_imfreq ) THEN
        imfreq_list_integrate(2,ifreq) = dtemp(n_imfreq)
     ELSE
        imfreq_list_integrate(2,ifreq) = ( dtemp(glob_ifreq) + dtemp(glob_ifreq+1) ) * 0.5_DP
     ENDIF
     !
  ENDDO
  !
  DEALLOCATE( dtemp )
  !
  ! TEMP
  ! 
  ALLOCATE( ibody1( aband%nloc, ifr%nloc, qp_bandrange(1):qp_bandrange(2), nks ) )
  ALLOCATE( rbody( aband%nloc, rfr%nloc, qp_bandrange(1):qp_bandrange(2), nks ) )
  IF( l_enable_lanczos ) THEN
     ALLOCATE( ibody2( n_lanczos, pert%nloc, ifr%nloc, qp_bandrange(1):qp_bandrange(2), nks ) )
     ALLOCATE( idiago( n_lanczos, pert%nloc, qp_bandrange(1):qp_bandrange(2), nks ) )
  ENDIF
  !
  ibody1 = 0._DP
  rbody = 0._DP
  IF( l_enable_lanczos ) THEN
     ibody2 = 0._DP
     idiago = 0._DP
  ENDIF
  !
  ! ibody1, ibody2, rdiago, idiago
  !
  WRITE(stdout,'(5x,"Collecting results from W and G")') 
  CALL io_push_bar
  !
  barra_load = nks * ( qp_bandrange(2) - qp_bandrange(1) + 1 )
  CALL start_bar_type( barra, 'coll_gw', barra_load )
  !
  ! LOOP 
  !
  DO iks = 1, nks   ! KPOINT-SPIN
     !
     nbndval = nbnd_occ(iks)
     !
     DO ib = qp_bandrange(1), qp_bandrange(2)
        !
        IF(ALLOCATED(overlap)) DEALLOCATE(overlap) 
        ALLOCATE(overlap(pert%nglob, nbnd ) )
        !
        CALL readin_overlap( 'g', iks_l2g(iks), ib, overlap, pert%nglob, nbnd )
        !
        ! ------
        ! ibody1
        ! ------
        !
        ALLOCATE( dtemp(nbnd) )
        !
        DO ifreq = 1, ifr%nloc 
           !
           dtemp = 0._DP
           !
           DO im = 1, nbnd 
              !
              DO glob_jp = 1, n_pdep_eigen_to_use
                 DO ip = 1, pert%nloc
                    glob_ip = pert%l2g(ip)
                    dtemp( im ) = dtemp( im ) + overlap(glob_jp,im)*overlap(glob_ip,im)*dmat_inv(glob_jp,ip,ifreq) 
                 ENDDO
              ENDDO 
              !
           ENDDO ! im
           !
           CALL mp_sum( dtemp, inter_image_comm ) 
           !
           DO im = 1, aband%nloc
              glob_im = aband%l2g(im)
              ibody1(im,ifreq,ib,iks) = dtemp(glob_im)
           ENDDO
           !
        ENDDO ! ifreq
        !
        DEALLOCATE( dtemp ) 
        !
        ! -----
        ! rbody
        ! -----
        !
        ALLOCATE( ztemp(nbnd) )
        !
        DO ifreq = 1, rfr%nloc 
           !
           ztemp = 0._DP
           !
           DO im = 1, nbnd 
              !
              DO glob_jp = 1, n_pdep_eigen_to_use
                 DO ip = 1, pert%nloc
                    glob_ip = pert%l2g(ip)
                    ztemp( im ) = ztemp( im ) + overlap(glob_jp,im)*overlap(glob_ip,im)*zmat_inv(glob_jp,ip,ifreq) 
                 ENDDO
              ENDDO 
              !
           ENDDO ! im
           !
           CALL mp_sum( ztemp, inter_image_comm ) 
           !
           DO im = 1, aband%nloc
              glob_im = aband%l2g(im)
              rbody(im,ifreq,ib,iks) = ztemp(glob_im)
           ENDDO
           !
        ENDDO ! ifreq
        !
        DEALLOCATE( ztemp ) 
        !
        ! -----------------------------
        ! LANCZOS part : idiago, ibody2
        ! -----------------------------
        !
        IF( l_enable_lanczos ) THEN
           !
           ALLOCATE( braket(pert%nglob,n_lanczos,pert%nloc) )
           ALLOCATE( diago(n_lanczos,pert%nloc) )
           CALL readin_solvegfreq( iks_l2g(iks), ib, diago, braket, io_comm, pert%nloc, pert%nglob, pert%myoffset )
           !
           DO ip = 1, pert%nloc
              DO il = 1, n_lanczos
                 idiago(il,ip,ib,iks) = diago(il,ip) 
              ENDDO
           ENDDO
           !
           DO ifreq = 1,ifr%nloc
              DO ip = 1, pert%nloc
                 DO il = 1, n_lanczos
                    DO glob_jp = 1, pert%nglob
                       ibody2(il,ip,ifreq,ib,iks) = ibody2(il,ip,ifreq,ib,iks) + braket(glob_jp,il,ip)*dmat_inv(glob_jp,ip,ifreq)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           !
           DEALLOCATE(diago)
           DEALLOCATE(braket)
           !
        ENDIF
        !
        CALL update_bar_type( barra, 'coll_gw', 1 )
        !
     ENDDO ! ibnd
     !
  ENDDO ! iks
  !
  CALL stop_bar_type( barra, 'coll_gw' )
  !
  !
  DEALLOCATE( dmat_inv )
  DEALLOCATE( zmat_inv )
  !
  ! Get Sigma_X
  !
  CALL readin_solvehf( sigma_hf(qp_bandrange(1),1), qp_bandrange(2)-qp_bandrange(1)+1, nks )
  !
  ! For CORR
  !
  ALLOCATE(en(qp_bandrange(1):qp_bandrange(2),nks,2)) 
  ALLOCATE(sc(qp_bandrange(1):qp_bandrange(2),nks,2)) 
  ALLOCATE(l_conv(qp_bandrange(1):qp_bandrange(2),nks)) 
  !
  ! 1st step of secant solver : E_KS - 0.5 * eshift
  !
  DO iks = 1, nks
     DO ib = qp_bandrange(1), qp_bandrange(2)
        en(ib,iks,1) = et(ib,iks) - eshift*0.5_DP
     ENDDO
  ENDDO
  CALL calc_corr( sc(:,:,1), en(:,:,1) )  
  !
  DO iks = 1, nks
     DO ib = qp_bandrange(1), qp_bandrange(2)
        WRITE(stdout,"(5X,1f14.6,' : ',2f14.6)") en(ib,iks,1) * rytoev, sc(ib,iks,1) * rytoev
     ENDDO
  ENDDO
  !
  ! 1st step of secant solver : E_KS + 0.5 * eshift
  !
  DO iks = 1, nks
     DO ib = qp_bandrange(1), qp_bandrange(2)
        en(ib,iks,2) = et(ib,iks) + eshift*0.5_DP
     ENDDO
  ENDDO
  CALL calc_corr( sc(:,:,2), en(:,:,2) ) 
  DO iks = 1, nks
     DO ib = qp_bandrange(1), qp_bandrange(2)
        WRITE(stdout,"(5X,1f14.6,' : ',2f14.6)") en(ib,iks,2) * rytoev, sc(ib,iks,2) * rytoev
     ENDDO
  ENDDO
  !
  ! Stage sigma_corr_in
  !
  DO iks = 1, nks
     DO ib = qp_bandrange(1), qp_bandrange(2)
        sigma_cor_in(ib,iks) = ( sc(ib,iks,2) + sc(ib,iks,1) ) * 0.5_DP 
     ENDDO
  ENDDO
  !
  ! Stage z_in
  !
  DO iks = 1, nks
     DO ib = qp_bandrange(1), qp_bandrange(2)
        z_in(ib,iks) = 1._DP / ( 1._DP - REAL( sc(ib,iks,2) - sc(ib,iks,1), KIND=DP ) / eshift  ) 
     ENDDO
  ENDDO
  !
  ! en 1 = EKS, sc 1 = sigma_corr_in
  ! en 2 = EKS + Z * ( sigmax - vxc + sigmacorrin)
  !
  DO iks = 1, nks
     DO ib = qp_bandrange(1), qp_bandrange(2)
        en(ib,iks,1) = et(ib,iks) 
        sc(ib,iks,1) = sigma_cor_in(ib,iks)
        en(ib,iks,2) = et(ib,iks) + z_in(ib,iks) * ( sigma_hf(ib,iks) + REAL( sigma_cor_in(ib,iks),KIND=DP) ) 
     ENDDO
  ENDDO
  CALL output_eqp_report(0,en(:,:,1),en(:,:,2),sigma_cor_in(:,:))
  !
  WRITE(stdout,"(/,5X)") 
  CALL io_push_bar()
  WRITE(stdout,"(5X,'Iter: ',i6.6,', QP energies [eV]')") 0
  CALL io_push_bar()
  WRITE(stdout,"(5X,6f14.6)") en(:,:,2) * rytoev
  CALL io_push_bar()
  !
  ! nth step of the secan solver
  !
  l_conv = .FALSE.
  notconv = nks * ( qp_bandrange(2) - qp_bandrange(1) + 1 )
  DO ifixed = 1, n_secant_maxiter
     ! 
     WRITE(cifixed,"(i6.6)") ifixed
     CALL io_push_title("Iter : "//cifixed)
     !
     CALL calc_corr( sc(:,:,2), en(:,:,2) )  
     DO iks = 1, nks
        DO ib = qp_bandrange(1), qp_bandrange(2)
           WRITE(stdout,"(5X,1f14.6,' : ',2f14.6)") en(ib,iks,2) * rytoev, sc(ib,iks,2) * rytoev
        ENDDO
     ENDDO
     !
     ! Resulting new energy:
     !
     DO iks = 1, nks
        DO ib = qp_bandrange(1), qp_bandrange(2)
            IF( .NOT. l_conv(ib,iks) ) THEN
               qp_energy(ib,iks) = en(ib,iks,2) + &
                         & ( et(ib,iks) + sigma_hf(ib,iks) + REAL(sc(ib,iks,2),KIND=DP) - en(ib,iks,2) ) / &
                         & ( 1._DP - REAL( sc(ib,iks,2) - sc(ib,iks,1), KIND=DP ) / ( en(ib,iks,2) - en(ib,iks,1) ) )
            ENDIF
        ENDDO
     ENDDO
     !
     WRITE(stdout,"(/,5X)") 
     CALL io_push_bar()
     WRITE(stdout,"(5X,'Iter: ',i6.6,', QP energies [eV]')") ifixed
     CALL io_push_bar()
     WRITE(stdout,"(5X,6f14.6)") qp_energy(:,:) * rytoev
     CALL io_push_bar()
     !
     ! Estimate l_conv
     !
     DO iks = 1, nks
        DO ib = qp_bandrange(1), qp_bandrange(2)
           l_conv(ib,iks) = ( ABS( qp_energy(ib,iks) - en(ib,iks,2) ) < trev_secant ) 
        ENDDO
     ENDDO
     !
     ! Count the number of notconverged QP energies
     !
     notconv = 0 
     DO iks = 1, nks
        DO ib = qp_bandrange(1), qp_bandrange(2)
           IF( .NOT.l_conv(ib,iks) ) notconv = notconv + 1
        ENDDO
     ENDDO
     !
     IF( notconv == 0 ) THEN 
        ! 
        CALL output_eqp_report(-1,en(:,:,2),qp_energy,sc(:,:,2))
        CALL io_push_title("CONVERGENCE ACHIEVED !!!")
        EXIT
     ELSE
        CALL io_push_bar()
        WRITE(stdout,'(5x,"Number of unconverged QP. = ",i12)') notconv
        CALL io_push_bar()
        CALL output_eqp_report(ifixed,en(:,:,2),qp_energy,sc(:,:,2))
     ENDIF
     !
     ! Iterate
     !
     en(:,:,1) = en(:,:,2)    
     sc(:,:,1) = sc(:,:,2) 
     en(:,:,2) = qp_energy(:,:)   
     !
  ENDDO
  !
  sigma_cor_out(:,:) = sc(:,:,2)
  DEALLOCATE( en, sc, l_conv )
  !
  IF( notconv .NE. 0 ) THEN 
     CALL io_push_title("CONVERGENCE **NOT** ACHIEVED !!!")
  ENDIF
  !
  ! Output it per k-point
  !
  ALLOCATE(out_tab(qp_bandrange(2)-qp_bandrange(1)+1,7))
  !
  DO iks=1,nks
     DO ib = qp_bandrange(1), qp_bandrange(2)
        out_tab( ib - qp_bandrange(1) + 1, 1) = REAL( ib, KIND=DP) 
        out_tab( ib - qp_bandrange(1) + 1, 2) = et(ib,iks) * rytoev
        out_tab( ib - qp_bandrange(1) + 1, 3) = (et(ib,iks)+sigma_hf(ib,iks)) * rytoev
        out_tab( ib - qp_bandrange(1) + 1, 4) = qp_energy(ib,iks) * rytoev
        out_tab( ib - qp_bandrange(1) + 1, 5) = (qp_energy(ib,iks) - et(ib,iks) ) * rytoev
        out_tab( ib - qp_bandrange(1) + 1, 6) = REAL(  sigma_cor_out(ib,iks), KIND=DP ) * rytoev
        out_tab( ib - qp_bandrange(1) + 1, 7) = AIMAG( sigma_cor_out(ib,iks) ) * rytoev
     ENDDO
     WRITE(myglobk,'(i5.5)') iks_l2g(iks)
     !
     CALL serial_table_output(mpime==root,4000,'eqp_K'//myglobk,out_tab,&
        & qp_bandrange(2)-qp_bandrange(1)+1,7,&
        & (/'      band','    E0[eV]','   EHF[eV]','   Eqp[eV]','Eqp-E0[eV]','Sc_Eqp[eV]',' Width[eV]'/))
  ENDDO
  !
  DEALLOCATE( out_tab )  
  !
  DEALLOCATE( sigma_hf )
  DEALLOCATE( sigma_cor_in )
  DEALLOCATE( sigma_cor_out )
  DEALLOCATE( z_in )
  DEALLOCATE( qp_energy )
  !
  CALL stop_clock( "solve_qp" )
  !
  CALL io_push_title('Done, take a look at the o-eqp_K*.tab file(s) .')
  !
  !
END SUBROUTINE
!
!
SUBROUTINE output_eqp_report(iteration,en1,en2,sc1)
  !
  USE kinds,                ONLY : DP
  USE westcom,              ONLY : qp_bandrange
  USE pwcom,                ONLY : nks,et
  USE constants,            ONLY : rytoev
  USE west_io,              ONLY : serial_table_output
  USE mp_world,             ONLY : mpime,root
  ! 
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER,INTENT(IN) :: iteration
  REAL(DP),INTENT(IN) :: en1(qp_bandrange(1):qp_bandrange(2),nks)
  REAL(DP),INTENT(IN) :: en2(qp_bandrange(1):qp_bandrange(2),nks)
  COMPLEX(DP),INTENT(IN) :: sc1(qp_bandrange(1):qp_bandrange(2),nks)
  !
  ! Workspace
  !
  CHARACTER(LEN=9) :: prefisso
  INTEGER :: contatore
  REAL(DP) :: out_tabella(nks*(qp_bandrange(2)-qp_bandrange(1)+1),7)
  INTEGER :: ib, iks
  !
  contatore=0
  DO iks=1,nks
     DO ib = qp_bandrange(1), qp_bandrange(2)
        contatore = contatore + 1
        out_tabella(contatore,1) = REAL(iks,DP)
        out_tabella(contatore,2) = REAL(ib,DP)
        out_tabella(contatore,3) = et(ib,iks)*rytoev
        out_tabella(contatore,4) = en1(ib,iks)*rytoev
        out_tabella(contatore,5) = REAL( sc1(ib,iks), KIND=DP) *rytoev
        out_tabella(contatore,6) = en2(ib,iks)*rytoev
        out_tabella(contatore,7) = (en2(ib,iks)-en1(ib,iks))*rytoev
     ENDDO
  ENDDO
  !
  IF(iteration>=0) THEN 
     WRITE(prefisso,"('itr_',i5.5)") iteration
  ELSE
     prefisso='converged'
  ENDIF
  CALL serial_table_output(mpime==root,4000,'eqp.'//TRIM(ADJUSTL(prefisso)),out_tabella,&
  & nks*(qp_bandrange(2)-qp_bandrange(1)+1),7,&
  & (/'       iks','        ib','   Eks[eV]','   Ein[eV]','Sc_Ein[eV]','  Eout[eV]',' Diff.[eV]'/))
  !
END SUBROUTINE
