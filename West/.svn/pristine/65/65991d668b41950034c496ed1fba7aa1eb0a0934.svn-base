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
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
!----------------------------------------------------------------------------
SUBROUTINE davidson_diago ( )
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problem:
  ! ... ( chi - ev ) * dvg = 0
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE distribution_center,  ONLY : index_distr_init,pert,destroy_container
  USE pwcom,                ONLY : npwx  
  USE io_push,              ONLY : io_push_title,io_push_bar
  USE westcom,              ONLY : dvg,dng,n_pdep_eigen,trev_pdep,n_pdep_maxiter,n_pdep_basis,wstat_calculation,ev,conv,&
                                   & n_pdep_restart_from_itr,n_pdep_read_from_file,n_steps_write_restart,n_pdep_times 
  USE pdep_db,              ONLY : pdep_db_write,pdep_db_read
  USE wstat_restart,        ONLY : wstat_restart_write, wstat_restart_clear, wstat_restart_read
  USE mp_world,             ONLY : mpime
  USE mp_global,            ONLY : inter_image_comm
  USE mp,                   ONLY : mp_sum
  USE gvect,                ONLY : gstart
  !
  IMPLICIT NONE
  !
  ! ... LOCAL variables
  !
  INTEGER :: nvec, nvecx
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
  INTEGER :: dav_iter, notcnv
    ! integer  number of iterations performed
    ! number of unconverged roots
  INTEGER :: kter, nbase, np, n, m
    ! counter on iterations
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! do-loop counters
    ! counter on the bands
  INTEGER :: ierr,mloc,mstart,mend
  INTEGER,ALLOCATABLE :: ishift(:)
  REAL(DP), ALLOCATABLE :: ew(:), hr_distr(:,:), vr_distr(:,:)
    ! Hamiltonian on the reduced basis
    ! S matrix on the reduced basis
    ! eigenvectors of the Hamiltonian
    ! eigenvalues of the reduced hamiltonian
  INTEGER :: il1,il2,ig1,ig2,i
  REAL(DP) :: time_spent(2)
  CHARACTER(LEN=8) :: iter_label
  !
  REAL(DP), EXTERNAL :: GET_CLOCK
  !
  ! ... INITIALIZATION
  !
  nvec = n_pdep_eigen
  nvecx = n_pdep_basis
  !
  CALL start_clock( 'chidiago' )
  time_spent(1)=get_clock( 'chidiago' )
  !
  ! ... DISTRIBUTE nvecx
  !
  CALL index_distr_init(pert,nvecx,'i','nvecx',.TRUE.) 
  CALL wstat_memory_report() ! Before allocating I report the memory required. 
  !
  ! ... MEMORY ALLOCATION
  !
  IF ( nvec > nvecx / 2 ) CALL errore( 'chidiago', 'nvecx is too small', 1 )
  !
  ALLOCATE( dvg( npwx, pert%nlocx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'chidiago',' cannot allocate dvg ', ABS(ierr) )
  !
  ALLOCATE( dng( npwx, pert%nlocx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'chidiago',' cannot allocate dng ', ABS(ierr) )
  !
  ALLOCATE( hr_distr( nvecx, pert%nlocx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'chidiago',' cannot allocate hr_distr ', ABS(ierr) )
  !
  ALLOCATE( vr_distr( nvecx, pert%nlocx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'chidiago',' cannot allocate vr_distr ', ABS(ierr) )
  !
  ALLOCATE( ew( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'chidiago',' cannot allocate ew ', ABS(ierr) )
  !
  ALLOCATE( ev( nvec  ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'chidiago',' cannot allocate ev ', ABS(ierr) )
  !
  ALLOCATE( conv( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'chidiago',' cannot allocate conv ', ABS(ierr) )
  !
  nbase  = nvec
  conv   = .FALSE.
  ev     = 0._DP
  ew     = 0._DP
  dng  = 0._DP
  dvg  = 0._DP
  hr_distr(:,:) = 0._DP
  vr_distr(:,:) = 0._DP
  notcnv  = nvec 
  dav_iter = -2
  !
  ! KIND OF CALCULATION
  !
  SELECT CASE(wstat_calculation)
  CASE('r','R')
     !
     ! RESTART
     !
     CALL wstat_restart_read( dav_iter, notcnv, nbase, ew, hr_distr, vr_distr )
     !
  CASE('s','S')
     !
     ! FROM SCRATCH
     !
     ! ... Eventually read from file
     !
     IF(n_pdep_read_from_file>0) CALL pdep_db_read( n_pdep_read_from_file )
     !
     ! ... Eventually randomize
     !
     IF(n_pdep_read_from_file<nvec) CALL do_randomize ( dvg, n_pdep_read_from_file+1, nvec  )
     !
     ! ... MGS
     !
     CALL do_mgs( dvg, 1, nvec)
     !
     WRITE(stdout, "( /,5x,'                  *----------*              *----------*               *----------*') ")
     WRITE(stdout, &
         & "(   5x,'#     Iteration = | ', a8,' |','   ','DFPT_dim = | ', i8,' |', '   ','Diago_dim = | ', i8,' |  x 1/2')")&
         & 'starting', nbase, nbase
     WRITE(stdout, "(   5x,'                  *----------*              *----------*               *----------*') ")
     !
     ! Apply operator with DFPT
     !
     mloc = 0
     mstart = 1 
     DO il1 = 1, pert%nloc
        ig1 = pert%l2g(il1)
        IF( ig1 < 1 .OR. ig1 > nvec ) CYCLE
        IF( mloc==0 ) mstart = il1
        mloc = mloc + 1
     ENDDO
     !
     CALL dfpt ( mloc, dvg(1,mstart), dng(1,mstart) ) 
     dav_iter = -1
     CALL wstat_restart_write( dav_iter, notcnv, nbase, ew, hr_distr, vr_distr)
     !
  CASE DEFAULT
     CALL errore('chidiago', 'Wrong wstat_calculation',1)
  END SELECT
  !
  IF( dav_iter == -2 ) CALL errore( 'chidiago','Cannot find the 1st starting loop',1) 
  !
  IF( dav_iter == -1 ) THEN
     !
     ! < EXTRA STEP > 
     !
     dvg = dng
     CALL do_mgs( dvg, 1, nvec)
     !
     WRITE(stdout, "( /,5x,'                  *----------*              *----------*               *----------*') ")
     WRITE(stdout, &
         & "(   5x,'#     Iteration = | ', a8,' |','   ','DFPT_dim = | ', i8,' |', '   ','Diago_dim = | ', i8,' |  x 2/2')")&
         & 'starting', nbase, nbase
     WRITE(stdout, "(   5x,'                  *----------*              *----------*               *----------*') ")
     !
     ! Apply operator with DFPT
     !
     mloc = 0
     mstart = 1 
     DO il1 = 1, pert%nloc
        ig1 = pert%l2g(il1)
        IF( ig1 < 1 .OR. ig1 > nvec ) CYCLE
        IF( mloc==0 ) mstart = il1
        mloc = mloc + 1
     ENDDO
     !
     CALL dfpt ( mloc, dvg(1,mstart), dng(1,mstart) ) 
     ! 
     ! </ EXTRA STEP >
     !
     ! hr = <dvg|dng>
     !
     CALL build_hr( dvg, dng, mstart, mstart+mloc-1, hr_distr, 1, nvec ) 
     !
     ! ... diagonalize the reduced hamiltonian
     !
     CALL diagox( nbase, nvec, hr_distr, nvecx, ew, vr_distr )
     time_spent(2)=get_clock( 'chidiago' )
     ev(1:nvec) = ew(1:nvec)
     !
     ! Write the eigenvalues & time spent
     !
     CALL output_ev_and_time(nvec,ev,time_spent)
     dav_iter = 0
     CALL wstat_restart_write( dav_iter, notcnv, nbase, ew, hr_distr, vr_distr)
     !
  ENDIF
  !
  ! --------------------
  !
  ! ... iterate
  !
  ! --------------------
  !
  iterate: DO kter = 1, n_pdep_maxiter
     !
     time_spent(1) = get_clock( 'chidiago' )
     !
     dav_iter = dav_iter + 1
     !
     WRITE(stdout, "( /,5x,'                  *----------*              *----------*               *----------*') ")
     WRITE(stdout, "(   5x,'#     Iteration = | ', i8,' |','   ','DFPT_dim = | ', i8,' |', '   ','Diago_dim = | ', i8,' |')") &
         &dav_iter, notcnv, nbase+notcnv
     WRITE(stdout, "(   5x,'                  *----------*              *----------*               *----------*') ")
     !
     ALLOCATE( ishift( nvecx ), STAT=ierr )
     IF( ierr /= 0 ) CALL errore( 'chidiago',' cannot allocate ishift ', ABS(ierr) )
     ishift=0
     np = 0
     !
     DO n = 1, nvec
        !
        IF ( .NOT. conv(n) ) THEN
           !
           ! ... this root not yet converged ... 
           !
           np = np + 1
           !
           ! ... reorder eigenvectors so that coefficients for unconverged
           ! ... roots come first. This allows to use quick matrix-matrix 
           ! ... multiplications to set a new basis vector (see below)
           !
           !IF ( np /= n ) vr(:,np) = vr(:,n)
           ishift(nbase+np) = n
           !
           ew(nbase+np) = ev(n)
           !   
        END IF
        !
     END DO
     !
     ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
     !
     CALL redistribute_vr_distr( notcnv, nbase, nvecx, vr_distr, ishift )
     DEALLOCATE(ishift)
     CALL update_with_vr_distr( dvg, dng, notcnv, nbase, nvecx, vr_distr, ew )
     !
     ! ... MGS
     !
     CALL do_mgs(dvg,nbase+1,nbase+notcnv)
     !
     ! apply the response function to new vectors
     !
     ! determine image that actually compute dng first
     !
     mloc = 0 
     mstart = 1
     DO il1 = 1, pert%nloc
        ig1 = pert%l2g(il1)
        IF( ig1 < nbase+1 .OR. ig1 > nbase+notcnv ) CYCLE
        IF( mloc==0 ) mstart = il1
        mloc = mloc + 1
     ENDDO
     !
     ! Apply operator with DFPT
     !
     CALL dfpt ( mloc, dvg(1,mstart), dng(1,mstart) ) 
     !
     ! ... update the reduced hamiltonian
     !
     !
     ! hr = <dvg|dng>
     !
     CALL build_hr( dvg, dng, mstart, mstart+mloc-1, hr_distr, nbase+1, nbase+notcnv ) 
     !
     nbase = nbase + notcnv
     !
     CALL symm_hr_distr(hr_distr,nbase,nvecx)
     !
     ! ... diagonalize the reduced hamiltonian
     !
     CALL diagox( nbase, nvec, hr_distr, nvecx, ew, vr_distr )
     time_spent(2)=get_clock( 'chidiago' )
     !
     ! ... test for convergence
     !
     conv(1:nvec) = ( ( ABS( ew(1:nvec) - ev(1:nvec) ) < trev_pdep ) )
     !
     notcnv = COUNT( .NOT. conv(:) )
     !
     ev(1:nvec) = ew(1:nvec)
     !
     ! Write the eigenvalues & time spent
     !
     CALL output_ev_and_time(nvec,ev,time_spent)
     !
     ! ... if overall convergence has been achieved, or the dimension of
     ! ... the reduced basis set is becoming too large, or in any case if
     ! ... we are at the last iteration refresh the basis set. i.e. replace
     ! ... the first nvec elements with the current estimate of the
     ! ... eigenvectors;  set the basis dimension to nvec.
     !
     IF ( notcnv == 0 .OR. nbase+notcnv > nvecx .OR. kter == n_pdep_maxiter ) THEN
        !
        CALL start_clock( 'chidiago:last' )
        !
        IF ( notcnv == 0 ) THEN
           !
           ! ... all roots converged: return
           !
           CALL stop_clock( 'chidiago:last' )
           !
           CALL refresh_with_vr_distr( dvg, nvec, nbase, nvecx, vr_distr )
           !
           CALL pdep_db_write( )
           CALL wstat_restart_clear()
           CALL output_a_report(-1)
           !
           WRITE(iter_label,'(i8)') kter
           CALL io_push_title("Convergence achieved !!! in "//TRIM(iter_label)//" steps") 
           !
           EXIT iterate
           !
        ELSE IF ( kter == n_pdep_maxiter ) THEN
           !
           ! ... last iteration, some roots not converged: return
           !
           WRITE( stdout, '(5X,"WARNING: ",I5, &
                &   " eigenvalues not converged in chidiago")' ) notcnv
           !
           CALL stop_clock( 'chidiago:last' )
           !
           EXIT iterate
           !
        END IF
        !
        ! ... refresh psi, H*psi and S*psi
        !
        WRITE(stdout,'(/,7x,"Refresh the basis set")')
        !
        CALL refresh_with_vr_distr( dvg, nvec, nbase, nvecx, vr_distr )
        CALL refresh_with_vr_distr( dng, nvec, nbase, nvecx, vr_distr )
        !
        ! ... refresh the reduced hamiltonian
        !
        nbase = nvec
        !
        hr_distr = 0._DP
        vr_distr = 0._DP
        !
        DO il1 = 1, pert%nloc
           ig1 = pert%l2g(il1)
           hr_distr(ig1,il1) = ev(ig1)
           vr_distr(ig1,il1) = 1._DP
        ENDDO
        !
        CALL stop_clock( 'chidiago:last' )
        !
     END IF
     !
     CALL wstat_restart_write( dav_iter, notcnv, nbase, ew, hr_distr, vr_distr)
     CALL output_a_report(dav_iter)
     !
  END DO iterate
  !
  DEALLOCATE( conv )
  DEALLOCATE( ew )
  DEALLOCATE( ev )
  DEALLOCATE( hr_distr )
  DEALLOCATE( vr_distr )
  !
  !
  DEALLOCATE( dng )
  DEALLOCATE( dvg )  
  !
  CALL stop_clock( 'chidiago' )
  !
  RETURN
  !
END SUBROUTINE
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
SUBROUTINE update_with_vr_distr( ag, bg, nselect, n, lda, vr_distr, ew )
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : my_image_id,inter_image_comm,intra_bgrp_comm,nimage,world_comm
  USE mp,                   ONLY : mp_sum,mp_circular_shift_left,mp_barrier
  USE distribution_center,  ONLY : pert
  USE westcom,              ONLY : npwq0
  USE pwcom,                ONLY : npwx
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP) :: ag(npwx,pert%nlocx)
  COMPLEX(DP) :: bg(npwx,pert%nlocx)
  INTEGER,INTENT(IN) :: nselect, n, lda
  REAL(DP),INTENT(IN) :: vr_distr(lda,pert%nlocx), ew(lda)
  !
  ! Workspace
  !
  COMPLEX(DP),ALLOCATABLE :: hg(:,:)
  INTEGER,ALLOCATABLE :: tmp_l2g(:)
  INTEGER :: il1, il2, ig1, ig2, icycl
  COMPLEX(DP) :: zconst
  !
  CALL mp_barrier( world_comm )
  !  
  CALL start_clock( 'update_vr' )
  !
  ALLOCATE( hg(npwx,pert%nlocx) )
  hg = 0._DP 
  !
  ALLOCATE( tmp_l2g(1:pert%nlocx) )
  !
  tmp_l2g = 0
  !
  DO il1 = 1, pert%nloc 
     tmp_l2g(il1) = pert%l2g(il1)
  ENDDO
  !
  DO icycl=0,nimage-1
     !
     DO il1=1,pert%nlocx
        !
        ig1 = tmp_l2g(il1)
        IF( ig1 == 0 .OR. ig1 > n ) CYCLE
        !
        DO il2=1,pert%nloc
           !
           ig2 = pert%l2g(il2)
           IF( ig2 <= n .OR. ig2 > n+nselect ) CYCLE 
           !
           zconst = CMPLX( vr_distr(ig1,il2), 0_DP, KIND=DP )
           CALL ZAXPY(npwq0,zconst,ag(1,il1),1,hg(1,il2),1)
           !dhg(:,il2) = dhg(:,il2) + amat(:,il1) * z(ig1,ig2-n) 
           !
        ENDDO
     ENDDO
     !
     ! Cycle the amat array
     ! 
     CALL mp_circular_shift_left( ag,      icycl,        inter_image_comm)
     CALL mp_circular_shift_left( tmp_l2g, icycl+nimage, inter_image_comm)
     !
  ENDDO
  !
  DO il2=1,pert%nloc
     ig2 = pert%l2g(il2)
     IF( ig2 <= n .OR. ig2 > n+nselect ) CYCLE
     ag(:,il2) = - ew(ig2) * hg(:,il2)
  ENDDO
  !
  hg = 0._DP
  !
  DO icycl=0,nimage-1
     !
     DO il1=1,pert%nlocx
        !
        ig1 = tmp_l2g(il1)
        IF( ig1 == 0 .OR. ig1 > n ) CYCLE
        !
        DO il2=1,pert%nloc
           !
           ig2 = pert%l2g(il2)
           IF( ig2 <= n .OR. ig2 > n+nselect ) CYCLE 
           !
           zconst = CMPLX( vr_distr(ig1,il2), 0_DP, KIND=DP )
           CALL ZAXPY(npwq0,zconst,bg(1,il1),1,hg(1,il2),1)
           !dhg(:,il2) = dhg(:,il2) + bmat(:,il1) * z(ig1,ig2-n) 
           !
        ENDDO
     ENDDO
     !
     ! Cycle the bmat array 
     ! 
     CALL mp_circular_shift_left( bg,      icycl,        inter_image_comm)
     CALL mp_circular_shift_left( tmp_l2g, icycl+nimage, inter_image_comm)
     !
  ENDDO
  !
  DO il2=1,pert%nloc
     ig2 = pert%l2g(il2)
     IF( ig2 <= n .OR. ig2 > n+nselect ) CYCLE
     ag(:,il2) = ag(:,il2) + hg(:,il2)
  ENDDO
  !
  DEALLOCATE( tmp_l2g )
  DEALLOCATE( hg )
  !
  CALL stop_clock( "update_vr" )
  !
END SUBROUTINE
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
SUBROUTINE refresh_with_vr_distr( ag, nselect, n, lda, vr_distr )
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : my_image_id,inter_image_comm,intra_bgrp_comm,nimage,world_comm
  USE mp,                   ONLY : mp_sum,mp_circular_shift_left,mp_barrier
  USE distribution_center,  ONLY : pert
  USE westcom,              ONLY : npwq0
  USE pwcom,                ONLY : npwx
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP) :: ag(npwx,pert%nlocx)
  INTEGER,INTENT(IN) :: nselect, n, lda
  REAL(DP),INTENT(IN) :: vr_distr(lda,pert%nlocx)
  !
  ! Workspace
  !
  COMPLEX(DP),ALLOCATABLE :: hg(:,:)
  INTEGER,ALLOCATABLE :: tmp_l2g(:)
  INTEGER :: il1, il2, ig1, ig2, icycl
  COMPLEX(DP) :: zconst
  !
  CALL mp_barrier( world_comm )
  !  
  CALL start_clock( 'refresh_vr' )
  !
  ALLOCATE( hg(npwx,pert%nlocx) )
  hg = 0._DP 
  ALLOCATE( tmp_l2g(1:pert%nlocx) )
  !
  tmp_l2g = 0
  !
  DO il1 = 1, pert%nloc 
     tmp_l2g(il1) = pert%l2g(il1)
  ENDDO
  !
  DO icycl=0,nimage-1
     !
     DO il1=1,pert%nlocx
        !
        ig1 = tmp_l2g(il1)
        IF( ig1 == 0 .OR. ig1 > n ) CYCLE
        !
        DO il2=1,pert%nloc
           !
           ig2 = pert%l2g(il2)
           IF( ig2 > nselect ) CYCLE 
           !
           zconst = CMPLX( vr_distr(ig1,il2), 0_DP, KIND=DP )
           CALL ZAXPY(npwq0,zconst,ag(1,il1),1,hg(1,il2),1)
           !dhg(:,il2) = dhg(:,il2) + amat(:,il1) * z(ig1,ig2) 
           !
        ENDDO
     ENDDO
     !
     ! Cycle the amat array 
     ! 
     CALL mp_circular_shift_left( ag  ,           icycl, inter_image_comm)
     CALL mp_circular_shift_left( tmp_l2g, icycl+nimage, inter_image_comm)
     !
  ENDDO
  !
  DO il2=1,pert%nloc
     ig2 = pert%l2g(il2)
     IF( ig2 > nselect ) THEN
        ag(:,il2) = 0._DP
     ELSE
        ag(:,il2) = hg(:,il2)
     ENDIF
  ENDDO
  !
  DEALLOCATE( hg )
  DEALLOCATE( tmp_l2g )
  !
  CALL stop_clock( "refresh_vr" )
  !
END SUBROUTINE
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
SUBROUTINE build_hr( ag, bg, l2_s, l2_e, c_distr, g_s, g_e )
  !
  !  c_distr = < ag | bg >
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : my_image_id,inter_image_comm,intra_bgrp_comm,nimage,world_comm
  USE mp,                   ONLY : mp_sum,mp_circular_shift_left,mp_barrier
  USE distribution_center,  ONLY : pert
  USE westcom,              ONLY : npwq0
  USE pwcom,                ONLY : npwx
  USE gvect,                ONLY : gstart
  USE control_flags,        ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  COMPLEX(DP),INTENT(INOUT) :: ag(npwx,pert%nlocx)
  COMPLEX(DP),INTENT(IN) :: bg(npwx,pert%nlocx)
  INTEGER,INTENT(IN) :: l2_s,l2_e
  REAL(DP),INTENT(INOUT) :: c_distr(pert%nglob,pert%nlocx)
  INTEGER,INTENT(IN) :: g_s, g_e
  !
  ! Workspace
  !
  REAL(DP),EXTERNAL :: DDOT
  COMPLEX(DP),EXTERNAL :: ZDOTC
  INTEGER,ALLOCATABLE :: tmp_l2g(:)
  INTEGER :: il1, il2, ig1, ig2, icycl
  !
  CALL mp_barrier(world_comm)
  !
  CALL start_clock ('build_hr')
  !
  ! Initialize to zero
  !
  c_distr(:,l2_s:l2_e)=0.0_DP
  !
  ALLOCATE( tmp_l2g(1:pert%nlocx) )
  !
  tmp_l2g = 0
  DO il1 = 1, pert%nloc 
     tmp_l2g(il1) = pert%l2g(il1)
  ENDDO
  !
  DO icycl=0,nimage-1
     !
     DO il1=1,pert%nlocx
        !
        ig1 = tmp_l2g(il1)
        IF( ig1 == 0 .OR. ig1 > g_e ) CYCLE
        !
        DO il2=l2_s,l2_e
           !
           !ig2 = pert%l2g(il2)
           !IF( ig2 < n1 .OR. ig2 > n2 ) CYCLE
           !
           IF( gamma_only ) THEN 
              IF(gstart==1) THEN
                 c_distr(ig1,il2) = 2.0_DP * DDOT(2*npwq0,ag(1,il1),1,bg(1,il2),1)
              ELSE
                 c_distr(ig1,il2) = 2.0_DP * DDOT(2*npwq0-2,ag(2,il1),1,bg(2,il2),1)
              ENDIF
           ELSE
              c_distr(ig1,il2) = REAL( ZDOTC(npwq0,ag(1,il1),1,bg(1,il2),1) )
           ENDIF
           !
        ENDDO
     ENDDO
     !
     ! Cycle the ag array 
     ! 
     CALL mp_circular_shift_left( ag,      icycl,        inter_image_comm)
     CALL mp_circular_shift_left( tmp_l2g, icycl+nimage, inter_image_comm)
     !
  ENDDO
  !
  ! Syncronize c_distr
  !
  CALL mp_sum(c_distr(:,l2_s:l2_e), intra_bgrp_comm)
  !
  DEALLOCATE( tmp_l2g )
  !
  CALL stop_clock( "build_hr" )
  !
END SUBROUTINE
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
SUBROUTINE do_mgs(amat,m_global_start,m_global_end) 
  !
  ! MGS of the vectors beloging to the interval [ m_global_start, m_global_end ]
  !    also with respect to the vectors belonging to the interval [ 1, m_global_start -1 ]
  !
  USE kinds,                  ONLY : DP
  USE io_global,              ONLY : stdout
  USE mp_global,              ONLY : intra_bgrp_comm,inter_image_comm,my_image_id,nimage,world_comm 
  USE gvect,                  ONLY : gstart
  USE mp,                     ONLY : mp_sum,mp_barrier,mp_bcast
  USE westcom,                ONLY : npwq0
  USE pwcom,                  ONLY : npwx
  USE control_flags,          ONLY : gamma_only
  USE io_push,                ONLY : io_push_title 
  USE distribution_center,    ONLY : pert,distr_g2l
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: m_global_start,m_global_end
  COMPLEX(DP) :: amat(npwx,pert%nlocx)
  !
  ! Workspace 
  !
  INTEGER :: ig, ip, j, ncol
  INTEGER :: k_global, k_local, j_local, k_id
  REAL(DP) :: anorm, braket(pert%nloc)
  REAL(DP) :: time_spent(2)
  COMPLEX(DP) :: zbraket(pert%nloc)
  COMPLEX(DP),ALLOCATABLE :: vec(:)
  LOGICAL :: unfinished
  COMPLEX(DP) :: za
  !INTEGER :: m_local_end
  INTEGER :: m_local_start,m_local_end
  !
  REAL(DP),EXTERNAL :: DDOT
  COMPLEX(DP),EXTERNAL :: ZDOTC
  !
  CALL mp_barrier(world_comm)
  !
  CALL start_clock ('paramgs')
  !
  ALLOCATE( vec(npwx) )
  !
  ! 1) Run some checks
  !
  IF( m_global_start < 1 .OR. m_global_start > m_global_end .OR. m_global_end > pert%nglob ) &
     & CALL errore( 'mgs', 'do_mgs problem', 1 )
  !
  ! 2) Localize m_global_start
  !
  m_local_start = 1 
  DO ip = 1, pert%nloc
     ig = pert%l2g(ip)
     IF( ig < m_global_start ) CYCLE
     m_local_start = ip 
     EXIT
  ENDDO
  !
  ! 3) Localize m_global_end
  !
  m_local_end = pert%nloc
  DO ip = pert%nloc, 1, -1
     ig = pert%l2g(ip)
     IF( ig > m_global_end ) CYCLE
     m_local_end = ip
     EXIT
  ENDDO
  !
  j_local=1
  unfinished=.true.
  !
  !DO k_global=m_global_start,m_global_end
  DO k_global=1,m_global_end
     !
     CALL distr_g2l(k_global,k_local,k_id,nimage)
     !
     IF(my_image_id==k_id) THEN
        !
        ! 4) Eventually, normalize the current vector
        !
        IF( k_global >= m_global_start ) THEN
           !
           ! anorm = < k_l | k_l >
           !
           IF(gamma_only) THEN
              anorm = 2.0_DP * DDOT(2*npwq0,amat(1,k_local),1,amat(1,k_local),1)
              IF(gstart==2) anorm = anorm - REAL(amat(1,k_local),KIND=DP) * REAL(amat(1,k_local),KIND=DP)
           ELSE
              anorm = DDOT(2*npwq0,amat(1,k_local),1,amat(1,k_local),1)
           ENDIF
           !
           CALL mp_sum(anorm,intra_bgrp_comm)
           !
           ! normalize | k_l >
           !
           za = CMPLX( 1._DP/SQRT(anorm), 0._DP,KIND=DP)
           CALL ZSCAL(npwq0,za,amat(1,k_local),1)
           !
        ENDIF
        !
        ! 5) Copy the current vector into V
        !
        CALL ZCOPY(npwx,amat(1,k_local),1,vec(1),1)
        !
        j_local=MAX(k_local+1,m_local_start)
        !
        !IF(j_local>pert%nloc) unfinished=.false.
        IF(j_local>m_local_end) unfinished=.false.
        !
     ENDIF
     !
     ! BCAST | vec >
     !
     CALL mp_bcast(vec,k_id,inter_image_comm)
     !
     ! Update when needed
     !
     IF(unfinished) THEN
        !
        ! IN the range ip=j_local:pert%nloc    = >    | ip > = | ip > - | vec > * < vec | ip >
        !
        IF(gamma_only) THEN
           DO ip = j_local,m_local_end !pert%nloc
              braket(ip) = 2._DP * DDOT(2*npwq0,vec(1),1,amat(1,ip),1)
           ENDDO
           !IF (gstart==2) FORALL( ip=j_local:pert%nloc ) braket(ip) = braket(ip) - REAL(vec(1),KIND=DP)*REAL(amat(1,ip),KIND=DP) 
           IF (gstart==2) FORALL( ip=j_local:m_local_end ) braket(ip) = braket(ip) - REAL(vec(1),KIND=DP)*REAL(amat(1,ip),KIND=DP) 
           !CALL mp_sum(braket(j_local:pert%nloc),intra_bgrp_comm)
           CALL mp_sum(braket(j_local:m_local_end),intra_bgrp_comm)
           !FORALL(ip=j_local:pert%nloc) zbraket(ip) = CMPLX( braket(ip), 0._DP, KIND=DP)
           FORALL(ip=j_local:m_local_end) zbraket(ip) = CMPLX( braket(ip), 0._DP, KIND=DP)
        ELSE
           DO ip = j_local,m_local_end !pert%nloc
              zbraket(ip) = ZDOTC(npwq0,vec(1),1,amat(1,ip),1)
           ENDDO
           !CALL mp_sum(zbraket(j_local:pert%nloc),intra_bgrp_comm)
           CALL mp_sum(zbraket(j_local:m_local_end),intra_bgrp_comm)
        ENDIF
        !
        ncol=m_local_end-j_local+1
        !ncol=m_global_end-j_local+1
        CALL ZGERU(npwx,ncol,(-1._DP,0._DP),vec(1),1,zbraket(j_local),1,amat(1,j_local),npwx)
        !
     ENDIF
     !
  ENDDO
  !
  DEALLOCATE( vec )
  !
  CALL mp_barrier(world_comm)
  !
  CALL stop_clock ('paramgs')
  !
END SUBROUTINE
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
SUBROUTINE do_randomize ( amat, mglobalstart, mglobalend  ) 
  !
  ! Randomize in dvg the vectors belonging to [ mglobalstart, mglobalend ]
  !
  USE kinds,                ONLY : DP
  USE random_numbers,       ONLY : randy
  USE gvect,                ONLY : g,gstart,ngm_g
  USE westcom,              ONLY : npwq0,q0ig_l2g,dvg
  USE constants,            ONLY : tpi
  USE mp,                   ONLY : mp_barrier
  USE mp_global,            ONLY : world_comm
  USE distribution_center,  ONLY : pert
  USE pwcom,                ONLY : npwx
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER,INTENT(IN) :: mglobalstart, mglobalend
  COMPLEX(DP) :: amat(npwx,pert%nlocx)
  !
  ! Workspace
  !
  REAL(DP),ALLOCATABLE :: random_num_debug(:,:)
  INTEGER :: il1,ig1,ig
  REAL(DP) :: aux_real
  REAL(DP) :: rr, arg
  !
  CALL mp_barrier(world_comm)
  !
  CALL start_clock ('randomize')
  !
  ! Random numbers are generated according to G-global
  !
  ALLOCATE(random_num_debug(2,ngm_g))
  !
  DO il1=1,pert%nloc
     !
     ig1=pert%l2g(il1)
     IF( ig1 < mglobalstart .OR. ig1 > mglobalend ) CYCLE
     !
     ! Initialize the sequence
     !
     aux_real=randy(ig1)
     !
     DO ig=1,ngm_g
        random_num_debug(1:2,ig) = (/ randy(), randy() /)
     ENDDO
     !
     amat(:,il1) = 0._DP
#ifdef __OPENMP
!$OMP PARALLEL private(ig,rr,arg)
!$OMP DO
#endif
     DO ig=gstart,npwq0
        rr = random_num_debug(1,q0ig_l2g(ig))
        arg = tpi * random_num_debug(2,q0ig_l2g(ig))
        amat(ig,il1) = CMPLX( rr*COS( arg ), rr*SIN( arg ), KIND=DP) / &
                      ( g(1,ig)*g(1,ig) + &
                        g(2,ig)*g(2,ig) + &
                        g(3,ig)*g(3,ig) + 1.0_DP )
     ENDDO
#ifdef __OPENMP
!$OMP ENDDO
!$OMP END PARALLEL
#endif
     !
  ENDDO
  !
  DEALLOCATE(random_num_debug)
  !
  CALL stop_clock ('randomize')
  !
END SUBROUTINE
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
SUBROUTINE output_a_report(iteration)
   !
   USE kinds,        ONLY : DP
   USE westcom,      ONLY : n_pdep_eigen,ev,conv
   USE west_io ,     ONLY : serial_table_output
   USE mp_world,     ONLY : mpime,root
   !
   IMPLICIT NONE
   !
   INTEGER,INTENT(IN) :: iteration
   CHARACTER(LEN=9) :: pref
   INTEGER :: ip
   REAL(DP) :: out_tab(n_pdep_eigen,3)
   !
   DO ip=1,n_pdep_eigen
      out_tab(ip,1) = REAL(ip,DP)
      out_tab(ip,2) = ev(ip)
      out_tab(ip,3) = 0._DP
      IF(conv(ip)) out_tab(ip,3) = 1._DP
   ENDDO
   IF(iteration>=0) THEN 
      WRITE(pref,"('itr_',i5.5)") iteration
   ELSE
      pref='converged'
   ENDIF
   CALL serial_table_output(mpime==root,4000,'wstat.'//TRIM(ADJUSTL(pref)),out_tab,n_pdep_eigen,3,(/'   iprt','eigenv.','  conv.'/))
   !
END SUBROUTINE
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
SUBROUTINE output_ev_and_time(nvec,ev,time)
   !
   USE kinds,                ONLY : DP
   USE io_global,            ONLY : stdout
   USE io_push,              ONLY : io_push_title,io_push_bar
   !
   IMPLICIT NONE
   !
   INTEGER,INTENT(IN) :: nvec
   REAL(DP),INTENT(IN) :: ev(nvec)
   REAL(DP),INTENT(IN) :: time(2)
   !
   INTEGER :: i,j
   CHARACTER(20),EXTERNAL :: human_readable_time
   !
   WRITE(stdout,'(7X,a)') ' '
   DO i = 1, INT( nvec / 9 )
      WRITE(stdout,'(6X, 9(f9.5,1x))') (ev(j), j=9*(i-1)+1,9*i)
   ENDDO
   IF( MOD(nvec,9) > 0 ) WRITE(stdout,'(6X, 9(f9.5,1x))') (ev(j), j=9*INT(nvec/9)+1,nvec)
   WRITE(stdout,'(7X,a)') ' '
   CALL io_push_bar()
   WRITE(stdout, "(5x,'Tot. elapsed time ',a,',  time spent in last iteration ',a) ") &
   TRIM(human_readable_time(time(2))), TRIM(human_readable_time(time(2)-time(1)))
   CALL io_push_bar()
   !
END SUBROUTINE
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
SUBROUTINE diagox( nbase, nvec, hr_distr, nvecx, ew, vr_distr  ) 
  !
  USE kinds,                 ONLY : DP
  USE io_global,             ONLY : stdout
  USE io_push,               ONLY : io_push_title,io_push_bar
  USE mp_global,             ONLY : nproc_bgrp,inter_image_comm,nimage 
  USE mp_world,              ONLY : nproc
  USE westcom,               ONLY : n_pdep_eigen
  USE linear_algebra_kernel, ONLY : matdiago_dsy
  USE io_push,               ONLY : io_push_title,io_push_bar
  USE distribution_center,   ONLY : pert
  USE mp,                    ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER  :: nbase, nvec, nvecx
  REAL(DP) :: hr_distr(nvecx,pert%nlocx)
  REAL(DP) :: ew(nvecx),vr_distr(nvecx,pert%nlocx)
  !
  REAL(DP) :: time_spent(2)
  REAL(DP), EXTERNAL :: GET_CLOCK
  CHARACTER(20),EXTERNAL :: human_readable_time
  LOGICAL :: l_parallel 
  INTEGER :: npur, npuc
  CHARACTER(LEN=8) :: aux_label_npur
  CHARACTER(LEN=8) :: aux_label_npuc
  !
  !
#if defined __SCALAPACK
  l_parallel = nproc > 3 
#else
  l_parallel = .FALSE.
#endif
  !
  CALL start_clock( 'diagox' )
  !
  IF( l_parallel ) THEN 
     !
     time_spent(1)=get_clock( 'diagox' )
#if defined __SCALAPACK
     !
     npur = MAX( MIN(INT( SQRT(DBLE(nproc))),INT(DBLE(nbase)/32.)), 2)
     npuc = npur
     CALL parallel_distributed_diago ( nvec, nbase, nvecx, pert%nlocx, hr_distr, vr_distr, ew(1:nvec), npur, npuc, pert)
     !
#endif
     time_spent(2)=get_clock( 'diagox' )
     WRITE(aux_label_npur,'(i8)') npur
     WRITE(aux_label_npuc,'(i8)') npuc
     WRITE(stdout, "(5x, 'p-DIAGOX done in ',a,' with a SCALAPACK grid ',a)") &
      & TRIM(ADJUSTL(human_readable_time(time_spent(2)-time_spent(1)))), &
      & "("//TRIM(ADJUSTL(aux_label_npur))//"x"//TRIM(ADJUSTL(aux_label_npuc))//")"
     !
  ELSE
     !
     time_spent(1)=get_clock( 'diagox' )
     !
     CALL serial_diagox( nvec, nbase, nvecx, hr_distr, ew(1:nvec), vr_distr )
     !
     time_spent(2)=get_clock( 'diagox' )
     WRITE(stdout, "(5x, 's-DIAGOX done in ',a20)") human_readable_time(time_spent(2)-time_spent(1))
     !
  ENDIF
  !
  CALL stop_clock( 'diagox' )
  !
  CALL io_push_bar()
  !
END SUBROUTINE
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
SUBROUTINE serial_diagox ( nselect, n, lda, hr_distr, e, vr_distr ) 
  !
  ! Diagox -- serial
  !   nselect : number of wanted ev
  !   n       : actual dimension of a
  !   lda     : leading dimension of a
  !   a       : matrix to be diagox
  !   e       : eigenval(1:nselsect), even though it is defined 1:lda
  !   z       : unitary trans. 
  !
  USE kinds,                 ONLY : DP
  USE mp,                    ONLY : mp_bcast,mp_sum
  USE mp_global,             ONLY : me_bgrp, root_bgrp, intra_bgrp_comm, inter_image_comm
  USE linear_algebra_kernel, ONLY : matdiago_dsy 
  USE distribution_center,   ONLY : pert
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER,INTENT(IN) :: nselect,lda,n
  REAL(DP),INTENT(IN) :: hr_distr(lda,pert%nlocx)
  REAL(DP),INTENT(OUT) :: e(nselect)
  REAL(DP),INTENT(OUT) :: vr_distr(lda,pert%nlocx)
  !
  ! Workspace
  !
  REAL(DP),ALLOCATABLE :: zz(:,:),ee(:)
  INTEGER :: i, j, il1, il2, ig1, ig2
  !
  ALLOCATE(zz(n,n))
  ALLOCATE(ee(n))
  !
  ! Create zz from hr_distr
  !
  zz = 0._DP
  DO il2 = 1, pert%nloc
     ig2 = pert%l2g(il2)
     IF(ig2>n) CYCLE
     DO ig1 = 1, pert%nglob
        IF(ig1>n) CYCLE
        zz( ig1,ig2 ) = hr_distr( ig1, il2)
     ENDDO
  ENDDO
  CALL mp_sum(zz,inter_image_comm)
  ee = 0._DP
  !
  IF(me_bgrp == root_bgrp ) THEN 
     !
     CALL matdiago_dsy(n,zz,ee,.FALSE.)
     !
  ENDIF
  !
  CALL mp_bcast( ee, root_bgrp, intra_bgrp_comm )
  CALL mp_bcast( zz, root_bgrp, intra_bgrp_comm )
  !
  DO il2 = 1, pert%nloc
     ig2 = pert%l2g(il2)
     IF(ig2>nselect) CYCLE
     DO ig1 = 1, pert%nglob
        IF(ig1>n) CYCLE
        vr_distr(ig1,il2) = zz(ig1,ig2)
     ENDDO
  ENDDO
  DO i = 1, nselect
     e(i) = ee(i)
  ENDDO
  !
  DEALLOCATE(zz)
  DEALLOCATE(ee)
  !
END SUBROUTINE
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
SUBROUTINE redistribute_vr_distr( nselect, n, lda, vr_distr, ishift)
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : my_image_id,inter_image_comm,intra_bgrp_comm,nimage,world_comm
  USE mp,                   ONLY : mp_sum,mp_circular_shift_left,mp_barrier
  USE distribution_center,  ONLY : pert
  USE westcom,              ONLY : npwq0
  USE pwcom,                ONLY : npwx
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER,INTENT(IN) :: nselect, n, lda
  REAL(DP),INTENT(INOUT) :: vr_distr(lda,pert%nlocx)
  INTEGER,INTENT(IN) :: ishift(lda)
  !
  ! Workspace
  !
  REAL(DP),ALLOCATABLE :: tmp_distr(:,:)
  INTEGER,ALLOCATABLE :: tmp_l2g(:)
  INTEGER :: il1, il2, ig1, ig2, icycl
  !
  CALL mp_barrier( world_comm )
  !  
  CALL start_clock( 'redistr_vr' )
  !
  ALLOCATE( tmp_distr(lda,pert%nlocx) )
  tmp_distr = vr_distr
  ALLOCATE( tmp_l2g(1:pert%nlocx) )
  !
  tmp_l2g = 0
  !
  DO il1 = 1, pert%nloc 
     tmp_l2g(il1) = pert%l2g(il1)
  ENDDO
  !
  DO icycl=0,nimage-1
     !
     DO il2=1,pert%nloc
        ig2 = pert%l2g(il2)
        IF( ig2 <= n .OR. ig2 > n+nselect ) CYCLE
        !
        DO il1=1,pert%nlocx
           ig1 = tmp_l2g(il1)
           !IF( ig1 .NE. ig2-n ) CYCLE
           IF( ig1 == 0 ) CYCLE
           IF( ig1 .NE. ishift(ig2) ) CYCLE
           !
           vr_distr(:,il2) = tmp_distr(:,il1)
           !
        ENDDO
     ENDDO
     !
     ! Cycle the tmp_distr array 
     ! 
     CALL mp_circular_shift_left( tmp_distr ,        icycl, inter_image_comm)
     CALL mp_circular_shift_left( tmp_l2g   , icycl+nimage, inter_image_comm)
     !
  ENDDO
  !
  DEALLOCATE( tmp_distr )
  DEALLOCATE( tmp_l2g )
  !
  CALL stop_clock( "redistr_vr" )
  !
END SUBROUTINE
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
SUBROUTINE symm_hr_distr( hr_distr, n, lda )
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : my_image_id,inter_image_comm,intra_bgrp_comm,nimage,world_comm
  USE mp,                   ONLY : mp_sum,mp_circular_shift_left,mp_barrier
  USE distribution_center,  ONLY : pert
  USE westcom,              ONLY : npwq0
  USE pwcom,                ONLY : npwx
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER,INTENT(IN) :: n, lda
  REAL(DP),INTENT(INOUT) :: hr_distr(lda,pert%nlocx)
  !
  ! Workspace
  !
  REAL(DP),ALLOCATABLE :: tmp_distr(:,:)
  INTEGER,ALLOCATABLE :: tmp_l2g(:)
  INTEGER :: il1, il2, ig1, ig2, icycl
  !
  CALL mp_barrier( world_comm )
  !  
  CALL start_clock( 'symm_hr' )
  !
  ALLOCATE( tmp_distr(lda,pert%nlocx) )
  tmp_distr = hr_distr
  ALLOCATE( tmp_l2g(1:pert%nlocx) )
  !
  tmp_l2g = 0
  !
  DO il1 = 1, pert%nloc 
     tmp_l2g(il1) = pert%l2g(il1)
  ENDDO
  !
  DO icycl=0,nimage-1
     !
     DO il1=1,pert%nloc
        ig1 = pert%l2g(il1)
        IF( ig1 > n ) CYCLE
        !
        DO il2=1,pert%nlocx
           ig2 = tmp_l2g(il2)
           IF( .NOT.ig2 > ig1 ) CYCLE
           !
           hr_distr(ig2,il1) = tmp_distr(ig1,il2)
           !
        ENDDO
     ENDDO
     !
     ! Cycle the tmp_distr array 
     ! 
     CALL mp_circular_shift_left( tmp_distr ,        icycl, inter_image_comm)
     CALL mp_circular_shift_left( tmp_l2g   , icycl+nimage, inter_image_comm)
     !
  ENDDO
  !
  DEALLOCATE( tmp_distr )
  DEALLOCATE( tmp_l2g )
  !
  CALL stop_clock( "symm_hr" )
  !
END SUBROUTINE
