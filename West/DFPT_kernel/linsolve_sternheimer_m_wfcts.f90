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
SUBROUTINE linsolve_sternheimer_m_wfcts (nbndval, m, b, x, h_diag, e, conv_root, anorm)
  !----------------------------------------------------------------------
  !
  !     iterative solution of the linear system:
  !
  !                 A * | x_i > = | b_i >       i=(1:m)
  !
  !     with :      A = ( H - E_i + alpha * P_v ) 
  !
  !                 where H is a complex hermitean matrix (H_{SCF}), E_v is a real scalar (energy of band v)
  !                 x and b are complex vectors
  !
  !     on input:
  !
  !                 e        real     unperturbed eigenvalue.
  !
  !                 x        contains an estimate of the solution
  !                          vector.
  !
  !                 b        contains the right hand side vector
  !                          of the system.
  !
  !                 ndmx     integer row dimension of x, ecc.
  !
  !                 ndm      integer actual row dimension of x
  !
  !     on output:  x        contains the refined estimate of the
  !                          solution vector.
  !
  !
  USE kinds,                ONLY : DP
  USE mp_global,            ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE gvect,                ONLY : gstart
  USE wavefunctions_module, ONLY : evc
  USE control_flags,        ONLY : gamma_only 
  USE westcom,              ONLY : tr2_dfpt,n_dfpt_maxiter,alphapv_dfpt
  USE pwcom,                ONLY : npw,npwx
  USE noncollin_module,     ONLY : npol, noncolin
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER,INTENT(IN)  :: nbndval, & ! input: the number of occupied bands
                         m ! the number of bands to be computed
  REAL(DP),INTENT(IN)  :: e(m), & ! input: the actual eigenvalue
                          h_diag(npwx*npol,m) ! input: an estimate of ( H - \epsilon )
  COMPLEX(DP),INTENT(IN)  :: b (npwx*npol, m)   ! input: the known term
  !
  REAL(DP),INTENT(OUT) :: anorm      ! output: the norm of the error in the solution
  COMPLEX(DP),INTENT(INOUT) :: x (npwx*npol, m) ! output: the solution of the linear syst
  LOGICAL,INTENT(OUT) :: conv_root ! output: if true the root is converged
  !
  ! Workspace
  !
  INTEGER :: iter, ibnd, lbnd
  ! counters on iteration, bands
  INTEGER :: ig
  INTEGER :: conv(m)
  ! if 1 the root is converged
  COMPLEX(DP),ALLOCATABLE :: g(:,:), t(:,:), h(:,:), hold(:,:)
  !  the gradient of psi
  !  the preconditioned gradient
  !  the delta gradient
  !  the conjugate gradient
  !  work space
  COMPLEX(DP) ::  dcgamma, dclambda
  !  the ratio between rho
  !  step length
  REAL(KIND=DP),EXTERNAL :: DDOT
  COMPLEX(KIND=DP),EXTERNAL :: ZDOTC
  !  the scalar product
  REAL(DP),ALLOCATABLE :: rho(:), rhoold(:), eu(:), a(:), c(:)
  ! the residue
  ! auxiliary for h_diag
  REAL(DP) :: energy, norma
  !
  CALL start_clock ('linstern')
  !
  ! Zeros
  !
  conv=0
  conv_root=.FALSE.
  !
  ALLOCATE( g(npwx*npol,m), t(npwx*npol,m), h(npwx*npol,m), hold(npwx*npol,m) )
  ALLOCATE( rho(m), rhoold(m), eu(m), a(m), c(m) )
  !
  g = 0.0_DP
  t = 0.0_DP
  h = 0.0_DP
  hold = 0.0_DP
  !
  ! Step 1, initialization of the loop
  !
  CALL apply_sternheimerop_to_m_wfcs(nbndval, x, g, e, alphapv_dfpt, m)
  !CALL cgoperator_at_gamma (nbndval, ndm, x, g, e, alpha, nbndval)
  !  FORALL(ibnd=1:nbndval,ig=1:ndm) g(ig,ibnd) = g(ig,ibnd) - b(ig,ibnd)
  DO ibnd=1,m
     CALL ZAXPY(npw,(-1._DP,0._DP),b(1,ibnd),1,g(1,ibnd),1)
  ENDDO
  IF(npol==2) THEN
     DO ibnd=1,m
        CALL ZAXPY(npw,(-1._DP,0._DP),b(npwx+1,ibnd),1,g(npwx+1,ibnd),1)
     ENDDO
  ENDIF
  !
  ! Loop
  !
  DO iter = 1, n_dfpt_maxiter
     !
     !    compute preconditioned residual vector and convergence check
     !
     lbnd = 0
     DO ibnd = 1, m
        IF (conv (ibnd) .EQ. 0) THEN
           lbnd = lbnd+1
           !
#ifdef __OPENMP
!$OMP PARALLEL DO
#endif
           DO ig=1,npwx*npol 
              h(ig,ibnd) = g(ig,ibnd) * h_diag(ig,ibnd)
           ENDDO
#ifdef __OPENMP
!$OMP END PARALLEL DO
#endif
           !
           IF (gamma_only) THEN
              rho(lbnd)=2.0_DP * DDOT(2*npw,h(1,ibnd),1,g(1,ibnd),1)
              IF(gstart==2) THEN
                 rho(lbnd)=rho(lbnd)-REAL(h(1,ibnd),KIND=DP)*REAL(g(1,ibnd),KIND=DP)
              ENDIF
           ELSE
              rho(lbnd) = ZDOTC (npw, h(1,ibnd), 1, g(1,ibnd), 1)
              IF(noncolin) rho(lbnd) = rho(lbnd) + ZDOTC (npw, h(1+npwx,ibnd), 1, g(1+npwx,ibnd), 1)
           ENDIF
           !
        ENDIF
     ENDDO
     !
     CALL mp_sum(  rho(1:lbnd) , intra_bgrp_comm )
     !
     DO ibnd = m, 1, -1
        IF (conv(ibnd) .EQ. 0) THEN
           rho(ibnd)=rho(lbnd)
           lbnd = lbnd -1
           anorm = SQRT (rho (ibnd) )
           IF (anorm .LT. tr2_dfpt) conv (ibnd) = 1
        ENDIF
     ENDDO
     !
     conv_root = .TRUE.
     DO ibnd = 1, m
        conv_root = conv_root .AND. (conv (ibnd) .EQ. 1)
     ENDDO
     IF (conv_root) GOTO 100
     !
     !        compute the step direction h. Conjugate it to previous step
     !
     lbnd = 0
     DO ibnd = 1, m
        IF (conv (ibnd) .EQ. 0) THEN
           !
           !          change sign to h
           !
           CALL DSCAL (2*npwx*npol, - 1._DP, h (1, ibnd), 1)
           IF (iter .NE. 1) THEN
              dcgamma = CMPLX( rho (ibnd) / rhoold (ibnd), 0.0_DP, KIND=DP) 
              CALL ZAXPY (npwx*npol, dcgamma, hold (1, ibnd), 1, h (1, ibnd), 1)
           ENDIF
           !
           ! here hold is used as auxiliary vector in order to efficiently compute t = A*h
           ! it is later set to the current (becoming old) value of h
           !
           lbnd = lbnd+1
           CALL ZCOPY(npwx*npol, h (1, ibnd), 1, hold (1, lbnd), 1)
           eu (lbnd) = e (ibnd)
        ENDIF
     ENDDO
     !
     !        compute t = A*h
     !
     CALL apply_sternheimerop_to_m_wfcs(nbndval, hold, t, eu, alphapv_dfpt, lbnd)
     !
     !        compute the coefficients a and c for the line minimization
     !        compute step length lambda
     lbnd=0
     DO ibnd = 1, m
        IF (conv (ibnd) .EQ. 0) THEN
           lbnd=lbnd+1
           IF (gamma_only) THEN
              a(lbnd) = 2.0_DP * DDOT(2*npw,h(1,ibnd),1,g(1,ibnd),1)
              c(lbnd) = 2.0_DP * DDOT(2*npw,h(1,ibnd),1,t(1,lbnd),1)
              IF (gstart == 2) THEN
                 a(lbnd)=a(lbnd)-REAL(h(1,ibnd),KIND=DP)*REAL(g(1,ibnd),KIND=DP)
                 c(lbnd)=c(lbnd)-REAL(h(1,ibnd),KIND=DP)*REAL(t(1,lbnd),KIND=DP)
              ENDIF
           ELSE
              a(lbnd) = ZDOTC (npw, h(1,ibnd), 1, g(1,ibnd), 1)
              c(lbnd) = ZDOTC (npw, h(1,ibnd), 1, t(1,lbnd), 1)
              IF(noncolin) THEN
                 a(lbnd) = a(lbnd) + ZDOTC (npw, h(1+npwx,ibnd), 1, g(1+npwx,ibnd), 1)
                 c(lbnd) = c(lbnd) + ZDOTC (npw, h(1+npwx,ibnd), 1, t(1+npwx,lbnd), 1)
              ENDIF
           ENDIF
        ENDIF
     ENDDO
     !
     CALL mp_sum(  a(1:lbnd), intra_bgrp_comm )
     CALL mp_sum(  c(1:lbnd), intra_bgrp_comm )
     !
     lbnd=0
     DO ibnd = 1, m
        IF (conv (ibnd) .EQ. 0) THEN
           lbnd=lbnd+1
           dclambda = CMPLX( - a(lbnd) / c(lbnd), 0.0_DP, KIND=DP)
           !
           !    move to new position
           !
           CALL ZAXPY(npwx*npol, dclambda, h(1,ibnd), 1, x(1,ibnd), 1)
           !
           !    update to get the gradient
           !
           CALL ZAXPY(npwx*npol, dclambda, t(1,lbnd), 1, g(1,ibnd), 1)
           !
           !    save current (now old) h and rho for later use
           !
           CALL ZCOPY(npwx*npol, h(1,ibnd), 1, hold(1,ibnd), 1)
           rhoold (ibnd) = rho (ibnd)
           !
        ENDIF
        !
     ENDDO ! on bands  
     !
  ENDDO ! on iterations
  !
  !
100 CONTINUE
  !
  DEALLOCATE( g, t, h, hold )
  DEALLOCATE( rho, rhoold, eu, a, c )
  !
  CALL stop_clock ('linstern')
  !
END SUBROUTINE
