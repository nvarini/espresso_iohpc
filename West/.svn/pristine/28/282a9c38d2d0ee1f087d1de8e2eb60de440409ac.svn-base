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
MODULE fft_at_gamma
  !-----------------------------------------------------------------------
  !
  ! Everything is done following dffts
  !
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : nl,nlm
  USE gvecs,                ONLY : nls,nlsm
  USE fft_interfaces,       ONLY : fwfft,invfft
  USE fft_base,             ONLY : dffts,dfftp
  !
  IMPLICIT NONE
  !
  !
  CONTAINS
  !
  !
  SUBROUTINE doubleband_invfft(n,a1,a2,lda,b,ldb)
    !
    ! INVFFT : G ---> R
    !
    ! INPUT  : n     = actual number of PW
    !          a1    = COMPLEX array containing ONE COMPLEX function in G space
    !          a2    = COMPLEX array containing ONE COMPLEX function in G space
    !          lda   = leading dimension of a1 or a2
    !          ldb   = leading dimension of b
    ! OUTPUT : b     = ONE COMPLEX array containing TWO REAL functions in R space
    !
    IMPLICIT NONE
    !
    INTEGER,INTENT(IN) :: lda, n, ldb
    COMPLEX(DP),INTENT(IN) :: a1(lda)
    COMPLEX(DP),INTENT(IN) :: a2(lda)
    COMPLEX(DP),INTENT(OUT) :: b(ldb)
    !
    INTEGER :: ig
    !
#ifdef __OPENMP
!$OMP PARALLEL private(ig)
!$OMP DO
#endif
    DO ig=1,ldb 
       b(ig)= (0.0_DP,0.0_DP) 
    ENDDO
#ifdef __OPENMP
!$OMP ENDDO
!$OMP DO
#endif
    DO ig=1,n 
       b(nls (ig))=        a1(ig) + (0.0_DP,1.0_DP) * a2(ig)
       b(nlsm(ig))= DCONJG(a1(ig) - (0.0_DP,1.0_DP) * a2(ig))
    ENDDO
#ifdef __OPENMP
!$OMP ENDDO
!$OMP END PARALLEL
#endif
    !
    CALL invfft('Wave', b, dffts)
    !
  END SUBROUTINE
  !
  !
  SUBROUTINE doubleband_fwfft(n,a,lda,b1,b2,ldb)
    !
    ! FWFFT : R ---> G
    !
    ! INPUT  : n     = actual number of PW
    !          a     = ONE COMPLEX array containing TWO REAL functions in R space
    !          lda   = leading dimension of a
    !          ldb   = leading dimension of b1 or b2 
    ! OUTPUT : b1     = ONE COMPLEX array containing ONE COMPLEX function in G space
    !          b2     = ONE COMPLEX array containing ONE COMPLEX function in G space
    !
    IMPLICIT NONE
    !
    INTEGER,INTENT(IN) :: lda, n, ldb
    COMPLEX(DP),INTENT(INOUT) :: a(lda)
    COMPLEX(DP),INTENT(OUT) :: b1(ldb)
    COMPLEX(DP),INTENT(OUT) :: b2(ldb)
    !
    INTEGER :: ig
    COMPLEX(DP) :: fm, fp 
    !
    ! FFT call
    CALL fwfft('Wave',a, dffts)
    ! Keep only G>=0
    !
#ifdef __OPENMP
!$OMP PARALLEL private(ig,fp,fm)
!$OMP DO
#endif
    DO ig = 1, n
       fp = ( a(nls (ig)) + a(nlsm(ig)) )*0.5_DP
       fm = ( a(nls (ig)) - a(nlsm(ig)) )*0.5_DP
       b1(ig) = CMPLX( DBLE(fp), DIMAG(fm), KIND=DP)
       b2(ig) = CMPLX(DIMAG(fp), -DBLE(fm), KIND=DP)
    ENDDO
#ifdef __OPENMP
!$OMP ENDDO
!$OMP END PARALLEL
#endif
    !
    DO ig = (n+1), ldb
       b1(ig) = (0.0_DP,0.0_DP)
       b2(ig) = (0.0_DP,0.0_DP)
    ENDDO
    !
  END SUBROUTINE
  !
  !
  SUBROUTINE singleband_invfft(n,a1,lda,b,ldb)
    !
    ! INVFFT : G ---> R
    !
    ! INPUT  : n     = actual number of PW
    !          a1     = ONE COMPLEX arrays containing ONE COMPLEX functions in G space
    !          lda   = leading dimension of a1
    !          ldb   = leading dimension of b          
    ! OUTPUT : b     = ONE COMPLEX array containing ONE REAL functions in R space + 0
    !
    IMPLICIT NONE
    !
    INTEGER,INTENT(IN) :: lda, n, ldb
    COMPLEX(DP),INTENT(IN) :: a1(lda)
    COMPLEX(DP),INTENT(OUT) :: b(ldb)
    !
    INTEGER :: ig
    !
    !
#ifdef __OPENMP
!$OMP PARALLEL private(ig)
!$OMP DO
#endif
    DO ig=1,ldb 
       b(ig)= (0.0_DP,0.0_DP) 
    ENDDO
#ifdef __OPENMP
!$OMP ENDDO
!$OMP DO
#endif
    DO ig=1,n
       b(nls (ig))=        a1(ig)
       b(nlsm(ig))= DCONJG(a1(ig))
    ENDDO
#ifdef __OPENMP
!$OMP ENDDO
!$OMP END PARALLEL
#endif
    !
    CALL invfft('Wave', b, dffts)
    !
  END SUBROUTINE
  !
  !
  SUBROUTINE singleband_fwfft(n,a,lda,b1,ldb)
    !
    ! FWFFT : R ---> G
    !
    ! INPUT  : n     = actual number of PW
    !          a     = ONE COMPLEX array containing ONE REAL functions in R space + 0
    !          lda   = leading dimension of a
    !          ldb   = leading dimension of b1 
    ! OUTPUT : b1     = ONE COMPLEX array containing ONE COMPLEX functions in G space
    !
    IMPLICIT NONE
    !
    INTEGER,INTENT(IN) :: lda, n, ldb
    COMPLEX(DP),INTENT(INOUT) :: a(lda)
    COMPLEX(DP),INTENT(OUT) :: b1(ldb)
    !
    INTEGER :: ig
    !
    !
    ! FFT call
    CALL fwfft('Wave',a, dffts)
    ! Keep only G>=0
    !
#ifdef __OPENMP
!$OMP PARALLEL private(ig)
!$OMP DO
#endif
    DO ig=1,n 
       b1(ig) = a(nls(ig))
    ENDDO
#ifdef __OPENMP
!$OMP ENDDO
!$OMP END PARALLEL
#endif
    !
    DO ig = (n+1), ldb
       b1(ig) = (0.0_DP,0.0_DP)
    ENDDO
    !
  END SUBROUTINE
  !
  !
  SUBROUTINE dense_fwfft(n,a,lda,b1,ldb)
    !
    ! FWFFT : R ---> G
    !
    ! INPUT  : n     = actual number of PW
    !          a     = ONE COMPLEX array containing ONE REAL functions in R space + 0
    !          lda   = leading dimension of a
    !          ldb   = leading dimension of b1 
    ! OUTPUT : b1     = ONE COMPLEX array containing ONE COMPLEX functions in G space
    !
    IMPLICIT NONE
    !
    INTEGER,INTENT(IN) :: lda, n, ldb
    COMPLEX(DP),INTENT(INOUT) :: a(lda)
    COMPLEX(DP),INTENT(OUT) :: b1(ldb)
    !
    INTEGER :: ig
    !
    !
    ! FFT call
    CALL fwfft('Dense',a, dfftp)
    ! Keep only G>=0
    !
#ifdef __OPENMP
!$OMP PARALLEL private(ig)
!$OMP DO
#endif
    DO ig=1,n 
       b1(ig) = a(nl(ig))
    ENDDO
#ifdef __OPENMP
!$OMP ENDDO
!$OMP END PARALLEL
#endif
    !
    DO ig = (n+1), ldb
       b1(ig) = (0.0_DP,0.0_DP)
    ENDDO
    !
  END SUBROUTINE
  !
  !
  SUBROUTINE smooth_fwfft(n,a,lda,b1,ldb)
    !
    ! FWFFT : R ---> G
    !
    ! INPUT  : n     = actual number of PW
    !          a     = ONE COMPLEX array containing ONE REAL functions in R space + 0
    !          lda   = leading dimension of a
    !          ldb   = leading dimension of b1 
    ! OUTPUT : b1     = ONE COMPLEX array containing ONE COMPLEX functions in G space
    !
    IMPLICIT NONE
    !
    INTEGER,INTENT(IN) :: lda, n, ldb
    COMPLEX(DP),INTENT(INOUT) :: a(lda)
    COMPLEX(DP),INTENT(OUT) :: b1(ldb)
    !
    INTEGER :: ig
    !
    !
    ! FFT call
    CALL fwfft('Smooth',a, dffts)
    ! Keep only G>=0
    !
#ifdef __OPENMP
!$OMP PARALLEL private(ig)
!$OMP DO
#endif
    DO ig=1,n 
       b1(ig) = a(nls(ig))
    ENDDO
#ifdef __OPENMP
!$OMP ENDDO
!$OMP END PARALLEL
#endif
    !
    DO ig = (n+1), ldb
       b1(ig) = (0.0_DP,0.0_DP)
    ENDDO
    !
  END SUBROUTINE
  !
  !
  SUBROUTINE dense_invfft(n,a1,lda,b,ldb)
    !
    ! INVFFT : G ---> R
    !
    ! INPUT  : n     = actual number of PW
    !          a1     = ONE COMPLEX arrays containing ONE COMPLEX functions in G space
    !          lda   = leading dimension of a1
    !          ldb   = leading dimension of b          
    ! OUTPUT : b     = ONE COMPLEX array containing ONE REAL functions in R space + 0
    !
    IMPLICIT NONE
    !
    INTEGER,INTENT(IN) :: lda, n, ldb
    COMPLEX(DP),INTENT(IN) :: a1(lda)
    COMPLEX(DP),INTENT(OUT) :: b(ldb)
    !
    INTEGER :: ig
    !
    !
#ifdef __OPENMP
!$OMP PARALLEL private(ig)
!$OMP DO
#endif
    DO ig=1,ldb 
       b(ig)= (0.0_DP,0.0_DP) 
    ENDDO
#ifdef __OPENMP
!$OMP ENDDO
!$OMP DO
#endif
    DO ig=1,n
       b(nl (ig))=        a1(ig)
       b(nlm(ig))= DCONJG(a1(ig))
    ENDDO
#ifdef __OPENMP
!$OMP ENDDO
!$OMP END PARALLEL
#endif
    !
    CALL invfft('Dense', b, dfftp)
    !
  END SUBROUTINE
  !
  !
  SUBROUTINE smooth_invfft(n,a1,lda,b,ldb)
    !
    ! INVFFT : G ---> R
    !
    ! INPUT  : n     = actual number of PW
    !          a1     = ONE COMPLEX arrays containing ONE COMPLEX functions in G space
    !          lda   = leading dimension of a1
    !          ldb   = leading dimension of b          
    ! OUTPUT : b     = ONE COMPLEX array containing ONE REAL functions in R space + 0
    !
    IMPLICIT NONE
    !
    INTEGER,INTENT(IN) :: lda, n, ldb
    COMPLEX(DP),INTENT(IN) :: a1(lda)
    COMPLEX(DP),INTENT(OUT) :: b(ldb)
    !
    INTEGER :: ig
    !
    !
#ifdef __OPENMP
!$OMP PARALLEL private(ig)
!$OMP DO
#endif
    DO ig=1,ldb 
       b(ig)= (0.0_DP,0.0_DP) 
    ENDDO
#ifdef __OPENMP
!$OMP ENDDO
!$OMP DO
#endif
    DO ig=1,n
       b(nls (ig))=        a1(ig)
       b(nlsm(ig))= DCONJG(a1(ig))
    ENDDO
#ifdef __OPENMP
!$OMP ENDDO
!$OMP END PARALLEL
#endif
    !
    CALL invfft('Smooth', b, dffts)
    !
  END SUBROUTINE
  !
  !
END MODULE
