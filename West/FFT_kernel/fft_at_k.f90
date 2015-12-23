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
MODULE fft_at_k
  !-----------------------------------------------------------------------
  !
  ! Everything is done following dffts
  !
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : nl,nlm
  USE gvecs,                ONLY : nls,nlsm
  USE fft_interfaces,       ONLY : fwfft,invfft
  USE fft_base,             ONLY : dffts
  USE pwcom,                ONLY : igk
  !
  IMPLICIT NONE
  !
  !
  CONTAINS
  !
  !
  SUBROUTINE singleband_invfft_k(n,a1,lda,b,ldb,luseigk)
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
    LOGICAL,INTENT(IN) :: luseigk
    !
    INTEGER :: ig
    !
    !
    b=0.0_DP
    IF(luseigk) THEN
       DO ig=1,n
          b(nls(igk(ig)))=a1(ig)
       ENDDO
    ELSE
       DO ig=1,n
          b(nls(ig))=a1(ig)
       ENDDO
    ENDIF
    !
    CALL invfft('Wave', b, dffts)
    !
  END SUBROUTINE
  !
  !
  SUBROUTINE singleband_fwfft_k(n,a,lda,b1,ldb,luseigk)
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
    LOGICAL,INTENT(IN) :: luseigk
    !
    INTEGER :: ig
    !
    !
    ! FFT call
    CALL fwfft('Wave',a, dffts)
    !
    b1=0.0_DP
    !
    IF(luseigk) THEN 
       DO ig=1,n 
          b1(ig) = a(nls(igk(ig)))
       ENDDO
    ELSE
       DO ig=1,n 
          b1(ig) = a(nls(ig))
       ENDDO
    ENDIF
    !
  END SUBROUTINE
  !
  !
END MODULE
