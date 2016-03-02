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
SUBROUTINE do_setup
  !-----------------------------------------------------------------------
  !
  USE pwcom,                  ONLY : npw,nbnd,nkstot,nspin,nelec,nelup,neldw,et,wk,wg,lspinorb,domag,lsda, &
                                     isk,nks,two_fermi_energies
  USE fixed_occ,              ONLY : tfixed_occ,f_inp
  USE kinds,                  ONLY : DP
  USE mp,                     ONLY : mp_sum
  USE mp_global,              ONLY : intra_bgrp_comm
  USE io_global,              ONLY : stdout
  USE lsda_mod,               ONLY : current_spin,lsda
  USE constants,              ONLY : rytoev
  USE control_flags,          ONLY : gamma_only
  USE noncollin_module,       ONLY : noncolin,npol
  USE cell_base,              ONLY : omega,celldm,at
  USE fft_base,               ONLY : dfftp,dffts
  USE gvecs,                  ONLY : ngms_g, ngms
  USE gvect,                  ONLY : ngm_g, ngm
  USE gvecw,                  ONLY : ecutwfc
  USE io_push
  !
  IMPLICIT NONE
  !
  INTEGER :: auxi,ib
  REAL(DP) :: alat
  INTEGER :: iks,spin
  !
  CALL start_clock('do_setup')
  !
  ! INIT PW
  !
  CALL init_pw_arrays(nbnd)
  CALL set_iks_l2g()
  !
  !
  IF ( lsda ) THEN
     IF ( INT( nelup ) == 0 .AND. INT( neldw ) == 0 ) THEN
     !IF ( .NOT. two_fermi_energies ) THEN
        DO iks = 1, nks
           spin = isk(iks)
           !
           SELECT CASE(spin)
           CASE(1)
              nelup = SUM( f_inp(:,1) )
           CASE(2)
              neldw = SUM( f_inp(:,2) )
           END SELECT
           !
        ENDDO
     ENDIF
     IF ( INT( nelup ) == 0 .AND. INT( neldw ) == 0 ) THEN
        CALL errore( 'do_setup','nelup = 0 and neldw = 0 ',1)
     ENDIF
  ENDIF
  !
  ! SYSTEM OVERVIEW
  !
  CALL io_push_title('System Overview')
  CALL io_push_value('gamma_only',gamma_only,20)
  CALL io_push_value('ecutwfc [Ry]',ecutwfc,20)
  CALL io_push_es0('omega [au^3]',omega,20)
  auxi = npw
  CALL mp_sum(auxi,intra_bgrp_comm)
  CALL io_push_value('glob. #G',auxi,20)
  CALL io_push_value('nbnd',nbnd,20)
  CALL io_push_value('nkstot',nkstot,20)
  CALL io_push_value('nspin',nspin,20)
  CALL io_push_value('nelec',nelec,20)
  IF(nspin == 2) THEN
     CALL io_push_value('nelup',nelup,20)
     CALL io_push_value('neldw',neldw,20)
  ENDIF
  CALL io_push_value('npol',npol,20)
  CALL io_push_value('lsda',lsda,20)
  CALL io_push_value('noncolin',noncolin,20)
  CALL io_push_value('lspinorb',lspinorb,20)
  CALL io_push_value('domag',domag,20)
  CALL io_push_bar
  !
  alat = celldm(1)
  !
  WRITE( stdout, '(/5x,"sFFT G-space: ",i8," G-vectors", 5x, &
       &               "R-space: (",i4,",",i4,",",i4,")")') &
       &         ngms_g, dffts%nr1, dffts%nr2, dffts%nr3
  WRITE( stdout, '( 5x,"dFFT G-space: ",i8," G-vectors", 5x, &
       &               "R-space: (",i4,",",i4,",",i4,")")') &
       &         ngm_g, dfftp%nr1, dfftp%nr2, dfftp%nr3
  WRITE( stdout, '(/5x,"Cell [a.u.]          = ",3f14.6)') alat*at(1,1:3)
  WRITE( stdout, '( 5x,"                     = ",3f14.6)') alat*at(2,1:3)
  WRITE( stdout, '( 5x,"                     = ",3f14.6)') alat*at(3,1:3)
  WRITE( stdout, '( 5x," ")')
  !
  CALL stop_clock('do_setup')
  !
END SUBROUTINE 
