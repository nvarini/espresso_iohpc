!
! Copyright (C) 2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE pw_restart_new
!----------------------------------------------------------------------------
#ifdef __XSD
  !
  ! ... New PWscf I/O using xml schema and hdf5 binaries
  !
  USE qes_module
  USE qexsd_module, ONLY: qexsd_init_schema, qexsd_openschema, qexsd_closeschema, &
                          qexsd_init_convergence_info, qexsd_init_algorithmic_info, & 
                          qexsd_init_atomic_species, qexsd_init_atomic_structure, &
                          qexsd_init_symmetries, qexsd_init_basis_set, qexsd_init_dft, &
                          qexsd_init_magnetization,qexsd_init_band_structure,  &
                          qexsd_init_total_energy,qexsd_init_forces,qexsd_init_stress, &
                          qexsd_init_outputElectricField
  USE iotk_module
  USE io_global, ONLY : ionode, ionode_id
  USE io_files,  ONLY : iunpun, xmlpun_schema, prefix, tmp_dir, &
       delete_if_present
  USE pw_restart,ONLY : gk_l2gmap_kdip

  IMPLICIT NONE
  !
  LOGICAL :: lcell_read   = .FALSE., &
             lpw_read     = .FALSE., &
             lions_read   = .FALSE., &
             lspin_read   = .FALSE., &
             lstarting_mag_read   = .FALSE., &
             lxc_read     = .FALSE., &
             locc_read    = .FALSE., &
             lbz_read     = .FALSE., &
             lbs_read     = .FALSE., &
             lefield_read = .FALSE., &
             lwfc_read    = .FALSE., &
             lsymm_read   = .FALSE.
  !
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  PRIVATE
  PUBLIC :: pw_write_schema, pw_write_binaries, &
       pw_readschema_file, init_vars_from_schema, read_collected_to_evc
  !
  CONTAINS
    !------------------------------------------------------------------------
    SUBROUTINE pw_write_schema( )
      !------------------------------------------------------------------------
      !
      USE control_flags,        ONLY : istep, twfcollect, conv_ions, &
                                       lscf, gamma_only, &
                                       tqr, tq_smoothing, tbeta_smoothing, &
                                       noinv, do_makov_payne, smallmem, &
                                       llondon, lxdm, ts_vdw, scf_error, n_scf_steps
      USE constants,            ONLY : e2
      USE realus,               ONLY : real_space
      USE uspp,                 ONLY : okvan
      USE paw_variables,        ONLY : okpaw
      USE uspp_param,           ONLY : upf
      USE global_version,       ONLY : version_number
      USE cell_base,            ONLY : at, bg, alat, tpiba, tpiba2, &
                                       ibrav, celldm
      USE gvect,                ONLY : ig_l2g
      USE ions_base,            ONLY : nsp, ityp, atm, nat, tau, if_pos
      USE noncollin_module,     ONLY : noncolin, npol
      USE io_files,             ONLY : nwordwfc, iunwfc, psfile
      USE buffers,              ONLY : get_buffer
      USE wavefunctions_module, ONLY : evc
      USE klist,                ONLY : nks, nkstot, xk, ngk, wk, qnorm, &
                                       lgauss, ngauss, degauss, nelec, &
                                       two_fermi_energies, nelup, neldw
      USE start_k,              ONLY : nk1, nk2, nk3, k1, k2, k3, &
                                       nks_start, xk_start, wk_start
      USE ktetra,               ONLY : ntetra, tetra, ltetra
      USE gvect,                ONLY : ngm, ngm_g, g, mill
      USE fft_base,             ONLY : dfftp
      USE basis,                ONLY : natomwfc
      USE gvecs,                ONLY : ngms_g, dual
      USE fft_base,             ONLY : dffts
      USE wvfct,                ONLY : npwx, et, wg, nbnd
      USE ener,                 ONLY : ef, ef_up, ef_dw, vtxc, etxc, ewld, etot, &
                                       ehart, eband, demet 
      USE gvecw,                ONLY : ecutwfc
      USE fixed_occ,            ONLY : tfixed_occ, f_inp
      USE ldaU,                 ONLY : lda_plus_u, lda_plus_u_kind, U_projection, &
                                       Hubbard_lmax, Hubbard_l, Hubbard_U, Hubbard_J, &
                                       Hubbard_alpha, Hubbard_J0, Hubbard_beta,&
                                       is_hubbard
      USE spin_orb,             ONLY : lspinorb, domag
      USE symm_base,            ONLY : nrot, nsym, invsym, s, ft, irt, &
                                       t_rev, sname, time_reversal, no_t_rev
      USE lsda_mod,             ONLY : nspin, isk, lsda, starting_magnetization, magtot, absmag
      USE noncollin_module,     ONLY : angle1, angle2, i_cons, mcons, bfield, magtot_nc, &
                                       lambda
      USE ions_base,            ONLY : amass
      USE funct,                ONLY : get_dft_name, get_inlc, get_nonlocc_name, dft_is_nonlocc
      USE kernel_table,         ONLY : vdw_table_name
      USE scf,                  ONLY : rho
      USE force_mod,            ONLY : lforce, sumfor, force, sigma, lstres
      USE extfield,             ONLY : tefield, dipfield, edir, etotefield, &
                                       emaxpos, eopreg, eamp, &
                                       monopole, zmon, relaxz, block, block_1,&
                                       block_2, block_height ! TB
      USE io_rho_xml,           ONLY : write_rho
      USE mp,                   ONLY : mp_sum
      USE mp_world,             ONLY : nproc
      USE mp_images,            ONLY : nproc_image
      USE mp_pools,             ONLY : kunit, nproc_pool, me_pool, root_pool, &
                                       intra_pool_comm, inter_pool_comm
      USE mp_bands,             ONLY : nproc_bgrp, me_bgrp, root_bgrp, &
                                       intra_bgrp_comm, nbgrp, ntask_groups
      USE mp_diag,              ONLY : nproc_ortho
      USE funct,                ONLY : get_exx_fraction, dft_is_hybrid, &
                                       get_gau_parameter, &
                                       get_screening_parameter, exx_is_active
      USE exx,                  ONLY : x_gamma_extrapolation, nq1, nq2, nq3, &
                                       exxdiv_treatment, yukawa, ecutvcut, ecutfock
      USE cellmd,               ONLY : lmovecell, cell_factor 
      USE martyna_tuckerman,    ONLY : do_comp_mt
      USE esm,                  ONLY : do_comp_esm, esm_nfit, esm_efield, esm_w, &
                                       esm_a, esm_bc
      USE london_module,        ONLY : scal6, lon_rcut, in_c6
      USE xdm_module,           ONLY : xdm_a1=>a1i, xdm_a2=>a2i
      USE tsvdw_module,         ONLY : vdw_isolated, vdw_econv_thr
      USE input_parameters,     ONLY : space_group, verbosity, calculation, ion_dynamics, starting_ns_eigenvalue, &
                                       vdw_corr, london, input_parameters_occupations => occupations
      USE bp,                   ONLY : lelfield, lberry, bp_mod_el_pol => el_pol, bp_mod_ion_pol => ion_pol
      !
      USE rap_point_group,      ONLY : elem, nelem, name_class
      USE rap_point_group_so,   ONLY : elem_so, nelem_so, name_class_so
      USE bfgs_module,          ONLY : bfgs_get_n_iter
      USE qexsd_module,         ONLY : qexsd_dipol_obj, qexsd_bp_obj
     
      !
      IMPLICIT NONE
      !
      CHARACTER(15)         :: subname="pw_write_schema"
      CHARACTER(LEN=20)     :: dft_name
      CHARACTER(LEN=256)    :: dirname
      CHARACTER(LEN=80)     :: vdw_corr_
      INTEGER               :: i, ig, ik, ngg, ierr, ipol, num_k_points
      INTEGER               :: npwx_g, npw_g, ispin, inlc
      INTEGER,  ALLOCATABLE :: ngk_g(:)
      INTEGER,  ALLOCATABLE :: igk_l2g(:,:), igk_l2g_kdip(:,:), mill_g(:,:)
      LOGICAL               :: lwfc, lrho, lxsd, occupations_are_fixed
      CHARACTER(iotk_attlenx)  :: attr
      INTEGER                  :: iclass, isym, ielem
      CHARACTER(LEN=15)        :: symop_2_class(48)
      LOGICAL                  :: opt_conv_ispresent
      INTEGER                  :: n_opt_steps, h_band
      REAL(DP)                 :: h_energy
      !
      TYPE(output_type) :: output
      
      !
      ! PW dimensions need to be properly computed 
      ! reducing across MPI tasks
      !
      ALLOCATE( ngk_g( nkstot ) )
      !
      ngk_g(1:nks) = ngk(:)
      CALL mp_sum( ngk_g(1:nks), intra_bgrp_comm )
      CALL ipoolrecover( ngk_g, 1, nkstot, nks )
      !
      ! ... compute the maximum number of G vector among all k points
      !
      npwx_g = MAXVAL( ngk_g(1:nkstot) )
      !
      !DEALLOCATE( ngk_g )
      ! do not deallocate ngk_g here, please, I need it for band_structure_init
      ! P. Delugas 
      !
      ! ... find out the global number of G vectors: ngm_g
      !
      ngm_g = ngm
      CALL mp_sum( ngm_g, intra_bgrp_comm )
      ! 
      IF (tefield .AND. dipfield ) THEN 
          CALL init_dipole_info(qexsd_dipol_obj, rho%of_r)   
          qexsd_dipol_obj%tagname = "dipoleInfo"
      END IF
      ! 
      !
      ! XML descriptor
      ! 
      dirname = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
      !
      CALL qexsd_init_schema( iunpun )
      !
      !
      IF ( ionode ) THEN  
         !
         ! ... here we init the variables and finally write them to file
         !
!-------------------------------------------------------------------------------
! ... HEADER
!-------------------------------------------------------------------------------
         !
         CALL qexsd_openschema(TRIM( dirname ) // '/' // TRIM( xmlpun_schema ))
         output%tagname="output"
         !
!-------------------------------------------------------------------------------
! ... CONVERGENCE_INFO
!-------------------------------------------------------------------------------
         SELECT CASE (TRIM( calculation )) 
            CASE ( "relax","vc-relax" ,"md")
                opt_conv_ispresent = .TRUE.
                IF (TRIM( ion_dynamics) == 'bfgs' ) THEN 
                    n_opt_steps = bfgs_get_n_iter('bfgs_iter ') 
                ELSE 
                    n_opt_steps = istep 
                END IF 
            CASE default
                opt_conv_ispresent = .FALSE.
                n_opt_steps        = 0 
         END SELECT
         ! 
         call qexsd_init_convergence_info(output%convergence_info, &
              n_scf_steps=n_scf_steps, scf_error=scf_error, &
              opt_conv_ispresent=lforce, &
              n_opt_steps=n_opt_steps, grad_norm=sumfor )
         !
!-------------------------------------------------------------------------------
! ... ALGORITHMIC_INFO
!-------------------------------------------------------------------------------
         !
         CALL qexsd_init_algorithmic_info(output%algorithmic_info, &
              real_space_q=real_space, uspp=okvan, paw=okpaw)
         !
!-------------------------------------------------------------------------------
! ... ATOMIC_SPECIES
!-------------------------------------------------------------------------------
         !
         ! while amass's are always present, starting_mag should not be passed
         ! for nspin==1 or contrained magnetization calculations
         !
         IF (noncolin) THEN
            CALL qexsd_init_atomic_species(output%atomic_species, nsp, atm, psfile, &
                 amass, STARTING_MAGNETIZATION = starting_magnetization, &
                 ANGLE1=angle1, ANGLE2=angle2)
         ELSE IF (nspin==2) THEN 
            CALL qexsd_init_atomic_species(output%atomic_species, nsp, atm, psfile, &
                 amass, STARTING_MAGNETIZATION=starting_magnetization)
         ELSE 
            CALL qexsd_init_atomic_species(output%atomic_species, nsp, atm,psfile, &
                 amass)
         END IF
         !
!-------------------------------------------------------------------------------
! ... ATOMIC_STRUCTURE
!-------------------------------------------------------------------------------
         !         
         CALL qexsd_init_atomic_structure(output%atomic_structure, nsp, atm, ityp, &
              nat, tau, 'Bohr', alat, alat*at(:,1), alat*at(:,2), alat*at(:,3), ibrav)
         !
!-------------------------------------------------------------------------------
! ... SYMMETRIES
!-------------------------------------------------------------------------------
         !
         symop_2_class="not found"
         IF (TRIM (verbosity) == 'medium' .OR. TRIM(verbosity) == 'high') THEN
            IF ( noncolin )  THEN 
               symmetries_so_loop:DO isym = 1, nrot 
                  classes_so_loop:DO iclass = 1, 24
                     elements_so_loop:DO ielem=1, nelem_so(iclass)
                        IF ( elem_so(ielem,iclass) == isym) THEN 
                           symop_2_class(isym) = name_class_so(iclass)
                           EXIT symmetries_so_loop
                        END IF
                     END DO elements_so_loop 
                  END DO classes_so_loop
               END DO symmetries_so_loop
            !
            ELSE
               symmetries_loop:DO isym = 1, nrot
                  classes_loop:DO iclass = 1, 12
                     elements_loop:DO ielem=1, nelem (iclass)
                        IF ( elem(ielem,iclass) == isym) THEN
                           symop_2_class(isym) = name_class(iclass)
                           EXIT classes_loop
                        END IF
                     END DO elements_loop
                  END DO classes_loop
               END DO symmetries_loop
            END IF
         END IF
         CALL qexsd_init_symmetries(output%symmetries, nsym, nrot, space_group, &
              s, ft, sname, t_rev, nat, irt,symop_2_class(1:nrot), verbosity, &
              noncolin)
         !
!-------------------------------------------------------------------------------
! ... BASIS SET
!-------------------------------------------------------------------------------
         !
         CALL qexsd_init_basis_set(output%basis_set, gamma_only, ecutwfc/e2, ecutwfc*dual/e2, &
              dfftp%nr1, dfftp%nr2, dfftp%nr3, dffts%nr1, dffts%nr2, dffts%nr3, &
              .FALSE., dfftp%nr1, dfftp%nr2, dfftp%nr3, ngm_g, ngms_g, npwx_g, &
              bg(:,1), bg(:,2), bg(:,3) )
         !
!-------------------------------------------------------------------------------
! ... DFT
!-------------------------------------------------------------------------------
         !
         dft_name = get_dft_name()
         inlc = get_inlc()
         !
         IF ( lda_plus_u .AND. noncolin) CALL errore(subname,"LDA+U and non-collinear case not implemented in qexsd",10)
         !
         vdw_corr_ = vdw_corr
         IF ( london ) vdw_corr_ = 'grimme-d2'
         CALL qexsd_init_dft(output%dft, dft_name, .TRUE., &
              dft_is_hybrid(), nq1, nq2, nq3, ecutfock, &
              get_exx_fraction(), get_screening_parameter(), exxdiv_treatment, &
              x_gamma_extrapolation, ecutvcut, lda_plus_u, lda_plus_u_kind, 2*Hubbard_lmax+1, &
              nspin, nsp, 2*Hubbard_lmax+1, nat, atm, ityp, Hubbard_U, Hubbard_J0, Hubbard_alpha, &
              Hubbard_beta, Hubbard_J, starting_ns_eigenvalue, rho%ns, rho%ns_nc, U_projection, &
              dft_is_nonlocc(), TRIM(vdw_corr_), TRIM ( get_nonlocc_name()), scal6, in_c6, lon_rcut, xdm_a1, xdm_a2,&
              vdw_econv_thr, vdw_isolated, is_hubbard, upf(1:nsp)%psd)
         !
!-------------------------------------------------------------------------------
! ... MAGNETIZATION
!-------------------------------------------------------------------------------
         !
         CALL qexsd_init_magnetization(output%magnetization, lsda, noncolin, lspinorb, &
              magtot, magtot_nc, absmag, domag )
         !

!--------------------------------------------------------------------------------------
! ... BAND STRUCTURE
!-------------------------------------------------------------------------------------
         !
         IF (TRIM(input_parameters_occupations) == 'fixed') THEN 
            occupations_are_fixed = .TRUE. 
            IF ( nspin == 1 ) THEN 
               h_band = NINT ( nelec/2.d0) 
            ELSE 
               h_band = NINT ( nelec ) 
            END IF  
            h_energy =MAXVAL (et(h_band, 1:nkstot))
         ELSE 
            occupations_are_fixed = .FALSE. 
            h_energy  = ef 
         END IF 
         CALL  qexsd_init_band_structure(output%band_structure,lsda,noncolin,lspinorb, &
              nbnd,nelec, natomwfc, occupations_are_fixed, & 
              h_energy,two_fermi_energies, [ef_up,ef_dw], et,wg,nkstot,xk,ngk_g,wk)
         !
!-------------------------------------------------------------------------------------------
! ... TOTAL ENERGY
!-------------------------------------------------------------------------------------------
         !
         IF (tefield) THEN
            CALL  qexsd_init_total_energy(output%total_energy,etot,eband,ehart,vtxc,etxc, &
                 ewld,degauss,demet, etotefield)
         ELSE 
            CALL  qexsd_init_total_energy(output%total_energy,etot,eband,ehart,vtxc,etxc, &
                 ewld,degauss,demet)
         END IF
         !
!---------------------------------------------------------------------------------------------
! ... FORCES
!----------------------------------------------------------------------------------------------
         !
         IF ( lforce ) THEN 
            output%forces_ispresent = .TRUE.
            CALL qexsd_init_forces(output%forces,nat,force,lforce)
         ELSE 
            output%forces_ispresent = .FALSE.
            output%forces%lwrite = .FALSE.  
         END IF 
         !
!------------------------------------------------------------------------------------------------
! ... STRESS 
!------------------------------------------------------------------------------------------------
         IF ( lstres) THEN
            output%stress_ispresent=.TRUE.
            CALL qexsd_init_stress(output%stress, sigma, lstres ) 
         ELSE 
            output%stress_ispresent=.FALSE.
            output%stress%lwrite=.FALSE.
         END IF
!-------------------------------------------------------------------------------------------------
! ... ELECTRIC FIELD
!-------------------------------------------------------------------------------------------------
         IF ( lelfield ) THEN
            output%electric_field_ispresent = .TRUE. 
            CALL qexsd_init_outputElectricField(output%electric_field, lelfield, tefield, dipfield, &
                 lberry, el_pol = bp_mod_el_pol, ion_pol = bp_mod_ion_pol) 
         ELSE IF ( lberry ) THEN 
            output%electric_field_ispresent = .TRUE.
            CALL qexsd_init_outputElectricField(output%electric_field, lelfield, tefield, dipfield, & 
                 lberry, bp_obj=qexsd_bp_obj) 
         ELSE IF ( tefield .AND. dipfield  ) THEN 
            output%electric_field_ispresent = .TRUE.
            CALL  qexsd_init_outputElectricField(output%electric_field, lelfield, tefield, dipfield, &
                 lberry, dipole_obj = qexsd_dipol_obj )                     
         ELSE 
            output%electric_field_ispresent = .FALSE.
         ENDIF
!------------------------------------------------------------------------------------------------
! ... ACTUAL WRITING
!-------------------------------------------------------------------------------
         !
         CALL qes_write_output(iunpun,output)
         CALL qes_reset_output(output)
         !
!-------------------------------------------------------------------------------
! ... CLOSING
!-------------------------------------------------------------------------------
         !
         CALL qexsd_closeschema()
         !
      END IF
      DEALLOCATE (ngk_g)
      !
      RETURN
       !
    END SUBROUTINE pw_write_schema
    !
    !------------------------------------------------------------------------
    SUBROUTINE pw_write_binaries( )
      !------------------------------------------------------------------------
      !
      USE iotk_module
      USE mp,                   ONLY : mp_bcast, mp_sum, mp_max
      USE io_global,            ONLY : ionode, ionode_id
      USE xml_io_base,          ONLY : write_wfc
      USE control_flags,        ONLY : gamma_only, smallmem
      USE gvect,                ONLY : ig_l2g
      USE noncollin_module,     ONLY : noncolin, npol

      USE buffers,              ONLY : get_buffer
      USE wavefunctions_module, ONLY : evc
      USE klist,                ONLY : nks, nkstot, xk, ngk, igk_k, wk
      USE gvect,                ONLY : ngm, ngm_g, g, mill
      USE fft_base,             ONLY : dfftp
      USE basis,                ONLY : natomwfc
      USE gvecs,                ONLY : ngms_g, dual
      USE wvfct,                ONLY : npw, npwx, et, wg, nbnd
      USE lsda_mod,             ONLY : nspin, isk, lsda
      USE mp_world,             ONLY : nproc
      USE mp_images,            ONLY : nproc_image, intra_image_comm
      USE mp_pools,             ONLY : kunit, nproc_pool, me_pool, root_pool, &
                                       intra_pool_comm, inter_pool_comm
      USE mp_bands,             ONLY : nproc_bgrp, me_bgrp, root_bgrp, &
                                       intra_bgrp_comm, nbgrp, ntask_groups
      USE mp_diag,              ONLY : nproc_ortho
      !
      IMPLICIT NONE
      !
      INTEGER               :: i, ig, ik, ngg, ierr, ipol, num_k_points
      INTEGER               :: nkl, nkr, npwx_g
      INTEGER               :: ike, iks, npw_g, ispin, inlc
      INTEGER,  ALLOCATABLE :: ngk_g(:)
      INTEGER,  ALLOCATABLE :: igk_l2g(:,:), igk_l2g_kdip(:,:), mill_g(:,:)
      CHARACTER(LEN=256)    :: dirname
      CHARACTER(LEN=320)    :: filename
      CHARACTER(iotk_attlenx)  :: attr
      !
      !
      dirname = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
      !
      IF ( nspin == 2 ) THEN
         num_k_points = nkstot / 2
      ELSE
         num_k_points = nkstot
      END IF
      !
      CALL  kpoint_global_indices (nkstot, iks, ike)
      !
      ! ... find out the global number of G vectors: ngm_g
      !
      ngm_g = ngm
      !
      CALL mp_sum( ngm_g, intra_bgrp_comm )
      !
      ! ... collect all G-vectors across processors within the pools
      !
      ALLOCATE( mill_g( 3, ngm_g ) )
      !
      mill_g = 0
      !
      DO ig = 1, ngm
         !
         mill_g(1,ig_l2g(ig)) = mill(1,ig)
         mill_g(2,ig_l2g(ig)) = mill(2,ig)
         mill_g(3,ig_l2g(ig)) = mill(3,ig)
         !
      END DO
      !
      CALL mp_sum( mill_g, intra_bgrp_comm )
      !
      ! ... build the igk_l2g array, yielding the correspondence between
      ! ... the local k+G index and the global G index, from previously
      ! ... computed arrays igk_k (k+G indices) and ig_l2g
      ! ... Beware: for variable-cell case, one has to use starting G and 
      ! ... k+G vectors
      !
      ALLOCATE ( igk_l2g( npwx, nks ) )
      igk_l2g = 0
      !
      DO ik = 1, nks
         DO ig = 1, ngk (ik)
            igk_l2g(ig,ik) = ig_l2g(igk_k(ig,ik))
         END DO
      END DO
      !
      ! ... compute the global number of G+k vectors for each k point
      !
      ALLOCATE( ngk_g( nkstot ) )
      !
      ngk_g = 0
      ngk_g(iks:ike) = ngk(1:nks)
      !
      CALL mp_sum( ngk_g, inter_pool_comm)
      CALL mp_sum( ngk_g, intra_pool_comm)
      !
      ngk_g = ngk_g / nbgrp
      !
      ! ... compute the maximum G vector index among all G+k an processors
      !
      npw_g = MAXVAL( igk_l2g(:,:) )
      !
      CALL mp_max( npw_g, inter_pool_comm )
      CALL mp_max( npw_g, intra_pool_comm )
      !
      ! ... compute the maximum number of G vector among all k points
      !
      npwx_g = MAXVAL( ngk_g(1:nkstot) )
      !
      ! ... define a further l2g map to write gkvectors and wfc coherently
      !
      ALLOCATE ( igk_l2g_kdip( npwx_g, nks ) )
      !
      igk_l2g_kdip = 0
      !
      DO ik = iks, ike
         !
         CALL gk_l2gmap_kdip( npw_g, ngk_g(ik), ngk(ik-iks+1), &
                              igk_l2g(1,ik-iks+1), igk_l2g_kdip(1,ik-iks+1) )
      END DO
      !
      IF ( ionode ) THEN  
         !
         ! ... write the G-vectors - backwards compatibility
         !
         filename = TRIM( dirname ) // '/gvectors.dat'
         CALL write_gvecs(iunpun, filename)
         !
      END IF
      !
      k_points_loop2: DO ik = 1, num_k_points
         !
         IF ( ionode ) filename = TRIM(dirname)//'/gkvectors'//TRIM(int_to_char(ik))//'.dat'
         IF (.NOT.smallmem) CALL write_gk( iunpun, ik, filename )
         !
         CALL write_this_wfc ( iunpun, ik )
         !
      END DO k_points_loop2
      !
      DEALLOCATE ( igk_l2g )
      DEALLOCATE ( igk_l2g_kdip )
      !
!-------------------------------------------------------------------------------
! ... END RESTART SECTIONS
!-------------------------------------------------------------------------------
      !
      DEALLOCATE( mill_g )
      DEALLOCATE( ngk_g )
      !
      RETURN
      !
      CONTAINS
        !
        !--------------------------------------------------------------------
        SUBROUTINE write_gvecs( iun, filename )
          !--------------------------------------------------------------------
          IMPLICIT NONE
          !
          INTEGER,            INTENT(IN) :: iun
          CHARACTER(LEN=*),   INTENT(IN) :: filename
          !
          CALL iotk_open_write( iun, FILE = TRIM( filename ), &
               BINARY = .true. )
          !
          CALL iotk_write_begin( iun, "G-VECTORS" )
          CALL iotk_write_attr( attr, "nr1s", dfftp%nr1, FIRST = .true. )
          CALL iotk_write_attr( attr, "nr2s", dfftp%nr2 )
          CALL iotk_write_attr( attr, "nr3s", dfftp%nr3 )
          CALL iotk_write_attr( attr, "gvect_number", ngm_g )
          CALL iotk_write_attr( attr, "gamma_only", gamma_only )
          CALL iotk_write_attr( attr, "units", "crystal" )
          CALL iotk_write_empty( iun, "INFO", ATTR = attr )
          !
          CALL iotk_write_dat  ( iun, "g", mill_g(1:3,1:ngm_g), COLUMNS=3)
          CALL iotk_write_end  ( iun, "G-VECTORS" )
          !
          CALL iotk_close_write( iun )
          !
        END SUBROUTINE write_gvecs
        !
        !--------------------------------------------------------------------
        SUBROUTINE write_gk( iun, ik, filename )
          !--------------------------------------------------------------------
          !
          IMPLICIT NONE
          !
          INTEGER,            INTENT(IN) :: iun, ik
          CHARACTER(LEN=*),   INTENT(IN) :: filename
          !
          INTEGER, ALLOCATABLE :: igwk(:,:)
          INTEGER, ALLOCATABLE :: itmp(:)
          !
          !
          ALLOCATE( igwk( npwx_g, nkstot ) )
          !
          igwk(:,ik) = 0
          !
          ALLOCATE( itmp( npw_g ) )
          !
          itmp = 0
          !
          IF ( ik >= iks .AND. ik <= ike ) THEN
             !
             DO ig = 1, ngk(ik-iks+1)
                !
                itmp(igk_l2g(ig,ik-iks+1)) = igk_l2g(ig,ik-iks+1)
                !
             END DO
             !
          END IF
          !
          CALL mp_sum( itmp, inter_pool_comm )
          CALL mp_sum( itmp, intra_pool_comm )
          !
          ngg = 0
          !
          DO ig = 1, npw_g
             !
             if ( itmp(ig) == ig ) THEN
                !
                ngg = ngg + 1
                !
                igwk(ngg,ik) = ig
                !
             END IF
             !
          END DO
          !
          DEALLOCATE( itmp )
          !
          IF ( ionode ) THEN
             !
             CALL iotk_open_write( iun, FILE = TRIM( filename ), &
                                   BINARY = .TRUE. )
             !
             CALL iotk_write_dat( iun, "NUMBER_OF_GK-VECTORS", ngk_g(ik) )
             CALL iotk_write_dat( iun, "MAX_NUMBER_OF_GK-VECTORS", npwx_g )
             CALL iotk_write_dat( iun, "GAMMA_ONLY", gamma_only )
             !
             CALL iotk_write_attr ( attr, "UNITS", "2 pi / a", FIRST = .TRUE. )
             CALL iotk_write_dat( iun, "K-POINT_COORDS", xk(:,ik), ATTR = attr )
             !
             CALL iotk_write_dat( iun, "INDEX", igwk(1:ngk_g(ik),ik) )
             CALL iotk_write_dat( iun, "GRID", mill_g(1:3,igwk(1:ngk_g(ik),ik)), &
                                  COLUMNS = 3 )
             !
             CALL iotk_close_write( iun )
             !
          END IF
          !
          DEALLOCATE( igwk )
          !
        END SUBROUTINE write_gk
        !
        !
        !--------------------------------------------------------------------
        SUBROUTINE write_this_wfc ( iun, ik )
          !--------------------------------------------------------------------
          !
          USE io_files, ONLY : iunwfc, nwordwfc
          IMPLICIT NONE
          !
          INTEGER, INTENT(IN) :: iun, ik
          CHARACTER(LEN=256)  :: filename
          INTEGER :: ispin,ik_eff
          !
          ! ... wavefunctions - do not read if already in memory (nsk==1)
          ! ...                 read only if on this pool (iks <= ik <= ike )
          !
          IF ( ( nks > 1 ) .AND. ( ik >= iks ) .AND. ( ik <= ike ) ) THEN
             CALL get_buffer ( evc, nwordwfc, iunwfc, (ik-iks+1) )
          END IF
          !
          IF ( nspin == 2 ) THEN
             !
             ! ... beware: with pools, isk(ik) has the correct value for 
             ! ... all k-points only on first pool (ionode proc is ok)
             !
             ispin = isk(ik)
             !
             IF ( ionode ) filename = TRIM(dirname)//'/wfc'//TRIM(int_to_char(ik))//'dat'
             !
             CALL write_wfc( iun, ik, nkstot, kunit, ispin, nspin, &
                  evc, npw_g, gamma_only, nbnd, igk_l2g_kdip(:,ik-iks+1),   &
                  ngk(ik-iks+1), filename, 1.D0, &
                  ionode, root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm )
             !
             ik_eff = ik + num_k_points
             !
             ispin = isk(ik_eff)
             !
             ! ... LSDA: now read minority wavefunctions (if not already
             ! ... in memory and if they are on this pool)
             !
             IF ( ( nks > 1 ) .AND. ( ik_eff >= iks ) .AND. ( ik_eff <= ike ) ) THEN
                !
                CALL get_buffer ( evc, nwordwfc, iunwfc, (ik_eff-iks+1) )
                !
             END IF
             !
             IF ( ionode ) THEN
                !
                IF ( ispin == 1 ) THEN
                   filename = TRIM(dirname)//'/wfcup'//TRIM(int_to_char(ik))//'.dat'
                ELSE
                   filename = TRIM(dirname)//'/wfcdw'//TRIM(int_to_char(ik))//'.dat'
                END IF
                !
             END IF
             !
             CALL write_wfc( iun, ik_eff, nkstot, kunit, ispin, nspin, &
                  evc, npw_g, gamma_only, nbnd, igk_l2g_kdip(:,ik_eff-iks+1), &
                  ngk(ik_eff-iks+1), filename, 1.D0, &
                  ionode, root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm )
             !
          ELSE
             !
             IF ( noncolin ) THEN
                !
                DO ipol = 1, npol
                   !
                   IF ( ipol == 1 ) THEN
                      filename = TRIM(dirname)//'/wfcup'//TRIM(int_to_char(ik))//'.dat'
                   ELSE
                      filename = TRIM(dirname)//'/wfcdw'//TRIM(int_to_char(ik))//'.dat'
                   END IF
                   !
                   ! TEMP  spin-up and spin-down spinor components are written
                   ! TEMP  to different files, like in LSDA - not a smart way
                   !
                   nkl=(ipol-1)*npwx+1
                   nkr= ipol   *npwx
                   CALL write_wfc( iun, ik, nkstot, kunit, ipol, npol,   &
                        evc(nkl:nkr,:), npw_g, gamma_only, nbnd, &
                        igk_l2g_kdip(:,ik-iks+1), ngk(ik-iks+1), &
                        filename, 1.D0, &
                        ionode, root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm )
                   !
                END DO
                !
             ELSE
                !
                ispin = 1
                !
                IF ( ionode ) filename = TRIM(dirname)//'/wfc'//TRIM(int_to_char(ik))//'.dat'

                !
                CALL write_wfc( iun, ik, nkstot, kunit, ispin, nspin, &
                     evc, npw_g, gamma_only, nbnd,            &
                     igk_l2g_kdip(:,ik-iks+1),                &
                     ngk(ik-iks+1), filename, 1.D0, &
                     ionode, root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm )
                !
             END IF
             !
          END IF
          !
        END SUBROUTINE write_this_wfc
        !
    END SUBROUTINE pw_write_binaries

    !------------------------------------------------------------------------
    SUBROUTINE pw_readschema_file(ierr, restart_output, restart_input, &
         restart_parallel_info, restart_general_info)
      !------------------------------------------------------------------------
      USE qes_types_module,     ONLY : input_type, output_type, general_info_type, parallel_info_type    
      USE qexsd_reader_module,  ONLY : qexsd_get_output, qexsd_get_input, qexsd_get_general_info, &
                                       qexsd_get_parallel_info
      !
      USE qes_libs_module,      ONLY : qes_write_input, qes_write_output, qes_write_parallel_info, &
                                       qes_write_general_info
      IMPLICIT NONE 
      ! 
      INTEGER                                            :: ierr, iotk_err  
      TYPE( output_type ),OPTIONAL,      INTENT(OUT)   :: restart_output
      TYPE(input_type),OPTIONAL,         INTENT(OUT)   :: restart_input
      TYPE(parallel_info_type),OPTIONAL, INTENT(OUT)   :: restart_parallel_info
      TYPE(general_info_type ),OPTIONAL, INTENT(OUT)   :: restart_general_info
      ! 
      LOGICAL                                   :: found
      CHARACTER(LEN=80)                         :: errmsg = "" 
      CHARACTER(LEN=256)                        :: dirname
      !  
      ! 
      ierr = 0 
      ! 
      dirname = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
      CALL iotk_free_unit( iunpun, iotk_err )
      !
      CALL errore( 'pw_readschema_file', &
                   'no free units to read xsd output', iotk_err )
      CALL qexsd_init_schema( iunpun )
      CALL iotk_open_read( iunpun, TRIM(dirname)//'/'//TRIM(xmlpun_schema))
      !
      IF ( PRESENT ( restart_general_info ) ) THEN 
         CALL qexsd_get_general_info ( iunpun, restart_general_info , found)
         IF (.NOT. found ) ierr = ierr + 1
         IF ( ierr /=0 ) THEN
            errmsg='error header of xml data file'
            GOTO 100
         END IF
         ! CALL qes_write_general_info( 82, restart_general_info) 
      END IF 
      ! 
      IF ( PRESENT ( restart_parallel_info ) ) THEN 
         CALL qexsd_get_parallel_info ( iunpun, restart_parallel_info, found ) 
         IF ( .NOT. found ) THEN 
            ierr = ierr + 1  
            errmsg='error parallel_info  of xsd data file' 
            GOTO 100
         END IF
         ! CALL qes_write_parallel_info ( 82, restart_parallel_info )
      END IF  
      ! 
      !
      IF ( PRESENT ( restart_input ) ) THEN   
         CALL qexsd_get_input ( iunpun, restart_input, found ) 
         IF ( .NOT. found ) THEN 
            ierr = ierr + 1  
            errmsg='error input of xsd data file' 
            GOTO 100
         END IF
         ! CALL qes_write_input( 82, restart_input )  
      END IF 
      ! 
      IF ( PRESENT ( restart_output ) ) THEN 
         CALL qexsd_get_output ( iunpun, restart_output, found ) 
         IF ( .NOT. found ) THEN 
            ierr = ierr + 1 
            errmsg = 'error output of xsd data file' 
            GOTO 100 
         END IF 
         ! CALL qes_write_output ( 82, restart_output ) 
      END IF 
      CALL iotk_close_read (iunpun)
      RETURN
 100  CALL errore('pw_readschemafile',TRIM(errmsg),ierr)
    END SUBROUTINE pw_readschema_file
    !  
    !------------------------------------------------------------------------
    SUBROUTINE init_vars_from_schema( what, ierr, output_obj, input_obj, par_info, gen_info )
      !------------------------------------------------------------------------
      !
      USE control_flags,        ONLY : twfcollect
      USE io_rho_xml,           ONLY : read_rho
      USE scf,                  ONLY : rho
      USE lsda_mod,             ONLY : nspin
      USE mp_world,             ONLY : mpime
      USE mp_bands,             ONLY : intra_bgrp_comm
      USE mp,                   ONLY : mp_sum, mp_barrier
      USE qes_types_module,     ONLY : input_type, output_type, general_info_type, parallel_info_type    
!      !
      IMPLICIT NONE
!      !
      CHARACTER(LEN=*), INTENT(IN)           :: what
      TYPE ( output_type), INTENT(IN)        :: output_obj
      TYPE ( input_type ), INTENT(IN)        :: input_obj
      TYPE ( parallel_info_type), INTENT(IN) :: par_info
      TYPE ( general_info_type ), INTENT(IN) :: gen_info
      INTEGER,INTENT (OUT)                   :: ierr 
      !
      CHARACTER(LEN=256) :: dirname
      LOGICAL            :: lcell, lpw, lions, lspin, linit_mag, &
                            lxc, locc, lbz, lbs, lwfc, lheader,          &
                            lsymm, lrho, lefield, ldim, &
                            lef, lexx, lesm
      !
      LOGICAL            :: need_qexml, found, electric_field_ispresent
      INTEGER            :: tmp, iotk_err 
      
      !    
      !
      ierr = 0 
      dirname = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
      !
      !
      !
      !
      ldim    = .FALSE.
      lcell   = .FALSE.
      lpw     = .FALSE.
      lions   = .FALSE.
      lspin   = .FALSE.
      linit_mag = .FALSE.
      lxc     = .FALSE.
      locc    = .FALSE.
      lbz     = .FALSE.
      lbs     = .FALSE.
      lwfc    = .FALSE.
      lsymm   = .FALSE.
      lrho    = .FALSE.
      lefield = .FALSE.
      lef     = .FALSE.
      lexx    = .FALSE.
      lesm    = .FALSE.
      !
     
         
      SELECT CASE( what )
      CASE( 'header' )
         !
         lheader = .TRUE.
         need_qexml = .TRUE.
         !
      CASE ( 'wf_collect' ) 
         ! 
         twfcollect = input_obj%control_variables%wf_collect
         !
      CASE( 'dim' )
         !
         ldim =       .TRUE.
         lbz  =       .TRUE.
         need_qexml = .TRUE.
         !
      CASE( 'pseudo' )
         !
         lions = .TRUE.
         need_qexml = .TRUE.
         !
      CASE( 'config' )
         !
         lcell = .TRUE.
         lions = .TRUE.
         need_qexml = .TRUE.
         !
      CASE( 'rho' )
         !
         lrho  = .TRUE.
         !
      CASE( 'wave' )
         !
         lpw   = .TRUE.
         lwfc  = .TRUE.
         need_qexml = .TRUE.
         !
      CASE( 'nowavenobs' )
         !
         lcell   = .TRUE.
         lpw     = .TRUE.
         lions   = .TRUE.
         lspin   = .TRUE.
         linit_mag   = .TRUE.
         lxc     = .TRUE.
         locc    = .TRUE.
         lbz     = .TRUE.
         lsymm   = .TRUE.
         lefield = .TRUE.
         need_qexml = .TRUE.
                                                                  
      CASE( 'nowave' )
         !
         lcell   = .TRUE.
         lpw     = .TRUE.
         lions   = .TRUE.
         lspin   = .TRUE.
         linit_mag   = .TRUE.
         lxc     = .TRUE.
         locc    = .TRUE.
         lbz     = .TRUE.
         lbs     = .TRUE.
         lsymm   = .TRUE.
         lefield = .TRUE.
         need_qexml = .TRUE.
         !
      CASE( 'all' )
         !
         lcell   = .TRUE.
         lpw     = .TRUE.
         lions   = .TRUE.
         lspin   = .TRUE.
         linit_mag  = .TRUE.
         lxc     = .TRUE.
         locc    = .TRUE.
         lbz     = .TRUE.
         lbs     = .TRUE.
         lwfc    = .TRUE.
         lsymm   = .TRUE.
         lefield = .TRUE.
         lrho    = .TRUE.
         need_qexml = .TRUE.
         !
      CASE( 'reset' )
         !
         lcell_read   = .FALSE.
         lpw_read     = .FALSE.
         lions_read   = .FALSE.
         lspin_read   = .FALSE.
         lstarting_mag_read   = .FALSE.
         lxc_read     = .FALSE.
         locc_read    = .FALSE.
         lbz_read     = .FALSE.
         lbs_read     = .FALSE.
         lwfc_read    = .FALSE.
         lsymm_read   = .FALSE.
         lefield_read = .FALSE.
         !
      CASE( 'ef' )
         !
         lef        = .TRUE.
         need_qexml = .TRUE.
         !
      CASE( 'exx' )
         !
         lexx       = .TRUE.
         need_qexml = .TRUE.
         !
      CASE( 'esm' )
         !
         lesm       = .TRUE.
         need_qexml = .TRUE.
         !
      END SELECT
      !
   
      electric_field_ispresent = input_obj%electric_field_ispresent 
      !
      IF ( lheader ) THEN 
         CALL readschema_header( gen_info )
      END IF 
      IF ( ldim ) THEN
         !         ! 

         ! 
         CALL readschema_dim(par_info, input_obj%atomic_species, output_obj%atomic_structure, output_obj%symmetries, &
                             output_obj%basis_set, output_obj%band_structure, input_obj) 
                                                                                                           
      ENDIF
      !
      IF ( lcell .AND. ( .NOT. lcell_read) ) THEN
         CALL readschema_cell( output_obj%atomic_structure,  input_obj )
      END IF
      !
      IF ( lpw ) THEN
         twfcollect = input_obj%control_variables%wf_collect
         CALL readschema_planewaves( output_obj%basis_set) 
      END IF
      IF ( lions ) THEN
         CALL readschema_ions( output_obj%atomic_structure, output_obj%atomic_species, input_obj, dirname)
      END IF
      IF ( lspin ) THEN

         CALL readschema_spin( output_obj%magnetization )
      END IF
      IF (linit_mag) THEN
         CALL readschema_magnetization (  output_obj%band_structure,  input_obj%atomic_species, input_obj)
      END IF
      IF ( lxc ) THEN
         CALL readschema_xc (  input_obj%atomic_species, output_obj%dft )
      END IF
      IF ( locc ) THEN
         CALL readschema_occupations( input_obj, output_obj%band_structure )
      END IF
      IF ( lbz ) THEN
         CALL readschema_brillouin_zone( input_obj%k_points_IBZ, input_obj%bands%occupations, &
                                         output_obj%symmetries,  output_obj%band_structure )
      END IF
      IF ( lbs ) THEN
         CALL readschema_band_structure( output_obj%band_structure )
      END IF
      IF ( lwfc ) THEN
         !
         IF (input_obj%control_variables%wf_collect) THEN 
            twfcollect =  input_obj%control_variables%wf_collect 
            CALL read_collected_to_evc(dirname, ierr ) 
         END IF
      END IF
      IF ( lsymm ) THEN
         CALL readschema_symmetry ( output_obj%symmetries, output_obj%basis_set, input_obj%symmetry_flags)
      END IF
      IF ( lefield ) THEN
         IF ( electric_field_ispresent ) THEN 
             CALL readschema_efield( input_obj%electric_field )
         ELSE 
             CALL readschema_efield()
         END IF
      END IF
      IF ( lrho ) THEN
         !
         ! ... to read the charge-density we use the routine from io_rho_xml 
         ! ... it also reads ns for ldaU and becsum for PAW
         !
         CALL read_rho( rho, nspin )
         !
      END IF
      IF ( lef ) THEN
               CALL readschema_ef ( output_obj%band_structure) 
         !
      END IF
      !
      IF ( lexx .AND. input_obj%dft%hybrid_ispresent  ) CALL readschema_exx ( input_obj%dft%hybrid )
      !
      IF ( lesm .AND. input_obj%boundary_conditions_ispresent ) THEN 
         IF ( input_obj%boundary_conditions%esm_ispresent ) &
                        CALL readschema_esm ( input_obj%boundary_conditions%esm) 
      END IF 
      !
      !
      !
      RETURN
      !
      ! uncomment to continue execution after an error occurs
      ! 100 IF (ionode .AND. need_qexml) THEN
      !        CALL qexml_closefile( 'read', IERR=tmp)
      !     ENDIF
      !     RETURN
      ! comment to continue execution after an error occurs
      !
    END SUBROUTINE init_vars_from_schema
    !-------------------------------------------------------------------------------
    SUBROUTINE readschema_header (gen_info_obj) 
    !-------------------------------------------------------------------------------
    ! 
    USE io_files,            ONLY: qexsd_fmt, qexsd_version, qexsd_init
    USE qes_types_module,    ONLY: general_info_type
    IMPLICIT NONE 
    !
    TYPE (general_info_type ),INTENT(IN)  :: gen_info_obj   
    ! 
    IF ( qexsd_init ) RETURN 
    qexsd_fmt = TRIM (gen_info_obj%xml_format%NAME)
    qexsd_version = TRIM ( gen_info_obj%xml_format%VERSION)
    qexsd_init = .TRUE. 
    END SUBROUTINE readschema_header 
    ! 
    !--------------------------------------------------------------------------
    SUBROUTINE readschema_dim(par_info_obj, atomic_species, atomic_structure, &
         symmetries, basis_set, band_structure, input_obj) 
      !
    USE constants,        ONLY : e2
    USE ions_base,        ONLY : nat, nsp
    USE symm_base,        ONLY : nsym
    USE gvect,            ONLY : ngm_g, ecutrho
    USE fft_base,         ONLY : dfftp
    USE gvecs,            ONLY : ngms_g, dual
    USE fft_base,         ONLY : dffts
    USE lsda_mod,         ONLY : lsda
    USE noncollin_module, ONLY : noncolin
    USE ktetra,           ONLY : ntetra
    USE klist,            ONLY : nkstot, nelec
    USE wvfct,            ONLY : nbnd, npwx
    USE gvecw,            ONLY : ecutwfc
    USE control_flags,    ONLY : gamma_only
    USE mp_pools,         ONLY : kunit
    USE mp_global,        ONLY : nproc_file, nproc_pool_file, &
                                 nproc_image_file, ntask_groups_file, &
                                 nproc_bgrp_file, nproc_ortho_file
    !
    USE qes_types_module, ONLY : parallel_info_type, atomic_species_type, atomic_structure_type, &
                                 symmetries_type, basis_set_type, band_structure_type, input_type  
    IMPLICIT NONE 
    !
    TYPE ( parallel_info_type ),INTENT(IN)     :: par_info_obj
    TYPE ( atomic_species_type ),INTENT(IN)    :: atomic_species
    TYPE ( atomic_structure_type ),INTENT(IN)  :: atomic_structure
    TYPE ( symmetries_type ),INTENT(IN)        :: symmetries
    TYPE ( basis_set_type ),INTENT(IN)         :: basis_set
    TYPE ( band_structure_type ),INTENT(IN)    :: band_structure 
    TYPE ( input_type ),INTENT(IN)             :: input_obj 
    ! 
    INTEGER                                    :: npwx_
    CALL readschema_cell ( atomic_structure, input_obj ) 
    ! 
    !---------------------------------------------------------------------
    !                                       PARALLEL  DIM 
    !----------------------------------------------------------------------
    nproc_file = par_info_obj%nprocs
    nproc_pool_file = nproc_file/par_info_obj%npool
    nproc_image_file = nproc_file 
    ntask_groups_file = par_info_obj%ntasks
    nproc_bgrp_file = nproc_image_file / par_info_obj%npool / par_info_obj%nbgrp 
    nproc_ortho_file = par_info_obj%ndiag
    !---------------------------------------------------------------------------
    !                                      ATOMS AND SPECIES 
    !--------------------------------------------------------------------------
    nsp = atomic_species%ntyp
    nat = atomic_structure%nat 
    !                                         SIMMETRIES 
    nsym = symmetries%nsym
    !-----------------------------------------------------------------------------
    !                                          BASIS SET 
    !-----------------------------------------------------------------------------
    ecutwfc = basis_set%ecutwfc*e2
    ecutrho = basis_set%ecutrho*e2
    dual = ecutrho/ecutwfc
    npwx_ = basis_set%npwx
    gamma_only= basis_set%gamma_only
    dfftp%nr1 = basis_set%fft_grid%nr1
    dfftp%nr2 = basis_set%fft_grid%nr2          
    dfftp%nr3 = basis_set%fft_grid%nr3
    dffts%nr1 = basis_set%fft_smooth%nr1
    dffts%nr2 = basis_set%fft_smooth%nr2
    dffts%nr3 = basis_set%fft_smooth%nr3
    ngm_g     = basis_set%ngm
    ngms_g    = basis_set%ngms
    !-------------------------------------------------------------------------
    !                                    BAND STRUCTURE  
    !-------------------------------------------------------------------------
    lsda  =    band_structure%lsda
    noncolin = band_structure%noncolin
    nelec =    band_structure%nelec
    nkstot =   band_structure%nks 
    IF ( lsda ) nkstot = nkstot * 2 
    nbnd = band_structure%nbnd
    IF ( input_obj%k_points_IBZ%monkhorst_pack_ispresent ) THEN 
       IF ( TRIM( input_obj%bands%occupations%occupations) == 'tetrahedra' ) THEN 
          ntetra = 6* input_obj%k_points_IBZ%monkhorst_pack%nk1* &
                      input_obj%k_points_IBZ%monkhorst_pack%nk2* &
                      input_obj%k_points_IBZ%monkhorst_pack%nk3 
       END IF 
    ELSE 
       ntetra = 0 
    END IF 
    ! 
    END SUBROUTINE readschema_dim
    !
    !-----------------------------------------------------------------------
    SUBROUTINE readschema_cell(atomic_structure, input_obj  )
    !-----------------------------------------------------------------------
    !
    USE constants,         ONLY : pi
    USE run_info,          ONLY : title
    USE cell_base,         ONLY : ibrav, alat, at, bg, celldm
    USE cell_base,         ONLY : tpiba, tpiba2, omega
    USE cellmd,            ONLY : lmovecell, cell_factor
    USE control_flags,     ONLY : do_makov_payne
    USE martyna_tuckerman, ONLY : do_comp_mt
    USE esm,               ONLY : do_comp_esm
    USE qes_types_module,  ONLY : input_type, atomic_structure_type
    !
    IMPLICIT NONE 
    ! 
    TYPE ( atomic_structure_type )            :: atomic_structure 
    TYPE ( input_type )                       :: input_obj
    !
    IF ( lcell_read ) RETURN 
    alat = atomic_structure%alat 
    IF ( atomic_structure%bravais_index_ispresent ) THEN 
       ibrav = atomic_structure%bravais_index 
    ELSE 
       ibrav = 0
    END IF
    celldm = 0.d0
    at(:,1) = atomic_structure%cell%a1
    at(:,2) = atomic_structure%cell%a2
    at(:,3) = atomic_structure%cell%a3
    SELECT CASE  (ibrav ) 
       CASE (1:3) 
          celldm(1) = alat
          celldm(2:6) = 0.d0
       CASE (4) 
          celldm(1) = alat
          celldm(2) = 0.d0
          celldm(3) = SQRT( DOT_PRODUCT(at(:,3),at(:,3)))/alat
          celldm(4:6) = 0.d0
       CASE (5) 
          celldm(1)= alat
          celldm(2:3) = 0.d0
          celldm(4) = DOT_PRODUCT(at(:,1),at(:,2))/(alat**2)
          celldm(5:6) = 0.d0
       CASE (6) 
          celldm(1)= alat 
          celldm(3)= SQRT( DOT_PRODUCT(at(:,2),at(:,2)))/alat
          celldm(2)= 1.d0
          celldm(4:6) = 0.d0
       CASE (7) 
          celldm(1) = alat
          celldm(3) = at(3,3) 
          celldm(2)=0.d0
          celldm(4:6) = 0.d0
       CASE (8)
          celldm(1) = alat
          celldm(2) = SQRT( DOT_PRODUCT (at(:,2),at(:,2)))/alat
          celldm(3) = SQRT( DOT_PRODUCT (at(:,3),at(:,3)))/alat 
          celldm(4:6) = 0.d0
       CASE (9) 
          celldm(1) = alat
          celldm(2) = ABS ( at(2,1)/at(1,1))
          celldm(3) = ABS ( at(3,3)/2.d0/at(1,1))
          celldm(4:6) = 0.d0 
       CASE (10) 
          celldm(1) = alat
          celldm(2) = ABS ( at(2,2)/at(2,1))
          celldm(3) = ABS ( at(3,1)/at(1,1))
          celldm(4:6) = 0.d0
       CASE (11) 
          celldm(1) = alat
          celldm(2) = ABS(at(2,1)/at(1,1))
          celldm(3) = ABS(at(3,1)/at(1,1))
          celldm(4:6) = 0.d0
       CASE (12) 
          celldm(1) = alat 
          celldm(2) = SQRT( DOT_PRODUCT(at(:,2),at(:,2))/DOT_PRODUCT(at(:,1),at(:,1)))
          celldm(3) = SQRT( DOT_PRODUCT(at(:,3),at(:,3))/DOT_PRODUCT(at(:,1),at(:,1)))
          celldm(4) = DOT_PRODUCT(at(:,1),at(:,2))/&
                      SQRT(DOT_PRODUCT(at(:,1),at(:,1))*DOT_PRODUCT(at(:,2),at(:,2)))
          celldm(5) =  DOT_PRODUCT(at(:,1),at(:,3))/&
                   SQRT(DOT_PRODUCT(at(:,1),at(:,1))*DOT_PRODUCT(at(:,3),at(:,3)))
          celldm(6) = 0.d0
       CASE (13) 
          celldm(1) = alat
          celldm(2) = SQRT( DOT_PRODUCT(at(:,2),at(:,2)))/(2.d0*at(1,1))
          celldm(3) = ABS (at(3,3)/at(1,3))
          celldm(4) = ATAN(at(2,2)/at(1,2))
          celldm(5:6) = 0.d0
       CASE (14) 
          celldm(1) = alat 
          celldm(2) = SQRT( DOT_PRODUCT(at(:,2),at(:,2))/DOT_PRODUCT(at(:,1),at(:,1)))
          celldm(3) = SQRT( DOT_PRODUCT(at(:,3),at(:,3))/DOT_PRODUCT(at(:,1),at(:,1)))
          celldm(4) = DOT_PRODUCT(at(:,3),at(:,2))/SQRT(DOT_PRODUCT(at(:,2),at(:,2))*&
                                                   DOT_PRODUCT(at(:,3),at(:,3)))
          celldm(5) = DOT_PRODUCT(at(:,3),at(:,1))/SQRT(DOT_PRODUCT(at(:,1),at(:,1))*&
                                                   DOT_PRODUCT(at(:,3),at(:,3)))
          celldm(6) = DOT_PRODUCT(at(:,1),at(:,2))/SQRT(DOT_PRODUCT(at(:,2),at(:,2))*&
                                                   DOT_PRODUCT(at(:,1),at(:,1)))
       CASE  default  
          celldm(1) = 1.d0
          IF (alat .GT. 0.d0 ) celldm(1) = alat
          celldm (2:6) = 0.d0
    END SELECT 
    tpiba = 2.d0*PI/alat
    tpiba2= tpiba**2
    omega = ABS (at(1,1)*at(2,2)*at(3,3)+at(1,2)*at(2,3)*at(3,1)+at(1,3)*at(2,1)*at(3,2)-&
                 at(3,1)*at(2,2)*at(1,3)-at(3,2)*at(2,3)*at(1,1)-at(3,3)*at(2,1)*at(1,2))
    at=at / alat
    CALL recips( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
    IF ( input_obj%boundary_conditions_ispresent )  THEN 
       SELECT CASE ( TRIM( input_obj%boundary_conditions%assume_isolated ))
         CASE ("makov-payne")
            do_makov_payne = .true.
            do_comp_mt     = .false.
            do_comp_esm    = .false. 
         CASE ("martyna-tuckerman")
            do_makov_payne = .false.
            do_comp_mt     = .true.
            do_comp_esm    = .false.
         CASE ("esm")
            do_makov_payne = .false.
            do_comp_mt     = .false.
            do_comp_esm    = .true.
         CASE ("none")
            do_makov_payne = .false.
            do_comp_mt     = .false.
            do_comp_esm    = .false.
       END SELECT 
    END IF         
    !
    title = TRIM(input_obj%control_variables%title)
    SELECT CASE ( TRIM ( input_obj%control_variables%calculation))
       CASE ('vc-relax', 'vc-md' ) 
           lmovecell = .TRUE. 
           IF ( input_obj%cell_control%cell_factor_ispresent ) THEN 
              cell_factor = input_obj%cell_control%cell_factor
           ELSE 
              cell_factor = 1.d0
           END IF 
       CASE default  
           lmovecell = .FALSE. 
    END SELECT 
    lcell_read = .TRUE.  
    !
    END SUBROUTINE readschema_cell
    ! 
    !------------------------------------------------------------------------
    SUBROUTINE readschema_ions( atomic_structure, atomic_species, input_obj, dirname ) 
    !------------------------------------------------------------------------
    ! 
    USE ions_base, ONLY : nat, nsp, ityp, amass, atm, tau, if_pos
    USE cell_base, ONLY : alat
    USE io_files,  ONLY : psfile, pseudo_dir, pseudo_dir_cur
    USE qes_types_module, ONLY: atomic_structure_type, atomic_species_type, input_type 
    ! 
    IMPLICIT NONE 
    ! 
    TYPE ( atomic_structure_type ),INTENT(IN) :: atomic_structure
    TYPE ( atomic_species_type ),INTENT(IN)        :: atomic_species 
    TYPE ( input_type ),INTENT(IN)            :: input_obj 
    CHARACTER(LEN=*), INTENT(IN)              :: dirname
    ! 
    INTEGER                                   :: iat, isp, idx
    CHARACTER(LEN = 3 ),ALLOCATABLE           :: symbols(:) 
    ! 
    IF ( lions_read ) RETURN 
    nat = atomic_structure%nat
    nsp = atomic_species%ntyp
    ALLOCATE ( symbols(nat) ) 
    DO isp = 1, nsp 
       amass(isp) = 0.d0 
       IF (atomic_species%species(isp)%mass_ispresent) amass(isp) = atomic_species%species(isp)%mass
       atm(isp) = TRIM ( atomic_species%species(isp)%name )
       psfile(isp) = TRIM ( atomic_species%species(isp)%pseudo_file) 
    END DO 
    !
    loop_on_atoms:DO iat = 1, nat
       idx = atomic_structure%atomic_positions%atom(iat)%index
       tau(:,idx) = atomic_structure%atomic_positions%atom(iat)%atom 
       symbols(idx)  = TRIM ( atomic_structure%atomic_positions%atom(idx)%name ) 
       loop_on_species:DO isp = 1, nsp
          IF ( TRIM(symbols(idx)) == TRIM (atm(isp))) THEN 
             ityp(iat) = isp 
             exit loop_on_species
          END IF 
       END  DO loop_on_species
    END DO loop_on_atoms
    !
    IF (ALLOCATED(if_pos) ) DEALLOCATE ( if_pos) 
    ALLOCATE (if_pos(3,nat) )
    IF ( input_obj%free_positions_ispresent ) THEN   
       if_pos = input_obj%free_positions%int_mat
    ELSE 
       if_pos = 1
    END IF 
    ! 
    IF ( .NOT. lcell_read .AND. ( atomic_structure%alat_ispresent )) alat = atomic_structure%alat 
    ! 
    pseudo_dir = TRIM(input_obj%control_variables%pseudo_dir)//'/'
    pseudo_dir_cur = TRIM ( dirname)//'/'  
    ! 
    lions_read = .TRUE.
    END SUBROUTINE readschema_ions
    !  
    !------------------------------------------------------------------------
    SUBROUTINE readschema_symmetry ( symms_obj, basis_obj  , symm_flags_obj) 
    !------------------------------------------------------------------------
      ! 
      USE symm_base,       ONLY : nrot, nsym, invsym, s, ft,ftau, irt, t_rev, &
                                 sname, sr, invs, inverse_s, s_axis_to_cart, &
                                 time_reversal, no_t_rev
      USE control_flags,   ONLY : noinv
      USE fft_base,        ONLY : dfftp
      USE qes_types_module,ONLY : symmetries_type, symmetry_type, basis_type
      ! 
      IMPLICIT NONE   
      ! 
      TYPE ( symmetries_type )               :: symms_obj 
      TYPE ( basis_set_type )                :: basis_obj
      TYPE ( symmetry_flags_type )           :: symm_flags_obj
      INTEGER                                :: isym 
      ! 
      IF (lsymm_read ) RETURN 
      nrot = symms_obj%nrot 
      nsym = symms_obj%nsym
      ! 
      noinv = symm_flags_obj%noinv 
      no_t_rev = symm_flags_obj%no_t_rev 
      !  
      invsym = .FALSE. 
      DO isym = 1, nrot
        s(:,:,isym) = symms_obj%symmetry(isym)%rotation%mat 
        sname(isym) = TRIM ( symms_obj%symmetry(isym)%info%name )  
        IF ( (TRIM(sname(isym)) == "inversion") .AND. (isym .LE. nsym) ) invsym = .TRUE.
        IF ( symms_obj%symmetry(isym)%fractional_translation_ispresent .AND. (isym .LE. nsym) ) THEN
           ft(1:3,isym)  =  symms_obj%symmetry(isym)%fractional_translation(1:3) 
           ftau(1,isym) = NINT( ft(1,isym)*DBLE( basis_obj%fft_grid%nr1 ) )
           ftau(2,isym) = NINT( ft(2,isym)*DBLE( basis_obj%fft_grid%nr2 ) )
           ftau(3,isym) = NINT( ft(3,isym)*DBLE( basis_obj%fft_grid%nr3 ) )
        END IF
        IF ( symms_obj%symmetry(isym)%info%time_reversal_ispresent ) THEN  
           IF (symms_obj%symmetry(isym)%info%time_reversal) THEN 
              t_rev( isym ) = 1
           ELSE
              t_rev( isym ) = 0 
           END IF 
        END IF
        IF ( symms_obj%symmetry(isym)%equivalent_atoms_ispresent .AND. (isym .LE. nsym) )   &
             irt(isym,:) = symms_obj%symmetry(isym)%equivalent_atoms%index_list(:)
      END DO
      CALL inverse_s()
      CALL s_axis_to_cart() 
      lsymm_read = .TRUE. 
    END SUBROUTINE readschema_symmetry 
    !
    !---------------------------------------------------------------------------
    SUBROUTINE readschema_efield( efield_obj  ) 
    !---------------------------------------------------------------------------
      !       
      USE extfield, ONLY : tefield, dipfield, edir, emaxpos, eopreg, eamp
      ! 
      IMPLICIT NONE 
      ! 
      TYPE ( electric_field_type),OPTIONAL, INTENT(IN)    :: efield_obj
      ! 
      !
      tefield = .FALSE. 
      dipfield = .FALSE. 
      IF ( .NOT. PRESENT( efield_obj) ) RETURN 
      IF (TRIM(efield_obj%electric_potential) == 'sawtooth_potential') THEN 
         tefield = .TRUE. 
         IF ( efield_obj%dipole_correction_ispresent ) THEN 
            dipfield = efield_obj%dipole_correction
         ELSE 
            dipfield = .FALSE. 
         END IF 
         IF ( efield_obj%electric_field_direction_ispresent ) THEN 
            edir = efield_obj%electric_field_direction
         ELSE 
            edir = 3 
         END IF
         IF ( efield_obj%potential_max_position_ispresent ) THEN 
            edir = efield_obj%electric_field_direction
         ELSE 
            emaxpos = 3 
         END IF 
         IF ( efield_obj%potential_decrease_width_ispresent ) THEN 
            eopreg = efield_obj%potential_decrease_width
         ELSE 
            eopreg = 1.d-1
         END IF 
         IF ( efield_obj%electric_field_amplitude_ispresent ) THEN 
            eamp = efield_obj%electric_field_amplitude
         ELSE 
            eamp = 1.d-3
         END IF
      END IF 
      !
  END SUBROUTINE readschema_efield  
    !-----------------------------------------------------------------------
    SUBROUTINE readschema_planewaves ( basis_set_obj ) 
    !-----------------------------------------------------------------------
    ! 
    USE constants,       ONLY : e2
    USE gvect,           ONLY : ngm_g, ecutrho
    USE gvecs,           ONLY : ngms_g, dual
    USE gvecw,           ONLY : ecutwfc
    USE fft_base,        ONLY : dfftp
    USE fft_base,        ONLY : dffts
    USE wvfct,           ONLY : npwx
    USE control_flags,   ONLY : gamma_only
    USE qes_types_module,ONLY : basis_set_type
    ! 
    IMPLICIT NONE 
    ! 
    TYPE ( basis_set_type )              :: basis_set_obj  
    !
    IF ( lpw_read ) RETURN 
    ecutwfc = basis_set_obj%ecutwfc*e2
    ecutrho = basis_set_obj%ecutrho*e2
    dual = ecutrho/ecutwfc
    !npwx = basis_set_obj%npwx
    gamma_only= basis_set_obj%gamma_only
    dfftp%nr1 = basis_set_obj%fft_grid%nr1
    dfftp%nr2 = basis_set_obj%fft_grid%nr2          
    dfftp%nr3 = basis_set_obj%fft_grid%nr3
    dffts%nr1 = basis_set_obj%fft_smooth%nr1
    dffts%nr2 = basis_set_obj%fft_smooth%nr2
    dffts%nr3 = basis_set_obj%fft_smooth%nr3
    ngm_g     = basis_set_obj%ngm
    ngms_g    = basis_set_obj%ngms
    !
    lpw_read = .TRUE.
    END SUBROUTINE readschema_planewaves 
    !--------------------------------------------------------------------------
    SUBROUTINE readschema_spin( magnetization_obj) 
    !--------------------------------------------------------------------------
      ! 
      USE spin_orb,         ONLY : lspinorb, domag
      USE lsda_mod,         ONLY : nspin, lsda
      USE noncollin_module, ONLY : noncolin, npol
      USE qes_types_module, ONLY : magnetization_type
      USE symm_base,        ONLY : time_reversal
      ! 
      IMPLICIT NONE 
      ! 
      TYPE ( magnetization_type ),INTENT(IN)         :: magnetization_obj
      ! 
      IF ( lspin_read ) RETURN
      lspinorb = magnetization_obj%spinorbit 
      domag =   magnetization_obj%do_magnetization 
      lsda  =   magnetization_obj%lsda
      noncolin = magnetization_obj%noncolin  
      IF ( noncolin .AND. domag ) time_reversal = .FALSE.
      IF ( lsda ) THEN  
        nspin = 2
        npol = 1
      ELSE IF (noncolin ) THEN 
        nspin = 4
        npol = 2
      ELSE 
        nspin =1
        npol = 1 
      END IF
      ! 
      lspin_read = .TRUE.
    END SUBROUTINE readschema_spin 
    !
    !------------------------------------------------------------------------
    SUBROUTINE readschema_magnetization( band_structure_obj, atomic_specs_obj, input_obj) 
      !----------------------------------------------------------------------
      ! 
      USE klist,            ONLY : two_fermi_energies, nelup, neldw
      USE ener,             ONLY : ef_up, ef_dw
      USE lsda_mod,         ONLY : starting_magnetization
      USE noncollin_module, ONLY : angle1, angle2, i_cons, mcons, bfield, &
                                   lambda
      USE electrons_base,   ONLY : set_nelup_neldw
      USE qes_types_module, ONLY : band_structure_type, atomic_species_type, input_type 
      !
      IMPLICIT NONE 
      ! 
      TYPE ( band_structure_type ),INTENT(IN)    :: band_structure_obj
      TYPE ( atomic_species_type ),INTENT(IN)    :: atomic_specs_obj
      TYPE ( input_type ), INTENT(IN)             :: input_obj
      !  
      REAL(DP)                   :: tot_mag_, nelec_, theta, phi, fixed_magnetization(3) 
      INTEGER                    :: nsp_, isp
      !
      IF ( lstarting_mag_read ) RETURN
      bfield = 0.d0
      nelec_ = band_structure_obj%nelec
      two_fermi_energies = band_structure_obj%two_fermi_energies_ispresent
      IF (two_fermi_energies) THEN 
         ef_up = band_structure_obj%two_fermi_energies(1)
         ef_dw = band_structure_obj%two_fermi_energies(2)
         IF ( input_obj%bands%tot_magnetization_ispresent )  THEN 
            tot_mag_ = input_obj%bands%tot_magnetization
         ELSE
            tot_mag_ = 0.d0
         END IF 
         CALL set_nelup_neldw( tot_mag_,  nelec_, nelup, neldw ) 
      END IF 
      nsp_ = atomic_specs_obj%ntyp
      !
      i_cons = 0
      IF (input_obj%spin_constraints_ispresent ) THEN
         lambda = input_obj%spin_constraints%lagrange_multiplier 
         SELECT CASE ( TRIM (input_obj%spin_constraints%spin_constraints ) )
            CASE ( 'atomic') 
               i_cons = 1 
            CASE ( 'atomic_direction' )
               i_cons =  2 
            CASE ( 'total' )
               i_cons = 3 
            CASE ( 'total_direction' ) 
               i_cons = 6 
         END SELECT  
         IF ( input_obj%spin_constraints%target_magnetization_ispresent ) THEN
            fixed_magnetization = input_obj%spin_constraints%target_magnetization 
            SELECT CASE ( i_cons) 
               CASE ( 3 ) 
                  mcons(1,1) = fixed_magnetization(1)
                  mcons(2,1) = fixed_magnetization(2)
                  mcons(3,1) = fixed_magnetization(3)
               CASE ( 6) 
                  mcons(3,1) = fixed_magnetization(3)
               CASE DEFAULT
                  CONTINUE
            END SELECT
         END IF
         !
      END IF
      DO isp = 1, nsp_
         IF ( atomic_specs_obj%species(isp)%starting_magnetization_ispresent) THEN
             starting_magnetization(isp) = atomic_specs_obj%species(isp)%starting_magnetization      
         END IF                                                                                      
         !                                                                                           
         IF ( band_structure_obj%noncolin ) THEN                                                         
            IF (    atomic_specs_obj%species(isp)%spin_teta_ispresent ) THEN 
               theta = atomic_specs_obj%species(isp)%spin_teta 
               angle1(isp) = theta 
            END IF                                                                  
            IF ( atomic_specs_obj%species(isp)%spin_phi_ispresent ) THEN                  
               phi = atomic_specs_obj%species(isp)%spin_phi
               angle2(isp) = phi
            END IF                                                                     
               !                                                                                     
            IF ( atomic_specs_obj%species(isp)%starting_magnetization_ispresent .AND. &
                                                                              i_cons == 1 ) THEN 
                !            
                mcons(1,isp) = starting_magnetization(isp) * sin( theta ) * cos( phi )
                mcons(2,isp) = starting_magnetization(isp) * sin( theta ) * sin( phi )
                mcons(3,isp) = starting_magnetization(isp) * cos( theta )
            ELSE IF ( i_cons == 2) THEN  
                mcons(3,isp) = cos(theta) 
            END IF
         ELSE IF ( atomic_specs_obj%species(isp)%starting_magnetization_ispresent .AND. &
                                                                                  i_cons == 1 ) THEN 
            mcons(1,isp) = starting_magnetization(isp)                                               
         END IF                                                                                      
      END DO   
      !
      lstarting_mag_read = .TRUE.
    END SUBROUTINE readschema_magnetization
    !-----------------------------------------------------------------------
    SUBROUTINE readschema_xc ( atomic_specs, dft_obj ) 
    !-----------------------------------------------------------------------
      ! 
      USE ions_base, ONLY : nsp
      USE funct,     ONLY : enforce_input_dft
      USE ldaU,      ONLY : lda_plus_u, lda_plus_u_kind, Hubbard_lmax, &
                            Hubbard_l, Hubbard_U, Hubbard_J, Hubbard_alpha, &
                            Hubbard_J0, Hubbard_beta, U_projection
      USE kernel_table,     ONLY : vdw_table_name
      USE control_flags,    ONLY : llondon, lxdm, ts_vdw
      USE london_module,    ONLY : scal6, lon_rcut
      USE tsvdw_module,     ONLY : vdw_isolated
      USE qes_types_module, ONLY : atomic_species_type, dft_type
      ! 
      IMPLICIT NONE 
      ! 
      TYPE ( atomic_species_type )    :: atomic_specs
      TYPE ( dft_type            )    :: dft_obj
      INTEGER                         :: ihub, nsp_, inlc, isp
      ! 
      CHARACTER(LEN = 256 )           :: label
      CHARACTER(LEN = 20  )           :: dft_name
      CHARACTER(LEN =  3 )            :: symbol
      ! 
      IF ( lxc_read ) RETURN
      dft_name = TRIM(dft_obj%functional) 
      CALL enforce_input_dft ( dft_name, .TRUE. ) 
      lda_plus_u = dft_obj%dftU_ispresent 
      IF ( lda_plus_u ) THEN 
         lda_plus_u_kind = dft_obj%dftU%lda_plus_u_kind
         U_projection = TRIM ( dft_obj%dftU%U_projection_type )
         Hubbard_l =-1 
         IF ( dft_obj%dftU%Hubbard_U_ispresent) THEN 
            loop_on_hubbardU:DO ihub =1, dft_obj%dftU%ndim_Hubbard_U
               symbol = TRIM(dft_obj%dftU%Hubbard_U(ihub)%specie)
               label  = TRIM(dft_obj%dftU%Hubbard_U(ihub)%label ) 
               loop_on_speciesU:DO isp = 1, nsp_
                  IF ( TRIM(symbol) == TRIM ( atomic_specs%species(isp)%name) ) THEN 
                     Hubbard_U(isp) = dft_obj%dftU%Hubbard_U(ihub)%HubbardCommon
                     SELECT CASE ( TRIM (label))
                        CASE ( '1s', '2s', '3s', '4s', '5s', '6s', '7s' ) 
                            Hubbard_l(isp) = 0 
                        CASE ( '2p', '3p', '4p', '5p', '6p' ) 
                            Hubbard_l(isp) = 1 
                        CASE ( '3d', '4d', '5d' ) 
                            Hubbard_l( isp ) = 2 
                        CASE ( '4f', '5f' )  
                            Hubbard_l(isp ) = 3
                        CASE  default 
                            CALL errore ("pw_readschema:", "unrecognized label for Hubbard "//label,&
                                          1 ) 
                     END SELECT   
                     EXIT loop_on_speciesU
                  END IF 
                END DO loop_on_speciesU
            END DO loop_on_hubbardU
         END IF 
         
         IF ( dft_obj%dftU%Hubbard_J0_ispresent ) THEN 
           loop_on_hubbardj0:DO ihub =1, dft_obj%dftU%ndim_Hubbard_J0
               symbol = TRIM(dft_obj%dftU%Hubbard_J0(ihub)%specie)
               loop_on_speciesj0:DO isp = 1, nsp_
                  IF ( TRIM(symbol) == TRIM ( atomic_specs%species(isp)%name) ) THEN
                     Hubbard_J0(isp) = dft_obj%dftU%Hubbard_J0(ihub)%HubbardCommon
                     EXIT loop_on_speciesj0
                  END IF
               END DO loop_on_speciesj0
            END DO loop_on_hubbardj0
         END IF
         IF ( dft_obj%dftU%Hubbard_alpha_ispresent) THEN 
            loop_on_hubbardAlpha:DO ihub =1, dft_obj%dftU%ndim_Hubbard_alpha
               symbol = TRIM(dft_obj%dftU%Hubbard_alpha(ihub)%specie)
               loop_on_speciesAlpha:DO isp = 1, nsp_
                  IF ( TRIM(symbol) == TRIM ( atomic_specs%species(isp)%name) ) THEN 
                     Hubbard_alpha(isp) = dft_obj%dftU%Hubbard_alpha(ihub)%HubbardCommon
                     EXIT loop_on_speciesAlpha
                  END IF 
                END DO loop_on_speciesAlpha
            END DO loop_on_hubbardAlpha
         END IF 
         IF ( dft_obj%dftU%Hubbard_beta_ispresent) THEN 
            loop_on_hubbardBeta:DO ihub =1, dft_obj%dftU%ndim_Hubbard_beta
               symbol = TRIM(dft_obj%dftU%Hubbard_beta(ihub)%specie)
               loop_on_speciesBeta:DO isp = 1, nsp_
                  IF ( TRIM(symbol) == TRIM ( atomic_specs%species(isp)%name) ) THEN 
                     Hubbard_beta(isp) = dft_obj%dftU%Hubbard_beta(ihub)%HubbardCommon
                     EXIT loop_on_speciesBeta
                  END IF 
                END DO loop_on_speciesBeta
            END DO loop_on_hubbardBeta
         END IF 
         IF ( dft_obj%dftU%Hubbard_J_ispresent) THEN 
            loop_on_hubbardJ:DO ihub =1, dft_obj%dftU%ndim_Hubbard_J
               symbol = TRIM(dft_obj%dftU%Hubbard_J(ihub)%specie)
               loop_on_speciesJ:DO isp = 1, nsp_
                  IF ( TRIM(symbol) == TRIM ( atomic_specs%species(isp)%name) ) THEN 
                     Hubbard_J(:,isp) = dft_obj%dftU%Hubbard_J(ihub)%HubbardJ
                     EXIT loop_on_speciesJ
                  END IF 
                END DO loop_on_speciesJ
            END DO loop_on_hubbardJ
         END IF 
         Hubbard_lmax = MAXVAL( Hubbard_l(1:nsp_) )
      END IF 
      ! 
      IF ( dft_obj%vdW_ispresent ) THEN 
           SELECT CASE( TRIM( dft_obj%vdW%vdw_corr ) )
    !
              CASE( 'grimme-d2', 'Grimme-D2', 'DFT-D', 'dft-d' )
                !
                llondon= .TRUE.
                ts_vdw= .FALSE.
                lxdm   = .FALSE.
                !
              CASE( 'TS', 'ts', 'ts-vdw', 'ts-vdW', 'tkatchenko-scheffler' )
                !
                llondon= .FALSE.
                ts_vdw= .TRUE.
                lxdm   = .FALSE.
                !
              CASE( 'XDM', 'xdm' )
                 !
                llondon= .FALSE.
                ts_vdw= .FALSE.
                lxdm   = .TRUE.
                !
              CASE DEFAULT
                !
                llondon= .FALSE.
                ts_vdw = .FALSE.
                lxdm   = .FALSE.
                !
           END SELECT
           SELECT CASE ( TRIM (dft_obj%vdW%non_local_term))
             CASE ('vdw1')  
                inlc = 1
             CASE ('vdw2') 
                inlc = 2
             CASE ('vv10' ) 
                inlc = 3 
             CASE ( 'vdW-DF-x') 
                inlc = 4
             CASE ( 'vdW-DF-y')
                inlc = 5
             CASE ( 'vdW-DF-z')
                inlc = 6
             CASE default 
                inlc = 0 
          END SELECT
          IF (inlc == 0 ) THEN 
             vdw_table_name = ' '
          ELSE IF ( inlc == 3 ) THEN 
             vdw_table_name = 'rVV10_kernel_table'
          ELSE
             vdw_table_name = 'vdW_kernel_table'
          END IF 
          IF (dft_obj%vdW%london_s6_ispresent ) THEN 
             scal6 = dft_obj%vdW%london_s6
          END IF 
          IF ( dft_obj%vdW%london_rcut_ispresent ) THEN 
             lon_rcut = dft_obj%vdW%london_rcut
          END IF 
          IF (dft_obj%vdW%ts_vdW_isolated_ispresent ) THEN 
             vdW_isolated = dft_obj%vdW%ts_vdW_isolated
          END IF 
      END IF 
      !         
      lxc_read = .TRUE.
    END SUBROUTINE readschema_xc
    !  
    !-----------------------------------------------------------------------------------------------------
    SUBROUTINE readschema_brillouin_zone( k_pointIBZ_obj , occupations_obj, symmetries_obj, band_struct_obj )
    !-----------------------------------------------------------------------------------------------------
       !
       USE lsda_mod, ONLY : lsda
       USE klist,    ONLY : nkstot, xk, wk, qnorm
       USE start_k,  ONLY : nks_start, xk_start, wk_start, &
                              nk1, nk2, nk3, k1, k2, k3
       USE symm_base,ONLY : nrot, s, sname
       USE qes_types_module, ONLY : k_points_IBZ_type, occupations_type, symmetries_type, band_structure_type
       !
       IMPLICIT NONE
       !
       TYPE ( k_points_IBZ_type ),  INTENT(IN)    :: k_pointIBZ_obj
       TYPE ( occupations_type ),   INTENT(IN)    :: occupations_obj
       TYPE ( symmetries_type ),    INTENT(IN)    :: symmetries_obj 
       TYPE ( band_structure_type ),INTENT(IN)    :: band_struct_obj 
       INTEGER                                    :: ik, isym, nks_
       ! 
       IF ( lbz_read ) RETURN 
       nks_ = band_struct_obj%nks
       nkstot = nks_
       IF ( band_struct_obj%lsda ) nkstot = nkstot * 2  
       ! 
       ! 
       DO ik = 1, nks_
          xk(:,ik) = band_struct_obj%ks_energies(ik)%k_point%k_point(:) 
       END DO 
       !   
       IF ( k_pointIBZ_obj%monkhorst_pack_ispresent ) THEN 
          nks_start = 0 
          nk1 = k_pointIBZ_obj%monkhorst_pack%nk1 
          nk2 = k_pointIBZ_obj%monkhorst_pack%nk2
          nk3 = k_pointIBZ_obj%monkhorst_pack%nk3 
           k1 = k_pointIBZ_obj%monkhorst_pack%k1
           k2 = k_pointIBZ_obj%monkhorst_pack%k2
           k3 = k_pointIBZ_obj%monkhorst_pack%k3
       ELSE IF (k_pointIBZ_obj%nk_ispresent .AND. &
                k_pointIBZ_obj%k_point_ispresent ) THEN 
           nks_start = k_pointIBZ_obj%nk
           ALLOCATE (xk_start(3,nks_start), wk_start(nks_start))
           DO ik =1, nks_start
               xk_start(:,ik) = k_pointIBZ_obj%k_point(ik)%k_point(:) 
               IF ( k_pointIBZ_obj%k_point(ik)%weight_ispresent) THEN 
                  wk_start(ik) = k_pointIBZ_obj%k_point(ik)%weight 
               ELSE 
                  wk_start(ik) = 0.d0
               END IF 
           END DO
       ELSE 
           CALL errore ("pw_readschema: ", &
                        " no information found for initializing brillouin zone information", 1)
       END IF  
       ! 
       IF ( .NOT. lsymm_read  ) THEN 
          nrot = symmetries_obj%nrot
          DO isym =1, symmetries_obj%ndim_symmetry
             s(:,:,isym)     = symmetries_obj%symmetry(isym)%rotation%mat
             sname(isym) = TRIM ( symmetries_obj%symmetry(isym)%info%name) 
         END DO 
       END IF    
       lbz_read = .TRUE.
    END SUBROUTINE readschema_brillouin_zone     
    !--------------------------------------------------------------------------------------------------
    SUBROUTINE readschema_occupations( input_obj, band_struct_obj ) 
      !------------------------------------------------------------------------------------------------
      ! 
      USE lsda_mod,         ONLY : lsda, nspin
      USE fixed_occ,        ONLY : tfixed_occ, f_inp
      USE ktetra,           ONLY : ntetra, ltetra
      USE klist,            ONLY : lgauss, ngauss, degauss, smearing
      USE electrons_base,   ONLY : nupdwn 
      USE wvfct,            ONLY : nbnd
      USE input_parameters, ONLY : input_parameters_occupations => occupations
      USE qes_types_module, ONLY : input_type, band_structure_type
      ! 
      IMPLICIT NONE 
      ! 
      TYPE ( input_type ),INTENT(IN)              :: input_obj
      TYPE ( band_structure_type ),INTENT(IN)     :: band_struct_obj 
      INTEGER                                     :: ispin, nk1, nk2, nk3, aux_dim1, aux_dim2 
      ! 
      IF ( locc_read ) RETURN 
      lsda= band_struct_obj%lsda
      nbnd = band_struct_obj%nbnd
      IF ( band_struct_obj%nbnd_up_ispresent ) nupdwn(1) = band_struct_obj%nbnd_up
      IF ( band_struct_obj%nbnd_dw_ispresent ) nupdwn(2) = band_struct_obj%nbnd_dw 
      IF ( lsda )  THEN 
         nspin = 2  
      ELSE IF ( band_struct_obj%noncolin) THEN 
         nspin = 4 
      ELSE 
         nspin = 1 
      END IF 
      !
      lgauss = .FALSE. 
      ltetra = .FALSE. 
      ngauss = 0
      input_parameters_occupations = TRIM ( input_obj%bands%occupations%occupations ) 
      IF (TRIM(input_obj%bands%occupations%occupations) == 'tetrahedra' ) THEN 
        ltetra = .TRUE. 
        nk1 = input_obj%k_points_IBZ%monkhorst_pack%nk1
        nk2 = input_obj%k_points_IBZ%monkhorst_pack%nk2
        nk3 = input_obj%k_points_IBZ%monkhorst_pack%nk3
        ntetra = 6* nk1 * nk2 * nk3 
      ELSE IF ( TRIM (input_obj%bands%occupations%occupations) == 'smearing') THEN 
        lgauss = .TRUE.  
        degauss = input_obj%bands%smearing%degauss
        SELECT CASE ( TRIM( input_obj%bands%smearing%smearing ) )
           CASE ( 'gaussian', 'gauss', 'Gaussian', 'Gauss' )
             ngauss = 0
             smearing  = 'gaussian'
           CASE ( 'methfessel-paxton', 'm-p', 'mp', 'Methfessel-Paxton', 'M-P', 'MP' )
             ngauss = 1
             smearing = 'Methfessel-Paxton'
           CASE ( 'marzari-vanderbilt', 'cold', 'm-v', 'mv', 'Marzari-Vanderbilt', 'M-V', 'MV')
             ngauss = -1
             smearing  = 'Marzari-Vanderbilt'
           CASE ( 'fermi-dirac', 'f-d', 'fd', 'Fermi-Dirac', 'F-D', 'FD')
             ngauss = -99
             smearing = 'Fermi-Dirac'
        END SELECT
      ELSE IF ( TRIM (input_obj%bands%occupations%occupations) == 'from_input' .AND. & 
                input_obj%bands%inputOccupations_ispresent ) THEN 
           tfixed_occ = .TRUE.
           IF ( .NOT. ALLOCATED(f_inp))  &
              aux_dim2 = input_obj%bands%ndim_inputOccupations
              aux_dim1 = MAXVAL(input_obj%bands%inputOccupations(1:aux_dim2)%ndim_vec)
              ALLOCATE (f_inp( aux_dim1, aux_dim2))
           DO ispin = 1, input_obj%bands%ndim_inputOccupations
              f_inp(:,ispin) = input_obj%bands%inputOccupations(ispin)%vec
           END DO
      END IF       
     !
     locc_read = .TRUE.
    END SUBROUTINE readschema_occupations
 !
    !------------------------------------------------------------------------
    SUBROUTINE readschema_band_structure( band_struct_obj )
      !------------------------------------------------------------------------
      !
      USE control_flags, ONLY : lkpoint_dir
      USE constants,     ONLY : e2
      USE basis,    ONLY : natomwfc
      USE lsda_mod, ONLY : lsda, isk
      USE klist,    ONLY : nkstot, wk, nelec
      USE wvfct,    ONLY : et, wg, nbnd
      USE ener,     ONLY : ef, ef_up, ef_dw
      USE qes_types_module, ONLY : band_structure_type
      !
      IMPLICIT NONE
      TYPE ( band_structure_type)         :: band_struct_obj
      INTEGER                             :: ik, nbnd_, nbnd_up_, nbnd_dw_
      ! 
      lkpoint_dir = .FALSE.
      lsda = band_struct_obj%lsda
      nkstot = band_struct_obj%nks 
      IF ( lsda) THEN 
        nkstot = nkstot * 2 
        isk(1:nkstot/2) = 1
        isk(nkstot/2+1:nkstot) = 2 
      ELSE 
        isk(1:nkstot)   = 1 
      END IF 
      ! 
      nelec = band_struct_obj%nelec
      nbnd  = band_struct_obj%nbnd 
      natomwfc = band_struct_obj%num_of_atomic_wfc
      IF ( band_struct_obj%fermi_energy_ispresent) THEN 
         ef = band_struct_obj%fermi_energy*e2 
      ELSE IF ( band_struct_obj%two_fermi_energies_ispresent ) THEN 
         ef = 0.d0 
         ef_up = band_struct_obj%two_fermi_energies(1)*e2
         ef_dw = band_struct_obj%two_fermi_energies(2)*e2
      ELSE 
         ef = 0.d0
         ef_up = 0.d0
         ef_dw = 0.d0
      END IF 
      DO ik =1, band_struct_obj%ndim_ks_energies
         IF ( band_struct_obj%lsda) THEN
            IF ( band_struct_obj%nbnd_up_ispresent .AND. band_struct_obj%nbnd_dw_ispresent) THEN
               nbnd_up_ = band_struct_obj%nbnd_up
               nbnd_dw_ = band_struct_obj%nbnd_dw 
            ELSE IF ( band_struct_obj%nbnd_up_ispresent ) THEN 
               nbnd_up_ = band_struct_obj%nbnd_up
               nbnd_dw_ = band_struct_obj%ks_energies(ik)%ndim_eigenvalues - nbnd_up_
            ELSE IF ( band_struct_obj%nbnd_dw_ispresent ) THEN 
               nbnd_dw_ = band_struct_obj%nbnd_dw
               nbnd_up_ = band_struct_obj%ks_energies(ik)%ndim_eigenvalues - nbnd_dw_ 
            ELSE 
               nbnd_up_ = band_struct_obj%ks_energies(ik)%ndim_eigenvalues/2  
               nbnd_dw_ = band_struct_obj%ks_energies(ik)%ndim_eigenvalues/2
            END IF
            wk(ik) = band_struct_obj%ks_energies(ik)%k_point%weight
            wk( ik + band_struct_obj%ndim_ks_energies ) = wk(ik) 
            et(1:nbnd_up_,ik) = band_struct_obj%ks_energies(ik)%eigenvalues(1:nbnd_up_)*e2
            et(1:nbnd_dw_,ik+band_struct_obj%ndim_ks_energies) =  &
                             band_struct_obj%ks_energies(ik)%eigenvalues(nbnd_up_+1:nbnd_up_+nbnd_dw_)*e2
            wg(1:nbnd_up_,ik) = band_struct_obj%ks_energies(ik)%occupations(1:nbnd_up_)*wk(ik)
            wg(1:nbnd_dw_,ik+band_struct_obj%ndim_ks_energies) =  &
                             band_struct_obj%ks_energies(ik)%occupations(nbnd_up_+1:nbnd_up_+nbnd_dw_)*wk(ik)
         ELSE 
            wk(ik) = band_struct_obj%ks_energies(ik)%k_point%weight
            nbnd_ = band_struct_obj%ks_energies(ik)%ndim_eigenvalues
            et (1:nbnd_,ik) = band_struct_obj%ks_energies(ik)%eigenvalues(1:nbnd_)*e2
            wg (1:nbnd_,ik) = band_struct_obj%ks_energies(ik)%occupations(1:nbnd_)*wk(ik)
         END IF  
      END DO 
    END SUBROUTINE readschema_band_structure 
    ! 
    !------------------------------------------------------------------------
    SUBROUTINE read_collected_to_evc( dirname, ierr )
      !------------------------------------------------------------------------
      !
      ! ... This routines reads wavefunctions from the new file format and
      ! ... writes them into the old format
      !
      USE control_flags,        ONLY : twfcollect, lkpoint_dir
      USE cell_base,            ONLY : tpiba2
      USE lsda_mod,             ONLY : nspin, isk
      USE klist,                ONLY : nkstot, wk, nks, xk, ngk
      USE wvfct,                ONLY : npw, npwx, g2kin, et, wg, nbnd
      USE gvecw,                ONLY : ecutwfc
      USE wavefunctions_module, ONLY : evc
      USE io_files,             ONLY : nwordwfc, iunwfc
      USE buffers,              ONLY : save_buffer
      USE gvect,                ONLY : ngm, ngm_g, g, ig_l2g
      USE noncollin_module,     ONLY : noncolin, npol
      USE mp_images,            ONLY : nproc_image, intra_image_comm
      USE mp_pools,             ONLY : kunit, nproc_pool, me_pool, root_pool, &
                                       intra_pool_comm, inter_pool_comm
      USE mp_bands,             ONLY : me_bgrp, nbgrp, root_bgrp, &
                                       intra_bgrp_comm
      USE mp,                   ONLY : mp_bcast, mp_sum, mp_max
      USE xml_io_base,          ONLY : read_wfc
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      CHARACTER(LEN=320)    :: filename
      INTEGER              :: ik, ig, ipol, ik_eff, num_k_points
      INTEGER, ALLOCATABLE :: kisort(:)
      INTEGER              :: nkl, nkr, npwx_g
      INTEGER              :: nupdwn(2), ike, iks, npw_g, ispin
      INTEGER, ALLOCATABLE :: ngk_g(:)
      INTEGER, ALLOCATABLE :: igk_l2g(:,:), igk_l2g_kdip(:,:)
      LOGICAL              :: opnd
      REAL(DP)             :: scalef
      REAL(DP),ALLOCATABLE :: gkin_aux(:)

      !
      ! The ierr output var is actually not given any value
      ! except this initialization
      !
      ierr = 0
      !
      IF ( iunwfc > 0 ) THEN
         !
         INQUIRE( UNIT = iunwfc, OPENED = opnd )
         !
         IF ( .NOT. opnd ) CALL errore( 'read_wavefunctions', &
                    & 'wavefunctions unit (iunwfc) is not opened', 1 )
      END IF
      ! 
      CALL  kpoint_global_indices (nkstot, iks, ike)
      !
      ! ... find out the global number of G vectors: ngm_g  
      !
      ngm_g = ngm
      !
      CALL mp_sum( ngm_g, intra_bgrp_comm )
      !
      ! ... build the igk_l2g array, yielding the correspondence between
      ! ... the local k+G index and the global G index - see also ig_l2g
      !
      ALLOCATE ( igk_l2g( npwx, nks ) )
      igk_l2g = 0
      !
      ! FIXME: is it really needed to re-generate and discard k+G indices here?
      !
      ALLOCATE( gkin_aux(ngm), kisort( npwx ) )
      DO ik = 1, nks
         kisort = 0
         npw    = npwx
         !
         CALL gk_sort( xk(1,ik+iks-1), ngm, g, &
                       ecutwfc/tpiba2, npw, kisort(1), gkin_aux )
         DO ig = 1, npw
            igk_l2g(ig,ik) = ig_l2g(kisort(ig))
         END DO
         !
      END DO
      DEALLOCATE ( kisort, gkin_aux )
      !
      ! ... compute the global number of G+k vectors for each k point
      !
      ALLOCATE( ngk_g( nkstot ) )
      !
      ngk_g = 0
      ngk_g(iks:ike) = ngk(1:nks)
      !
      CALL mp_sum( ngk_g, inter_pool_comm )
      CALL mp_sum( ngk_g, intra_pool_comm )
      ngk_g = ngk_g / nbgrp
      !
      ! ... compute the Maximum G vector index among all G+k an processors
      !
      npw_g = MAXVAL( igk_l2g(:,:) )
      !
      CALL mp_max( npw_g, inter_pool_comm )
      CALL mp_max( npw_g, intra_pool_comm )

      !
      ! ... compute the Maximum number of G vector among all k points
      !
      npwx_g = MAXVAL( ngk_g(1:nkstot) )
      !
      ! 
      ! ... define a further l2g map to read gkvectors and wfc coherently 
      ! 
      ALLOCATE( igk_l2g_kdip( npwx_g, nks ) )
      igk_l2g_kdip = 0
      !
      DO ik = iks, ike
         !
         CALL gk_l2gmap_kdip( npw_g, ngk_g(ik), ngk(ik-iks+1), &
                              igk_l2g(1,ik-iks+1), igk_l2g_kdip(1,ik-iks+1) )
      END DO
      !
      num_k_points = nkstot
      !
      IF ( nspin == 2 ) num_k_points = nkstot / 2
      !
      IF ( .NOT. twfcollect ) RETURN 
      !
      k_points_loop: DO ik = 1, num_k_points
         !
         IF ( nspin == 2 ) THEN
            !
            ispin = 1 
            evc=(0.0_DP, 0.0_DP)
            !
            ! ... no need to read isk here: they are read from band structure
            ! ... and correctly distributed across pools in read_file
            !!! isk(ik) = 1
            !
            IF ( ionode ) filename = TRIM(dirname)//'/wfcup'//TRIM(int_to_char(ik))//'.dat'
            !
            CALL read_wfc( iunpun, ik, nkstot, kunit, ispin, nspin,      &
                           evc, npw_g, nbnd, igk_l2g_kdip(:,ik-iks+1),   &
                           ngk(ik-iks+1), filename, scalef, &
                           ionode, root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm )
            !
            IF ( ( ik >= iks ) .AND. ( ik <= ike ) ) THEN
               !
               CALL save_buffer ( evc, nwordwfc, iunwfc, (ik-iks+1) )
               !
            END IF
            !
            ispin = 2
            ik_eff = ik + num_k_points
            evc=(0.0_DP, 0.0_DP)
            !
            ! ... no need to read isk here (see above why)
            !isk(ik_eff) = 2
            !
            IF ( ionode ) filename = TRIM(dirname)//'/wfcdw'//TRIM(int_to_char(ik))//'.dat'
            !
            CALL read_wfc( iunpun, ik_eff, nkstot, kunit, ispin, nspin,      &
                           evc, npw_g, nbnd, igk_l2g_kdip(:,ik_eff-iks+1),   &
                           ngk(ik_eff-iks+1), filename, scalef, &
                           ionode, root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm )
            !
            IF ( ( ik_eff >= iks ) .AND. ( ik_eff <= ike ) ) THEN
               !
               CALL save_buffer ( evc, nwordwfc, iunwfc, (ik_eff-iks+1) )
               !
            END IF
            !
         ELSE
            !
            ! ... no need to read isk here (see above why)
            !isk(ik) = 1
            !
            evc=(0.0_DP, 0.0_DP)
            IF ( noncolin ) THEN
               !
               DO ipol = 1, npol
                  !
                  IF ( ipol == 1 ) THEN
                     filename = TRIM(dirname)//'wfcup'//TRIM(int_to_char(ik))//'.dat'
                  ELSE
                     filename = TRIM(dirname)//'wfcdw'//TRIM(int_to_char(ik))//'.dat'
                  END IF
                  !
                  !!! TEMP
                  nkl=(ipol-1)*npwx+1
                  nkr= ipol   *npwx
                  CALL read_wfc( iunpun, ik, nkstot, kunit, ispin,          &
                                 npol, evc(nkl:nkr,:), npw_g, nbnd,         &
                                 igk_l2g_kdip(:,ik-iks+1), ngk(ik-iks+1),   &
                                 filename, scalef, & 
                                 ionode, root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm )
                  !
               END DO
               !
            ELSE
               !
               IF ( ionode ) filename = TRIM(dirname)//'/wfc'//TRIM(int_to_char(ik))//'.dat'
               !
               CALL read_wfc( iunpun, ik, nkstot, kunit, ispin, nspin,         &
                              evc, npw_g, nbnd, igk_l2g_kdip(:,ik-iks+1),      &
                              ngk(ik-iks+1), filename, scalef, &
                              ionode, root_pool, intra_pool_comm, inter_pool_comm, intra_image_comm )
               !
            END IF
            !
            IF ( ( ik >= iks ) .AND. ( ik <= ike ) ) THEN
               CALL save_buffer ( evc, nwordwfc, iunwfc, (ik-iks+1) )
            END IF
            !
         END IF
         !
      END DO k_points_loop
      !
      DEALLOCATE ( igk_l2g )
      DEALLOCATE ( igk_l2g_kdip )
      !
      RETURN
      !

    END SUBROUTINE read_collected_to_evc
    !
    !----------------------------------------------------------------------------------------
    SUBROUTINE readschema_ef ( band_struct_obj )
    !----------------------------------------------------------------------------------------
       !
       USE ener,  ONLY            : ef, ef_up, ef_dw
       USE klist, ONLY            : two_fermi_energies, nelec
       USE qes_types_module, ONLY : band_structure_type 
       ! 
       IMPLICIT NONE 
       ! 
       TYPE ( band_structure_type ),INTENT(IN)      :: band_struct_obj 
       ! 
       two_fermi_energies = band_struct_obj%two_fermi_energies_ispresent 
       nelec = band_struct_obj%nelec
       IF ( two_fermi_energies) THEN 
          ef_up = band_struct_obj%two_fermi_energies(1) 
          ef_dw = band_struct_obj%two_fermi_energies(2)
       ELSE IF ( band_struct_obj%fermi_energy_ispresent ) THEN 
          ef = band_struct_obj%fermi_energy
       END IF 
    END SUBROUTINE readschema_ef 
    !------------------------------------------------------------------------
    SUBROUTINE readschema_exx ( hybrid_obj) 
    !------------------------------------------------------------------------
      ! 
      USE funct,                ONLY : set_exx_fraction, set_screening_parameter, &
                                      set_gau_parameter, enforce_input_dft, start_exx
      USE exx,                  ONLY : x_gamma_extrapolation, nq1, nq2, nq3, &
                                       exxdiv_treatment, yukawa, ecutvcut, ecutfock
      ! 
      USE  qes_types_module,   ONLY : hybrid_type 
      IMPLICIT NONE
      ! 
      TYPE ( hybrid_type), INTENT(IN)               :: hybrid_obj 
      ! 
      x_gamma_extrapolation = hybrid_obj%x_gamma_extrapolation
      nq1 = hybrid_obj%qpoint_grid%nqx1
      nq2 = hybrid_obj%qpoint_grid%nqx2
      nq3 = hybrid_obj%qpoint_grid%nqx3
      CALL set_exx_fraction( hybrid_obj%exx_fraction) 
      CALL set_screening_parameter ( hybrid_obj%screening_parameter) 
      ecutvcut = hybrid_obj%ecutvcut
      ecutfock = hybrid_obj%ecutfock
      CALL start_exx() 
    END SUBROUTINE  readschema_exx 
    !-----------------------------------------------------------------------------------
    SUBROUTINE readschema_esm ( esm_obj ) 
    !-----------------------------------------------------------------------------------
       ! 
       USE esm, ONLY : esm_nfit, esm_efield, esm_w, esm_a, esm_bc  
       USE qes_types_module, ONLY : esm_type
       ! 
       IMPLICIT NONE 
       ! 
       TYPE ( esm_type ), INTENT(IN)    :: esm_obj 
       ! 
       esm_nfit =   esm_obj%nfit
       esm_efield = esm_obj%efield
       esm_w      = esm_obj%w 
       esm_bc     = esm_obj%bc 
       esm_a      = 0.d0 
    END SUBROUTINE readschema_esm 
    !   
#endif
  END MODULE pw_restart_new
