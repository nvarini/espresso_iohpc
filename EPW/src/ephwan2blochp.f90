  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! 
  !---------------------------------------------------------------------------
  subroutine ephwan2blochp ( nmodes, xxq, irvec, ndegen, nrr_q, cuf, epmatf, nbnd, nrr_k )
  !---------------------------------------------------------------------------
  !!
  !! even though this is for phonons, I use the same notations
  !! adopted for the electronic case (nmodes->nmodes etc)
  !!
  USE kinds,         only : DP
  USE epwcom,        only : parallel_k, parallel_q, etf_mem
  USE elph2,         only : epmatwp
  USE constants_epw, ONLY : twopi, ci, czero
  USE io_files,      ONLY : prefix, tmp_dir
  USE io_epw,        ONLY : iunepmatwp
  USE mp_global,     ONLY : mp_sum
  USE mp_world,      ONLY : world_comm
  USE parallel_include
  implicit none
  !
  !  input variables
  !
  INTEGER, INTENT (in) :: nmodes
  !! Total number of modes
  INTEGER, INTENT (in) :: nrr_q
  !! Number of WS points
  INTEGER, INTENT (in) :: irvec ( 3, nrr_q)
  !! Coordinates of WS points
  INTEGER, INTENT (in) :: ndegen (nrr_q)
  !! Number of degeneracy of WS points
  INTEGER, INTENT (in) :: nbnd
  !! Number of bands
  INTEGER, INTENT (in) ::  nrr_k
  !! Number of electronic WS points
  REAL(kind=DP) :: xxq(3)
  !! Kpoint for the interpolation (WARNING: this must be in crystal coord!)
  COMPLEX(kind=DP), INTENT (in) :: cuf (nmodes, nmodes)
  !! e-p matrix in Wanner representation
  COMPLEX(kind=DP), INTENT (out) :: epmatf (nbnd, nbnd, nrr_k, nmodes)
  !! e-p matrix in Bloch representation, fine grid
  ! 
  ! Local variables 
  !
  character (len=256) :: filint
  integer :: ir, ir_start, ir_stop, iunepmatwp2, ierr
  integer(kind=8) ::  lrepmatw,  lrepmatw2
  real(kind=DP) :: rdotk
  complex(kind=DP) :: eptmp( nbnd, nbnd, nrr_k, nmodes)
  complex(kind=DP) :: cfac(nrr_q)
  !complex(kind=DP):: aux( nbnd*nbnd*nrr_k*nmodes )
  complex(kind=DP), allocatable :: epmatw ( :,:,:,:)
  !
  CALL start_clock('ephW2Bp')
  !----------------------------------------------------------
  !  STEP 3: inverse Fourier transform of g to fine k mesh
  !----------------------------------------------------------
  !
  !  g~ (k') = sum_R 1/ndegen(R) e^{-ik'R} g (R)
  !
  !  g~(k') is epmatf (nmodes, nmodes, ik )
  !  every pool works with its own subset of k points on the fine grid
  !
  IF (parallel_k) THEN
     CALL para_bounds(ir_start, ir_stop, nrr_q)
  ELSEIF (parallel_q) THEN
     ir_start = 1
     ir_stop  = nrr_q
  ELSE 
     CALL errore ('ephwan2blochp', 'Problem with parallel_k/q scheme', nrr_q)
  ENDIF
  !
#ifdef __MPI
  IF (.NOT. etf_mem) then
    ! Check for directory given by "outdir"
    !      
    filint = trim(tmp_dir)//trim(prefix)//'.epmatwp1'
    CALL MPI_FILE_OPEN(world_comm,filint,MPI_MODE_RDONLY,MPI_INFO_NULL,iunepmatwp2,ierr)
    IF( ierr /= 0 ) CALL errore( 'ephwan2blochp', 'error in MPI_FILE_OPEN',1 )
    IF( parallel_q ) CALL errore( 'ephwan2blochp', 'q-parallel+etf_mem=.false. is not supported',1 ) 
    !CALL MPI_COMM_RANK(world_comm,my_id,ierr)
  ENDIF
#endif  
  !
  eptmp = czero
  cfac(:) = czero
  !
  DO ir = ir_start, ir_stop
     !   
     ! note xxq is assumed to be already in cryst coord
     !
     rdotk = twopi * dot_product ( xxq, dble(irvec(:, ir)) )
     cfac(ir) = exp( ci*rdotk ) / dble( ndegen(ir) )
  ENDDO
  ! 
  IF (etf_mem) then
    !DO ir = ir_start, ir_stop
    !  eptmp(:,:,:,:) = eptmp(:,:,:,:) +&
    !    cfac(ir)*epmatwp( :, :, :, :, ir)
    !ENDDO
    ! SP: This is faster by 20 % 
    Call zgemv( 'n',  nbnd * nbnd * nrr_k * nmodes, ir_stop - ir_start + 1, ( 1.d0, 0.d0 ),&
             epmatwp(1,1,1,1,ir_start), nbnd * nbnd * nrr_k * nmodes, cfac(ir_start),1,( 0.d0, 0.d0),eptmp, 1 )    
    !
  ELSE
    !
    ALLOCATE(epmatw ( nbnd, nbnd, nrr_k, nmodes))
    !
    lrepmatw2   = 2 * nbnd * nbnd * nrr_k * nmodes
    ! 
    DO ir = ir_start, ir_stop
#ifdef __MPI
      ! DEBUG: print*,'Process ',my_id,' do ',ir,'/ ',ir_stop
      !
      !  Direct read of epmatwp for this ir
      lrepmatw   = 2 * nbnd * nbnd * nrr_k * nmodes * 8 * (ir-1)
      ! SP: mpi seek is used to set the position at which we should start
      ! reading the file. It is given in bits. 
      ! Note : The process can be collective (=blocking) if using MPI_FILE_SET_VIEW & MPI_FILE_READ_ALL
      !        or noncollective (=non blocking) if using MPI_FILE_SEEK & MPI_FILE_READ. 
      !        Here we want non blocking because not all the process have the same nb of ir. 
      !
      CALL MPI_FILE_SEEK(iunepmatwp2,lrepmatw,MPI_SEEK_SET,ierr)
      IF( ierr /= 0 ) CALL errore( 'ephwan2blochp', 'error in MPI_FILE_SEEK',1 )
      !CALL MPI_FILE_READ(iunepmatwp2, aux, lrepmatw2, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierr)
      CALL MPI_FILE_READ(iunepmatwp2, epmatw, lrepmatw2, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierr)
      IF( ierr /= 0 ) CALL errore( 'ephwan2blochp', 'error in MPI_FILE_READ_ALL',1 )
      ! 
     ! i = 0
     ! DO imode = 1, nmodes
     !  DO ip = 1, nrr_k
     !   DO jbnd = 1, nbnd
     !    DO ibnd = 1, nbnd
     !      i = i + 1
     !      epmatw ( ibnd, jbnd, ip, imode ) = aux (i)
     !      ! 
     !    ENDDO
     !   ENDDO
     !  ENDDO
     ! ENDDO
#else      
      call rwepmatw ( epmatw, nbnd, nrr_k, nmodes, ir, iunepmatwp, -1)
#endif
      !
      !eptmp = eptmp + cfac(ir)*epmatw
      CALL ZAXPY(nbnd * nbnd * nrr_k * nmodes, cfac(ir), epmatw, 1, eptmp, 1)
      ! 
    ENDDO
    DEALLOCATE(epmatw)
  ENDIF
  !
#ifdef __MPI
  IF (parallel_k) CALL mp_sum(eptmp, world_comm)
  IF (.NOT. etf_mem) then
    CALL MPI_FILE_CLOSE(iunepmatwp2,ierr)
    IF( ierr /= 0 ) CALL errore( 'ephwan2blochp', 'error in MPI_FILE_CLOSE',1 )
  ENDIF  
#endif  
  !
  !----------------------------------------------------------
  !  STEP 4: un-rotate to Bloch space, fine grid
  !----------------------------------------------------------
  !
  ! epmatf(j) = sum_i eptmp(i) * uf(i,j)
  !
  Call zgemm( 'n', 'n', nbnd * nbnd * nrr_k, nmodes, nmodes, ( 1.d0, 0.d0 ),eptmp  , nbnd * nbnd * nrr_k, &
                                                                                    cuf, nmodes         , &
                                                             ( 0.d0, 0.d0 ),epmatf, nbnd * nbnd * nrr_k )
 ! DO ibnd = 1, nbnd
 !  DO jbnd = 1, nbnd
 !    !
 !    CALL zgemm ('n','n',nrr_k, nmodes, nmodes,(1.d0,0.d0),eptmp(ibnd,jbnd,:,:),nrr_k,&
 !        cuf,nmodes,(0.d0,0.d0), epmatf(ibnd,jbnd,:,:), nrr_k )
 !    !
 !  ENDDO
 ! ENDDO
  !
  CALL stop_clock('ephW2Bp')
  !
  end subroutine ephwan2blochp

