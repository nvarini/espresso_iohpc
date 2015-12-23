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
MODULE distribution_center
  !-----------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  SAVE
  !
  TYPE :: index_distribution
     INTEGER,ALLOCATABLE :: l2g(:)
     INTEGER :: nloc
     INTEGER :: nlocx
     INTEGER :: nglob
     INTEGER :: myoffset
     INTEGER :: mpicomm
     CHARACTER(LEN=1) :: info
  END TYPE index_distribution
  !
  !
  TYPE(index_distribution) :: pert
  TYPE(index_distribution) :: macropert
  TYPE(index_distribution) :: ifr
  TYPE(index_distribution) :: rfr
  TYPE(index_distribution) :: aband
  !
  CONTAINS
  !
  !
  !
  SUBROUTINE index_distr_init(dato,n,level,label,lverbose)
    !
    USE mp_global,   ONLY : &
                          & my_image_id,nimage,inter_image_comm, &
                          & my_pool_id,npool,inter_pool_comm, &
                          & my_bgrp_id,nbgrp,inter_bgrp_comm, &
                          & me_bgrp,nproc_bgrp,intra_bgrp_comm
    USE mp_world,    ONLY : mpime,root,world_comm,nproc
    USE io_global,   ONLY : stdout
    USE mp,          ONLY : mp_max,mp_min,mp_sum
    USE io_push,     ONLY : io_push_title,io_push_value,io_push_bar
    !
    IMPLICIT NONE
    !
    ! I/O
    !
    TYPE(index_distribution) :: dato 
    INTEGER,INTENT(IN) :: n
    CHARACTER(LEN=1),INTENT(IN) :: level
    CHARACTER(LEN=*),INTENT(IN) :: label
    LOGICAL,INTENT(IN) :: lverbose
    !
    ! Workspace
    ! 
    INTEGER :: mylevelid, nlevel, level_comm
    INTEGER :: nloc, nlocx, nglob, nloc_min, myoffset
    INTEGER :: iloc, iglob, ii
    INTEGER,ALLOCATABLE :: tmp_tab(:)
    CHARACTER(LEN=1) :: info
    !
    !
    SELECT CASE( level ) 
    CASE('W','w')
       !
       mylevelid  = mpime
       nlevel     = nproc
       level_comm = world_comm
       info       = 'w'
       !
    CASE('I','i')
       !
       mylevelid  = my_image_id
       nlevel     = nimage
       level_comm = inter_image_comm
       info       = 'i'
       !
    CASE('P','p')
       !
       mylevelid  = my_pool_id
       nlevel     = npool
       level_comm = inter_pool_comm
       info       = 'p'
       !
    CASE('B','b')
       !
       mylevelid  = my_bgrp_id
       nlevel     = nbgrp
       level_comm = inter_bgrp_comm
       info       = 'b'
       !
    CASE('Z','z')
       !
       mylevelid  = me_bgrp
       nlevel     = nproc_bgrp
       level_comm = intra_bgrp_comm
       info       = 'z'
       !
    CASE DEFAULT
       CALL errore( 'distr_center', 'Invalid level', 1 )
    END SELECT
    !
    ! Generate nloc, nlocx, nglob
    !
    nglob = n
    nloc  = n/nlevel
    IF( mylevelid .LT. MOD( n, nlevel ) ) THEN
       nloc = nloc + 1
    END IF
    !
    nlocx    = nloc
    nloc_min = nloc
    !
    CALL mp_max(nlocx   , level_comm )
    CALL mp_min(nloc_min, level_comm )
    ! 
    ! Report the distribution across groups
    !
    IF(lverbose) THEN
       CALL io_push_title('Parallelization for '//TRIM(label))
       CALL io_push_value('nglob',           nglob,  20)
       CALL io_push_value('nlevel',          nlevel,  20)
       CALL io_push_value('Min nglob/nlevel',nloc_min,20)
       CALL io_push_value('Max nglob/nlevel',nlocx   ,20)
       CALL io_push_bar()
    ENDIF
    !
    ! Staging distribution into container
    !
    IF(ALLOCATED(dato%l2g)) DEALLOCATE( dato%l2g )
    ALLOCATE( dato%l2g( nloc ) )
    !
    DO iloc = 1, nloc
       CALL distr_l2g(iglob,iloc,mylevelid,nlevel)
       dato%l2g(iloc) = iglob
    ENDDO
    !
    ! Get myoffset
    !
    ALLOCATE( tmp_tab(0:nlevel) )
    !
    tmp_tab = 0 
    tmp_tab( mylevelid ) = nloc
    CALL mp_sum( tmp_tab, level_comm)
    !
    myoffset = 0 
    DO ii = 0, mylevelid - 1
       myoffset = myoffset + tmp_tab(ii)
    ENDDO
    DEALLOCATE( tmp_tab )
    !
    dato%nloc     = nloc
    dato%nlocx    = nlocx
    dato%nglob    = nglob
    dato%myoffset = myoffset 
    dato%mpicomm  = level_comm
    dato%info     = info 
    !
  END SUBROUTINE
  !
  !
  !
  SUBROUTINE distr_g2l(iglob,iloc,who,nlevel)
    !
    IMPLICIT NONE
    !
    INTEGER,INTENT(IN) :: iglob,nlevel
    INTEGER,INTENT(OUT) :: iloc,who
    !
    iloc = ((iglob-1)/(nlevel))+1
    who = MOD( iglob-1, nlevel )
    !
  END SUBROUTINE
  !
  !
  !
  SUBROUTINE distr_l2g(iglob,iloc,mylevelid,nlevel)
    !
    IMPLICIT NONE
    !
    INTEGER,INTENT(IN) :: iloc,mylevelid,nlevel
    INTEGER,INTENT(OUT) :: iglob
    !
    iglob = nlevel*(iloc-1)+mylevelid+1
    !
  END SUBROUTINE
  !
  !
  !
  SUBROUTINE copy_container(dato_in,dato_out)
    !
    IMPLICIT NONE
    !
    TYPE(index_distribution),INTENT(IN)  :: dato_in
    TYPE(index_distribution),INTENT(OUT) :: dato_out
    !
    dato_out%nloc     = dato_in%nloc
    dato_out%nlocx    = dato_in%nlocx
    dato_out%nglob    = dato_in%nglob
    dato_out%myoffset = dato_in%myoffset
    dato_out%mpicomm  = dato_in%mpicomm
    dato_out%info     = dato_out%info
    !
    IF(ALLOCATED(dato_out%l2g)) DEALLOCATE( dato_out%l2g )
    ALLOCATE( dato_out%l2g( dato_in%nloc ) )
    !
    dato_out%l2g(:) = dato_in%l2g(:)
    !
  END SUBROUTINE
  !
  !
  !
  SUBROUTINE destroy_container(dato_in)
    !
    IMPLICIT NONE
    !
    TYPE(index_distribution),INTENT(INOUT)  :: dato_in
    !
    IF(ALLOCATED(dato_in%l2g)) DEALLOCATE( dato_in%l2g )
    !
  END SUBROUTINE
  !
  !
END MODULE
