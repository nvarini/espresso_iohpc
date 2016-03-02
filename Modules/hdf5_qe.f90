module hdf5_qe
  USE HDF5
  USE Kinds, ONLY : DP
  USE io_global, ONLY : ionode


  TYPE HDF5_type_2d
   INTEGER(HID_T) :: file_id       ! File identifier 
   INTEGER(HID_T) :: dset_id       ! Dataset identifier 
   INTEGER(HID_T) :: filespace     ! Dataspace identifier in file 
   INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
   INTEGER(HID_T) :: plist_id      ! Property list identifier 
   CHARACTER(LEN=40) :: dsetname  ! Dataset name
   INTEGER          :: rank
   INTEGER(HSIZE_T), DIMENSION(2) :: counts, counts_g, offset
   INTEGER(HSIZE_T), DIMENSION(1:2) :: size
   INTEGER(HID_T) :: crp_list      ! Dataset creation property identifier 
   INTEGER          :: comm
   INTEGER(HSIZE_T), DIMENSION(1:2) :: maxdims
   INTEGER(HSIZE_T), DIMENSION(1:2) :: chunk_dim
  END TYPE HDF5_type_2d

  TYPE(HDF5_type_2d), save :: evc_hdf5, evc_hdf5_write
  INTEGER, save ::  off_npw, npw_g

  INTERFACE write_data_hdf5
    MODULE procedure write_data_hdf5_2d,&
      write_data_hdf5_1d
  END INTERFACE write_data_hdf5

  INTERFACE extend_dataset_hdf5
    MODULE procedure extend_dataset_hdf5_1d, &
      extend_dataset_hdf5_2d
  END INTERFACE extend_dataset_hdf5

  INTERFACE set_index_hdf5
    MODULE procedure set_index_hdf5_1d, &
      set_index_hdf5_2d
  END INTERFACE set_index_hdf5

  INTERFACE read_data_hdf5
    MODULE procedure read_data_hdf5_2d,&
      read_data_hdf5_1d
  END INTERFACE read_data_hdf5



  contains

  subroutine initialize_hdf5()
    implicit none
    integer :: error
    call h5open_f(error)
  end subroutine initialize_hdf5

  subroutine finalize_hdf5()
    implicit none
    integer :: error
   
    call h5close_f(error)
  end subroutine finalize_hdf5

  
  subroutine setup_file_property_hdf5(hdf5desc ,filename, hyperslab, write)
   use parallel_include
   use mp_world,  only : mpime
   implicit none
   type(HDF5_type_2d), intent(inout) :: hdf5desc 
   character(len=*), intent(inout) :: filename
   logical,  intent(in) :: hyperslab, write
   integer(HID_T) :: plist_id
   integer :: error, info
   info = MPI_INFO_NULL
   if(hyperslab .eq. .true. ) then
     CALL h5pcreate_f(H5P_FILE_ACCESS_F, hdf5desc%plist_id, error) ! Properties for file creation
     CALL h5pset_fapl_mpio_f(hdf5desc%plist_id, hdf5desc%comm, info, error) ! Stores MPI IO communicator information to the file access property list
     CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, hdf5desc%file_id, error, access_prp = hdf5desc%plist_id) ! create the file collectively
     CALL h5pclose_f(hdf5desc%plist_id, error)
   else
    if(write.eq..true.)then
      CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, hdf5desc%file_id, error) ! create the file collectively
    else
      if(ionode)then
        CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, hdf5desc%file_id, error) ! create the file collectively
        CALL h5dopen_f(hdf5desc%file_id, hdf5desc%dsetname, hdf5desc%dset_id, error)
      endif
    endif
   endif
   

  end subroutine setup_file_property_hdf5


  subroutine define_dataset_hdf5(hdf5desc, hyperslab)
    USE mp_world,             ONLY : mpime
   implicit none
   type(HDF5_type_2d), intent(inout) :: hdf5desc
   logical, intent(in) :: hyperslab
   integer :: error

   if(hyperslab.eq..true.)then
     CALL h5screate_simple_f(hdf5desc%rank, hdf5desc%counts_g, hdf5desc%filespace, error) !define HDF5 dataset
     CALL h5dcreate_f(hdf5desc%file_id, hdf5desc%dsetname, H5T_NATIVE_DOUBLE, hdf5desc%filespace, &
                      hdf5desc%dset_id, error)
     CALL h5sclose_f(hdf5desc%filespace, error)

     CALL h5screate_simple_f(hdf5desc%rank, hdf5desc%counts, hdf5desc%memspace, error) 
     CALL h5dget_space_f(hdf5desc%dset_id, hdf5desc%filespace, error)
     CALL h5sselect_hyperslab_f(hdf5desc%filespace, H5S_SELECT_SET_F, hdf5desc%offset, hdf5desc%counts, error) ! create hyperslab to read from more than 1 proc
   else

     if(hdf5desc%rank.eq.2)then
       hdf5desc%maxdims = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
     else 
       hdf5desc%maxdims = (/H5S_UNLIMITED_F,1/)
     endif
     !Modify dataset creation properties, i.e. enable chunking
     !
  !
     CALL h5pcreate_f(H5P_DATASET_CREATE_F, hdf5desc%crp_list, error)

    CALL h5pset_chunk_f(hdf5desc%crp_list, hdf5desc%rank, hdf5desc%chunk_dim(1), error)

     CALL h5screate_simple_f(hdf5desc%rank, hdf5desc%counts(1), hdf5desc%filespace, error,hdf5desc%maxdims(1)) !define HDF5 dataset
     CALL h5dcreate_f(hdf5desc%file_id, hdf5desc%dsetname, H5T_NATIVE_DOUBLE, hdf5desc%filespace, &
                      hdf5desc%dset_id, error, hdf5desc%crp_list)
     CALL h5sclose_f(hdf5desc%filespace, error)
   endif 
  end subroutine define_dataset_hdf5

  subroutine extend_dataset_hdf5_2d(hdf5desc,var,extend,tsize)  
    USE kinds, only : DP
    USE mp_world, only : mpime
    implicit none
    COMPLEX(DP), intent(in) :: var(:,:) 
    type(HDF5_type_2d), intent(inout) :: hdf5desc
    INTEGER, intent(in)  ::  extend,tsize
    integer :: error
    
    hdf5desc%size(1) = hdf5desc%size(1) + extend*tsize
    hdf5desc%size(2) =  size(var(1,:))
    hdf5desc%counts(1)   = extend*tsize
    hdf5desc%counts(2)   = size(var(1,:)) 
    hdf5desc%offset(1)   =  hdf5desc%size(1) - hdf5desc%counts(1)
    hdf5desc%offset(2)   = 0
 
    call h5dset_extent_f(hdf5desc%dset_id,hdf5desc%size,error)
    CALL h5screate_simple_f (hdf5desc%rank, hdf5desc%counts, hdf5desc%memspace, error)
    CALL h5dget_space_f(hdf5desc%dset_id, hdf5desc%filespace, error)
    CALL h5sselect_hyperslab_f(hdf5desc%filespace, H5S_SELECT_SET_F, &
         hdf5desc%offset, hdf5desc%counts, error)

  end subroutine extend_dataset_hdf5_2d

  subroutine extend_dataset_hdf5_1d(hdf5desc,var,extend,tsize, write)  
    USE kinds, only : DP
    USE mp_world, only : mpime
    implicit none
    COMPLEX(DP), intent(in) :: var(:) 
    type(HDF5_type_2d), intent(inout) :: hdf5desc
    INTEGER, intent(in)  ::  extend,tsize
    logical, intent(in)  :: write
    integer :: error
    
    hdf5desc%size(1) = hdf5desc%size(1) + extend*tsize
    hdf5desc%counts(1)   = extend*tsize
    hdf5desc%offset(1)   =  hdf5desc%size(1) - hdf5desc%counts(1)
 
    if(write.eq..true.)then
      call h5dset_extent_f(hdf5desc%dset_id,hdf5desc%size,error)
      CALL h5screate_simple_f (hdf5desc%rank, hdf5desc%counts, hdf5desc%memspace, error)
    endif
    CALL h5dget_space_f(hdf5desc%dset_id, hdf5desc%filespace, error)
    CALL h5sselect_hyperslab_f(hdf5desc%filespace, H5S_SELECT_SET_F, &
         hdf5desc%offset, hdf5desc%counts, error)

  end subroutine extend_dataset_hdf5_1d


 
  subroutine  write_data_hdf5_2d(hdf5desc, data, hyperslab, niter)
   USE kinds, ONLY : DP
   USE ISO_C_BINDING
   USE mp_world, ONLY : mpime
   implicit none
   type(HDF5_type_2d), intent(inout) :: hdf5desc
   complex(kind=dp), intent(inout) :: data(:,:)
   logical, intent(in) :: hyperslab
   integer, intent(in), optional :: niter
   !INTEGER, ALLOCATABLE :: data (:,:)  ! Data to write
   integer :: error 
   real(kind=dp)   :: tmp
   integer(HID_T)     :: complex_id, double_id
   integer(HSIZE_T)   :: double_size, complex_size
   TYPE(C_PTR) :: f_ptr
   
   

   if(hyperslab.eq..true.)then 
     CALL h5pcreate_f(H5P_DATASET_XFER_F, hdf5desc%plist_id, error)
     CALL h5pset_dxpl_mpio_f(hdf5desc%plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    
     f_ptr = C_LOC(data(1,1))
     CALL h5dwrite_f(hdf5desc%dset_id, H5T_NATIVE_DOUBLE, f_ptr, error,&
                    file_space_id = hdf5desc%filespace, mem_space_id = hdf5desc%memspace, &
                    xfer_prp = hdf5desc%plist_id)
   else
     f_ptr = C_LOC(data(1,1))
     if(niter.eq.1)then
       CALL h5dwrite_f(hdf5desc%dset_id, H5T_NATIVE_DOUBLE, f_ptr, error)
     else
       CALL h5dwrite_f(hdf5desc%dset_id, H5T_NATIVE_DOUBLE, f_ptr, error,&
                    hdf5desc%memspace,hdf5desc%filespace)
     endif
 
   endif

  end subroutine write_data_hdf5_2d



  subroutine  write_data_hdf5_1d(hdf5desc, data, hyperslab, niter)
   USE kinds, ONLY : DP
   USE ISO_C_BINDING
   USE mp_world, ONLY : mpime
   implicit none
   type(HDF5_type_2d), intent(inout) :: hdf5desc
   complex(kind=dp), intent(inout) :: data(:)
   logical, intent(in) :: hyperslab
   integer, intent(in), optional :: niter
   !INTEGER, ALLOCATABLE :: data (:,:)  ! Data to write
   integer :: error 
   real(kind=dp)   :: tmp
   integer(HID_T)     :: complex_id, double_id
   integer(HSIZE_T)   :: double_size, complex_size
   TYPE(C_PTR) :: f_ptr
   
   

   if(hyperslab.eq..true.)then 
     CALL h5pcreate_f(H5P_DATASET_XFER_F, hdf5desc%plist_id, error)
     CALL h5pset_dxpl_mpio_f(hdf5desc%plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    
     f_ptr = C_LOC(data(1))
     CALL h5dwrite_f(hdf5desc%dset_id, H5T_NATIVE_DOUBLE, f_ptr, error,&
                    file_space_id = hdf5desc%filespace, mem_space_id = hdf5desc%memspace, &
                    xfer_prp = hdf5desc%plist_id)
   else
     f_ptr = C_LOC(data(1))
     if(niter.eq.1)then
       CALL h5dwrite_f(hdf5desc%dset_id, H5T_NATIVE_DOUBLE, f_ptr, error)
     else
       CALL h5dwrite_f(hdf5desc%dset_id, H5T_NATIVE_DOUBLE, f_ptr, error,&
                    hdf5desc%memspace,hdf5desc%filespace)
     endif
 
   endif

  end subroutine write_data_hdf5_1d




  subroutine  read_data_hdf5_2d(hdf5desc, data, hyperslab)
   type(HDF5_type_2d), intent(inout) :: hdf5desc
   complex(kind=dp), intent(inout) :: data(:,:)
   logical, intent(in) :: hyperslab
   integer :: error
   TYPE(C_PTR) :: f_ptr

   if(hyperslab.eq..true.)then
     f_ptr = C_LOC(data(1,1))
     CALL H5dread_f(hdf5desc%dset_id, H5T_NATIVE_DOUBLE, f_ptr, error, &
                     mem_space_id = hdf5desc%memspace, file_space_id = hdf5desc%filespace ,&
                    xfer_prp = hdf5desc%plist_id)
   else
     f_ptr = C_LOC(data(1,1))

     CALL h5dget_space_f(hdf5desc%dset_id, hdf5desc%filespace, error)
     CALL h5sselect_hyperslab_f(hdf5desc%filespace, H5S_SELECT_SET_F, &
         hdf5desc%offset, hdf5desc%counts, error)
     CALL h5sget_simple_extent_ndims_f(hdf5desc%filespace, hdf5desc%rank, error)
     CALL h5screate_simple_f(hdf5desc%rank, hdf5desc%counts, hdf5desc%memspace, error)

     CALL H5dread_f(hdf5desc%dset_id, H5T_NATIVE_DOUBLE, f_ptr,  &
       error, hdf5desc%memspace, hdf5desc%filespace)
   endif

  end subroutine read_data_hdf5_2d

  subroutine  read_data_hdf5_1d(hdf5desc, data, hyperslab)
   type(HDF5_type_2d), intent(inout) :: hdf5desc
   complex(kind=dp), intent(inout) :: data(:)
   logical, intent(in) :: hyperslab
   integer :: error
   TYPE(C_PTR) :: f_ptr

   if(hyperslab.eq..true.)then
     f_ptr = C_LOC(data(1))
     CALL H5dread_f(hdf5desc%dset_id, H5T_NATIVE_DOUBLE, f_ptr, error, &
                     mem_space_id = hdf5desc%memspace, file_space_id = hdf5desc%filespace ,&
                    xfer_prp = hdf5desc%plist_id)
   else
     f_ptr = C_LOC(data(1))

     CALL h5dget_space_f(hdf5desc%dset_id, hdf5desc%filespace, error)
     CALL h5sselect_hyperslab_f(hdf5desc%filespace, H5S_SELECT_SET_F, &
         hdf5desc%offset(1), hdf5desc%counts(1), error)
     CALL h5sget_simple_extent_ndims_f(hdf5desc%filespace, hdf5desc%rank, error)
     CALL h5screate_simple_f(hdf5desc%rank, hdf5desc%counts(1), hdf5desc%memspace, error)


     CALL H5dread_f(hdf5desc%dset_id, H5T_NATIVE_DOUBLE, f_ptr,  &
       error, hdf5desc%memspace, hdf5desc%filespace)
   endif

  end subroutine read_data_hdf5_1d


  subroutine close_fileandprop_hdf5(hdf5desc)
   implicit none
   type(HDF5_type_2d), intent(inout) :: hdf5desc
   integer error
   ! Close dataspaces.
   !
   CALL h5sclose_f(hdf5desc%filespace, error)
   CALL h5sclose_f(hdf5desc%memspace, error)

   !
   ! Close the dataset and property list.
   !
   CALL h5dclose_f(hdf5desc%dset_id, error)
   CALL h5pclose_f(hdf5desc%plist_id, error)

   !
   ! Close the file.
   !
   CALL h5fclose_f(hdf5desc%file_id, error)

  end subroutine close_fileandprop_hdf5

  SUBROUTINE prepare_index_hdf5(sendm,recm,globalm,comm,nproc)

   USE parallel_include
   USE mp, ONLY : mp_sum
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: comm, nproc
   INTEGER, INTENT(INOUT) :: sendm, recm, globalm
   INTEGER :: errore

   call mpi_scan(sendm,recm,1,MPI_INTEGER,MPI_SUM,comm,errore)
   recm=recm-sendm
   globalm=sendm
   call mp_sum(globalm,comm)


 END SUBROUTINE prepare_index_hdf5


  SUBROUTINE set_index_hdf5_2d(hdf5desc, var, offset, nglobal,tsize)
    
    USE kinds, only : DP
    implicit none
    COMPLEX(DP), intent(in) :: var(:,:) 
    type(HDF5_type_2d), intent(inout) :: hdf5desc
    INTEGER, intent(in)  :: offset, nglobal,tsize
    
    hdf5desc%counts(1)   = size(var(:,1))*tsize
    hdf5desc%counts(2)   = size(var(1,:)) 
    hdf5desc%counts_g(1) = nglobal*tsize
    hdf5desc%counts_g(2) = size(var(1,:)) 
    hdf5desc%offset(1)   = offset*tsize
    hdf5desc%offset(2)   = 0
 
  END SUBROUTINE set_index_hdf5_2d

  SUBROUTINE set_index_hdf5_1d(hdf5desc, var, offset, nglobal,tsize)
    
    USE kinds, only : DP
    implicit none
    COMPLEX(DP), intent(in) :: var(:) 
    type(HDF5_type_2d), intent(inout) :: hdf5desc
    INTEGER, intent(in)  :: offset, nglobal,tsize
    
    hdf5desc%counts(1)   = size(var)*tsize
    hdf5desc%counts_g(1) = nglobal*tsize
    hdf5desc%offset(1)   = offset*tsize
 
  END SUBROUTINE set_index_hdf5_1d


  subroutine prepare_for_writing(hdf5desc,comm,chunk_dim,filename_input)
    USE mp_world,             ONLY : mpime
    implicit none
    type(HDF5_type_2d), intent(inout) :: hdf5desc
    character(len=*), intent(in):: filename_input
    integer, intent(in) :: comm, chunk_dim
    character(len=256) filename
    character*4 mpimestring
    integer :: ik
 
    hdf5desc%dsetname="evcwrite"
    hdf5desc%comm=comm
    hdf5desc%rank =1 
    hdf5desc%chunk_dim=(/chunk_dim,1/)

    filename = trim(filename_input) //".wfchdf5"
    CALL setup_file_property_hdf5(hdf5desc,filename ,.false.,.true.)
    !call errore('','',1)
    !CALL prepare_index_hdf5(npwx,off_npw,npw_g,hdf5desc%comm,nproc)
    !CALL set_index_hdf5(evc_hdf5,evc,off_npw,npw_g,2)
    !evc_hdf5%size(1) = npwx*2
    hdf5desc%offset(1) = 0
    hdf5desc%counts(1) = chunk_dim
    hdf5desc%size(1) = chunk_dim
    CALL define_dataset_hdf5(hdf5desc,.false.)

  end subroutine prepare_for_writing

  subroutine prepare_for_reading(hdf5desc,comm,chunk_dim,filename_input)
    USE mp_world,             ONLY : mpime
    implicit none
    type(HDF5_type_2d), intent(inout) :: hdf5desc
    character(len=*), intent(in):: filename_input
    integer, intent(in) :: comm, chunk_dim
    character(len=256) filename
    character*4 mpimestring
    integer :: ik
 
    hdf5desc%dsetname="evcwrite"
    hdf5desc%comm=comm
    hdf5desc%rank =1 
    hdf5desc%chunk_dim=(/chunk_dim,1/)

    filename = trim(filename_input) //".wfchdf5"
    CALL setup_file_property_hdf5(hdf5desc,filename ,.false.,.false.)
    !CALL prepare_index_hdf5(npwx,off_npw,npw_g,hdf5desc%comm,nproc)
    !CALL set_index_hdf5(evc_hdf5,evc,off_npw,npw_g,2)
    !evc_hdf5%size(1) = npwx*2
    hdf5desc%offset(1) = 0
    hdf5desc%counts(1) = chunk_dim
    hdf5desc%size(1) = chunk_dim
    !CALL define_dataset_hdf5(hdf5desc,.false.)
    !call errore('','',1)

  end subroutine prepare_for_reading




end module hdf5_qe
