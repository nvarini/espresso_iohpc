module hdf5_qe
  USE HDF5
  USE Kinds, ONLY : DP
  USE io_global, ONLY : ionode


  TYPE HDF5_type
   INTEGER(HID_T) :: file_id       ! File identifier 
   INTEGER(HID_T) :: dset_id       ! Dataset identifier 
   INTEGER(HID_T) :: filespace     ! Dataspace identifier in file 
   INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
   INTEGER(HID_T) :: plist_id      ! Property list identifier 
   INTEGER(HID_T) :: group_id      ! Group identifier 
   CHARACTER(LEN=40) :: dsetname  ! Dataset name
   INTEGER          :: rank
   INTEGER(HSIZE_T), DIMENSION(2) :: counts, counts_g, offset
   INTEGER(HSIZE_T), DIMENSION(1:2) :: size
   INTEGER(HID_T) :: crp_list      ! Dataset creation property identifier 
   INTEGER          :: comm
   INTEGER(HSIZE_T), DIMENSION(1:2) :: maxdims
   INTEGER(HSIZE_T), DIMENSION(1:2) :: chunk_dim
   character(len=256) filename
  END TYPE HDF5_type

  TYPE(HDF5_type), save :: evc_hdf5, evc_hdf5_write
  INTEGER, save ::  off_npw, npw_g




  contains

  subroutine initialize_hdf5()
    implicit none
    integer :: error
    call h5open_f(error)
  end subroutine initialize_hdf5

  subroutine finalize_hdf5(hdf5desc)
    implicit none
    type(HDF5_type), intent(in) :: hdf5desc
    integer :: error
    
    call h5pclose_f(hdf5desc%plist_id,error)
    call h5close_f(error)
  end subroutine finalize_hdf5

  
  subroutine setup_file_property_hdf5(hdf5desc ,filename, hyperslab, write, kpoint)
   use parallel_include
   use mp_world,  only : mpime
   implicit none
   type(HDF5_type), intent(inout) :: hdf5desc 
   character(len=*), intent(inout) :: filename
   logical,  intent(in) :: hyperslab, write
   integer, intent(in) ::  kpoint
   integer(HID_T) :: plist_id
   integer :: error, info
   character*4 kstring
   write(kstring,'(I0)') kpoint


   info = MPI_INFO_NULL



   if(hyperslab .eq. .true. ) then
     CALL h5pcreate_f(H5P_FILE_ACCESS_F, hdf5desc%plist_id, error) ! Properties for file creation
     CALL h5pset_fapl_mpio_f(hdf5desc%plist_id, hdf5desc%comm, info, error) ! Stores MPI IO communicator information to the file access property list
     CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, hdf5desc%file_id, error, access_prp = hdf5desc%plist_id) ! create the file collectively
     CALL h5pclose_f(hdf5desc%plist_id, error)
   else

    if(write.eq..true.)then
      
      if(kpoint.eq.1)then
        !CALL h5pcreate_f(H5P_FILE_ACCESS_F, hdf5desc%plist_id, error)
        !CALL h5pset_fapl_mpio_f(hdf5desc%plist_id, hdf5desc%comm, info, error)
        !CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, hdf5desc%file_id, error, access_prp=hdf5desc%plist_id) ! create the file collectively
        CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, hdf5desc%file_id, error) 
      endif
    else
      if(ionode)then
        CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, hdf5desc%file_id, error) ! create the file collectively
        !CALL h5dopen_f(hdf5desc%file_id, hdf5desc%dsetname, hdf5desc%dset_id, error)
      endif
    endif
   endif
   
  end subroutine setup_file_property_hdf5


  subroutine define_dataset_hdf5_hyperslab(hdf5desc)
    USE mp_world,             ONLY : mpime
   implicit none
   type(HDF5_type), intent(inout) :: hdf5desc
   integer :: error

     CALL h5screate_simple_f(hdf5desc%rank, hdf5desc%counts_g, hdf5desc%filespace, error) !define HDF5 dataset
     CALL h5dcreate_f(hdf5desc%file_id, hdf5desc%dsetname, H5T_NATIVE_DOUBLE, hdf5desc%filespace, &
                      hdf5desc%dset_id, error)
     CALL h5sclose_f(hdf5desc%filespace, error)

     CALL h5screate_simple_f(hdf5desc%rank, hdf5desc%counts, hdf5desc%memspace, error) 
     CALL h5dget_space_f(hdf5desc%dset_id, hdf5desc%filespace, error)
     CALL h5sselect_hyperslab_f(hdf5desc%filespace, H5S_SELECT_SET_F, hdf5desc%offset, hdf5desc%counts, error) ! create hyperslab to read from more than 1 proc
   end subroutine define_dataset_hdf5_hyperslab


 
  subroutine  write_data_hdf5(hdf5desc, data,  niter)
   USE kinds, ONLY : DP
   USE ISO_C_BINDING
   USE mp_world, ONLY : mpime
   implicit none
   type(HDF5_type), intent(inout) :: hdf5desc
   complex(kind=dp), intent(inout) :: data(:,:)
   integer, intent(in), optional :: niter
   !INTEGER, ALLOCATABLE :: data (:,:)  ! Data to write
   integer :: error 
   real(kind=dp)   :: tmp
   integer(HID_T)     :: complex_id, double_id
   integer(HSIZE_T)   :: double_size, complex_size
   TYPE(C_PTR) :: f_ptr
   
   

     CALL h5pcreate_f(H5P_DATASET_XFER_F, hdf5desc%plist_id, error)
     CALL h5pset_dxpl_mpio_f(hdf5desc%plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    
     f_ptr = C_LOC(data(1,1))
     CALL h5dwrite_f(hdf5desc%dset_id, H5T_NATIVE_DOUBLE, f_ptr, error,&
                    file_space_id = hdf5desc%filespace, mem_space_id = hdf5desc%memspace, &
                    xfer_prp = hdf5desc%plist_id)
  end subroutine write_data_hdf5






  subroutine  read_data_hdf5(hdf5desc, data)
   type(HDF5_type), intent(inout) :: hdf5desc
   complex(kind=dp), intent(inout) :: data(:,:)
   integer :: error
   TYPE(C_PTR) :: f_ptr

     f_ptr = C_LOC(data(1,1))
     CALL H5dread_f(hdf5desc%dset_id, H5T_NATIVE_DOUBLE, f_ptr, error, &
                     mem_space_id = hdf5desc%memspace, file_space_id = hdf5desc%filespace ,&
                    xfer_prp = hdf5desc%plist_id)
  end subroutine read_data_hdf5

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




  subroutine prepare_for_writing_final(hdf5desc,comm,filename_input,kpoint)
    USE mp_world,             ONLY : mpime
    USE io_files, ONLY : wfc_dir, prefix, tmp_dir
    implicit none
    type(HDF5_type), intent(inout) :: hdf5desc
    character(len=*), intent(in):: filename_input
    integer, intent(in) :: comm, kpoint
    character(len=256) filename
    character*4 mpimestring
    integer :: ik, error
    character*4 kstring
    write(kstring,'(I0)') kpoint
    
    
 
    hdf5desc%comm=comm
    !hdf5desc%filename = trim(filename_input) //".wfchdf5"
    hdf5desc%filename=trim(tmp_dir) //"evc.hdf5"!//trim(mpimestring)
    
    CALL setup_file_property_hdf5(hdf5desc,hdf5desc%filename ,.false.,.true.,kpoint)
    
    if(kpoint>1) CALL h5fopen_f(hdf5desc%filename, H5F_ACC_RDWR_F, hdf5desc%file_id, error) ! create the file collectively
    CALL h5gcreate_f(hdf5desc%file_id, kstring, hdf5desc%group_id, error)
    CALL h5gclose_f(hdf5desc%group_id, error)

  end subroutine prepare_for_writing_final



  subroutine prepare_for_reading_final(hdf5desc,comm,filename_input,kpoint)
    USE mp_world,             ONLY : mpime
    USE io_files,             ONLY : wfc_dir, prefix, tmp_dir
    implicit none
    type(HDF5_type), intent(inout) :: hdf5desc
    character(len=*), intent(in):: filename_input
    integer, intent(in) :: comm,kpoint
    character(len=256) filename
    character*4 mpimestring
    integer :: ik
 
    hdf5desc%comm=comm
    hdf5desc%rank =1 
    !filename = trim(filename_input) //".wfchdf5"
    filename=trim(tmp_dir) //"evc.hdf5"!//trim(mpimestring)
    write(mpime+100,*) filename
    CALL setup_file_property_hdf5(hdf5desc,filename ,.false.,.false.,kpoint)

  end subroutine prepare_for_reading_final





  subroutine write_final_data(hdf5desc,dsetname,var,kpoint)
    USE kinds, ONLY : DP
    USE mp_world, ONLY : mpime
    implicit none
    type(HDF5_type), intent(inout) :: hdf5desc
    integer, intent(in) :: dsetname, kpoint
    complex(kind=DP), intent(in) :: var(:)
    INTEGER(HID_T) :: dspace_id, dset_id     ! Dataspace identifier
    integer :: error
    INTEGER(HSIZE_T), DIMENSION(1) :: counts
    character*4 dset_name
    TYPE(C_PTR) :: f_ptr
    character*4 kstring
    write(kstring,'(I0)') kpoint

    write(dset_name,'(I0)') dsetname
    counts=size(var)*2  
    CALL h5screate_simple_f(1, counts, dspace_id, error) !create the dataspace
    CALL h5gopen_f(hdf5desc%file_id,kstring,hdf5desc%group_id,error)
    CALL h5dcreate_f(hdf5desc%group_id, dset_name, H5T_NATIVE_DOUBLE, dspace_id, &
                      dset_id, error)
    f_ptr = C_LOC(var(1))
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, f_ptr, error)
    CALL h5dclose_f(dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5gclose_f(hdf5desc%group_id, error)
  end subroutine write_final_data

  subroutine read_final_data(hdf5desc,dsetname,var,kpoint)
    USE kinds, ONLY : DP
    USE mp_world, ONLY : mpime
    implicit none
    type(HDF5_type), intent(inout) :: hdf5desc
    integer, intent(in) :: dsetname, kpoint
    complex(kind=DP), intent(inout) :: var(:)
    INTEGER(HID_T) :: dspace_id, dset_id     ! Dataspace identifier
    integer :: error
    INTEGER(HSIZE_T), DIMENSION(1) :: counts
    character*4 dset_name
    TYPE(C_PTR) :: f_ptr
    character*4 kstring
    character*100 errmsg

    write(dset_name,'(I0)') dsetname
    write(kstring,'(I0)') kpoint

    write(mpime+200,*) dset_name, kstring
    counts=size(var)*2  
    CALL h5gopen_f(hdf5desc%file_id,kstring,hdf5desc%group_id,error)
    if(error.ne.0) call errore('error in h5gopen_f','',error)
    CALL h5dopen_f(hdf5desc%group_id, dset_name, dset_id, error)
    if(error.ne.0) call errore('error in h5dopen_f','',error)
    CALL h5dget_type_f(dset_id, dspace_id, error)
    f_ptr = C_LOC(var(1))
    !CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, f_ptr, error)
    CALL h5dread_f(dset_id, dspace_id, f_ptr, error)
    CALL h5dclose_f(dset_id, error)
    CALL h5gclose_f(hdf5desc%group_id, error)
  end subroutine read_final_data

  subroutine initialize_io_hdf5(hdf5desc,comm, data, write,kpoint)
    USE io_files, ONLY : wfc_dir, prefix, tmp_dir
    USE mp_world, ONLY : mpime
    USE kinds,    ONLY : dp
    USE mp_world, ONLY : nproc
    implicit none
    
    TYPE(HDF5_type), intent(inout) :: hdf5desc
    complex(kind=dp), intent(in) :: data(:,:)
    integer, intent(in) :: comm, kpoint
    logical, intent(in) :: write
    character(len=80) :: filename
    character*4 mpimestring
    integer :: npwx, nbnd
    write(mpimestring,'(I0)') mpime
    npwx=size(data(:,1))
    nbnd=size(data(1,:))
    
    filename=trim(tmp_dir) //TRIM(prefix) //".wfchdf5"
    call initialize_hdf5_array(hdf5desc,comm,npwx,nbnd)
    if(write.eq..true.)then
      CALL setup_file_property_hdf5(hdf5desc, filename,.true.,.true.,kpoint)
    else
      CALL setup_file_property_hdf5(hdf5desc, filename,.false.,.false.,kpoint)
    endif
    CALL prepare_index_hdf5(npwx,off_npw,npw_g,hdf5desc%comm,nproc)
    CALL set_index_hdf5(hdf5desc,data,off_npw,npw_g,2)
    if(write.eq..true.)CALL define_dataset_hdf5_hyperslab(hdf5desc)

  end subroutine initialize_io_hdf5

  subroutine initialize_hdf5_array(hdf5desc,comm,n1,n2)
    use mp_world, only : mpime 
    implicit none
    integer, intent(in) :: n1, n2, comm
    type(HDF5_type), intent(inout) ::  hdf5desc
    hdf5desc%dsetname="evc"
    hdf5desc%comm=comm
    hdf5desc%rank =2 
    hdf5desc%chunk_dim=(/n1,n2/)
    hdf5desc%size(1) = n1*2
    hdf5desc%size(2) = n2
    hdf5desc%offset(1) = 0
    hdf5desc%offset(2) = 0
   
  end subroutine initialize_hdf5_array

  SUBROUTINE set_index_hdf5(hdf5desc, var, offset, nglobal,tsize)
    
    USE kinds, only : DP
    implicit none
    COMPLEX(DP), intent(in) :: var(:,:) 
    type(HDF5_type), intent(inout) :: hdf5desc
    INTEGER, intent(in)  :: offset, nglobal,tsize
    
    hdf5desc%counts(1)   = size(var(:,1))*tsize
    hdf5desc%counts(2)   = size(var(1,:)) 
    hdf5desc%counts_g(1) = nglobal*tsize
    hdf5desc%counts_g(2) = size(var(1,:)) 
    hdf5desc%offset(1)   = offset*tsize
    hdf5desc%offset(2)   = 0
 
  END SUBROUTINE set_index_hdf5





end module hdf5_qe
