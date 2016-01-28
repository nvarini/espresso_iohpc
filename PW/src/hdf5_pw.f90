module hdf5_pw
  USE HDF5
  USE Kinds, ONLY : DP

  TYPE HDF5_type
   INTEGER(HID_T) :: file_id       ! File identifier 
   INTEGER(HID_T) :: dset_id       ! Dataset identifier 
   INTEGER(HID_T) :: filespace     ! Dataspace identifier in file 
   INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
   INTEGER(HID_T) :: plist_id      ! Property list identifier 
   CHARACTER(LEN=3) :: dsetname  ! Dataset name
   INTEGER(HSIZE_T), DIMENSION(2) :: counts, counts_g, offset
   INTEGER          :: comm
   INTEGER          :: rank
  END TYPE HDF5_type
  TYPE(HDF5_type), save :: evc_hdf5
  INTEGER off_npw, npw_g
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

  
  subroutine setup_file_property_hdf5( hdf5desc ,filename)
   use parallel_include
   implicit none
   type(HDF5_type), intent(inout) :: hdf5desc 
   character(len=*), intent(in) :: filename
   integer(HID_T) :: plist_id
   integer :: error, info
   info = MPI_INFO_NULL
   CALL h5pcreate_f(H5P_FILE_ACCESS_F, hdf5desc%plist_id, error) ! Properties for file creation
   CALL h5pset_fapl_mpio_f(hdf5desc%plist_id, hdf5desc%comm, info, error) ! Stores MPI IO communicator information to the file access property list
   CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, hdf5desc%file_id, error, access_prp = hdf5desc%plist_id) ! create the file collectively
   CALL h5pclose_f(hdf5desc%plist_id, error)

  end subroutine setup_file_property_hdf5


  subroutine define_dataset_hdf5( hdf5desc)
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
  end subroutine define_dataset_hdf5
 
  subroutine  write_data_hdf5(hdf5desc, data)
   USE kinds, ONLY : DP
   USE ISO_C_BINDING
   implicit none
   type(HDF5_type), intent(inout) :: hdf5desc
   complex(kind=dp), intent(inout) :: data(:,:)
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


  subroutine  read_data_hdf5(hdf5desc, data, rank)
   type(HDF5_type), intent(inout) :: hdf5desc
   complex(kind=dp), intent(inout) :: data(:,:)
   integer, intent(in) :: rank
   integer :: error
   TYPE(C_PTR) :: f_ptr

   f_ptr = C_LOC(data(1,1))
   CALL H5dread_f(hdf5desc%dset_id, H5T_NATIVE_DOUBLE, f_ptr, error, &
                   mem_space_id = hdf5desc%memspace, file_space_id = hdf5desc%filespace ,&
                  xfer_prp = hdf5desc%plist_id)

  end subroutine read_data_hdf5



  subroutine close_fileandprop_hdf5(hdf5desc)
   implicit none
   type(HDF5_type), intent(inout) :: hdf5desc
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


  subroutine set_index_hdf5(hdf5desc, var, offset, nglobal,tsize)
    
    USE kinds, only : DP
    implicit none
    COMPLEX(DP), intent(in) :: var(:,:) 
    type(HDF5_type), intent(inout) :: hdf5desc
    INTEGER, intent(in)  :: offset, nglobal,tsize
    
    hdf5desc%counts(1) = size(var(:,1))*tsize
    hdf5desc%counts(2) = size(var(1,:)) 
    hdf5desc%counts_g(1) = nglobal*tsize
    hdf5desc%counts_g(2) = size(var(1,:)) 
    hdf5desc%offset(1) = offset*tsize
    hdf5desc%offset(2) = 0
 
  end subroutine set_index_hdf5



end module hdf5_pw
