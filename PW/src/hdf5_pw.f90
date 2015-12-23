module hdf5_pw
  USE HDF5

  TYPE HDF5_type
   INTEGER(HID_T) :: file_id       ! File identifier 
   INTEGER(HID_T) :: dset_id       ! Dataset identifier 
   INTEGER(HID_T) :: filespace     ! Dataspace identifier in file 
   INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
   INTEGER(HID_T) :: plist_id      ! Property list identifier 
   INTEGER(HSIZE_T) :: counts(2)
   INTEGER(HSIZE_T) :: offset(2)
  END TYPE HDF5_type
  
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

  
  subroutine setup_file_property_hdf5(comm, plist_id, filename, file_id)
   implicit none
   integer, intent(in) :: comm
   integer(HID_T), intent(inout) :: plist_id
   character(len=*), intent(in) :: filename
   integer, intent(inout) :: file_id
   integer :: error, info
   CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error) ! Properties for file creation
   CALL h5pset_fapl_mpio_f(plist_id, comm, info, error) ! Stores MPI IO communicator information to the file access property list
   CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id) ! create the file collectively
   CALL h5pclose_f(plist_id, error)
  end subroutine setup_file_property_hdf5


  subroutine define_dataset_hdf5( rank, memspace, filespace, dset_id, counts ,offset)
   implicit none
   integer(HSIZE_T), intent(in) :: counts(2), offset(2)
   integer(HID_T), intent(inout) :: dset_id, filespace
   integer, intent(in) :: rank
   integer, intent(inout) :: memspace
   integer :: error
   CALL h5screate_simple_f(rank, counts, memspace, error) !define HDF5 dataset
   CALL h5dget_space_f(dset_id, filespace, error)
   CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, counts, error) ! create hyperslab to read from more than 1 proc
  end subroutine define_dataset_hdf5
 
  subroutine  write_data_hdf5(dset_id, data, dimsfi, filespace, memspace, plist_id, counts)
  !subroutine  write_data_hdf5(dset_id,  counts,dimsfi,filespace, memspace, plist_id)
   USE kinds, ONLY : DP
   implicit none
   integer(HID_T), intent(in) :: dset_id, filespace, memspace
   integer(HID_T), intent(inout) ::  plist_id
   integer(HSIZE_T), intent(in) :: counts(2)
   INTEGER(HSIZE_T), intent(in),DIMENSION(2) :: dimsfi 
   real(kind=dp), intent(inout) :: data(:,:)
   !INTEGER, ALLOCATABLE :: data (:,:)  ! Data to write
   integer :: error 
 
   CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
   CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

   
   CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dimsfi, error, &
                  file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)



  end subroutine write_data_hdf5


  subroutine close_fileandprop_hdf5(filespace, memspace, dset_id, plist_id, file_id)
   implicit none
   integer(hid_t), intent(in) :: filespace, memspace, dset_id, file_id, plist_id
   integer error
   ! Close dataspaces.
   !
   CALL h5sclose_f(filespace, error)
   CALL h5sclose_f(memspace, error)

   !
   ! Close the dataset and property list.
   !
   CALL h5dclose_f(dset_id, error)
   CALL h5pclose_f(plist_id, error)

   !
   ! Close the file.
   !
   CALL h5fclose_f(file_id, error)

  end subroutine close_fileandprop_hdf5

 


end module hdf5_pw
