module io_hpc

  USE hdf5_qe
  USE mp_world, ONLY : nproc
  USE io_files, ONLY : wfc_dir, prefix, tmp_dir

  contains

  subroutine initialize_io_hpc(which, comm, data, write)
    
    USE mp_world, ONLY : mpime
    USE kinds,    ONLY : dp
    implicit none
    
    complex(kind=dp), intent(in) :: data(:,:)
    integer, intent(in) :: which, comm
    logical, intent(in) :: write
    character(len=80) :: filename
    character*4 mpimestring
    integer :: npwx, nbnd
    
    write(mpimestring,'(I0)') mpime
    npwx=size(data(:,1))
    nbnd=size(data(1,:))
   
    if(which.eq.1)then
      call initialize_hdf5()
      call initialize_hdf5_array(evc_hdf5,comm,npwx,nbnd)
      !filename=trim(wfc_dir) //TRIM(prefix) //".wfchdf5"//trim(mpimestring)
      filename=trim(tmp_dir) //TRIM(prefix) //".wfchdf5"!//trim(mpimestring)
      if(write.eq..true.)then
        CALL setup_file_property_hdf5(evc_hdf5, filename,.true.,.true.)
      else
        CALL setup_file_property_hdf5(evc_hdf5, filename,.false.,.false.)
      endif
      CALL prepare_index_hdf5(npwx,off_npw,npw_g,evc_hdf5%comm,nproc)
      CALL set_index_hdf5(evc_hdf5,data,off_npw,npw_g,2)
      if(write.eq..true.)CALL define_dataset_hdf5(evc_hdf5,.true.)
    else
    endif
  end subroutine initialize_io_hpc

  subroutine finalize_io_hpc(which)

    implicit none
    integer, intent(in) :: which

    if(which.eq.1)then
      call finalize_hdf5()
    else
    endif
  
  end subroutine finalize_io_hpc  
   
  subroutine initialize_hdf5_array(hdf5desc,comm,n1,n2)
  
    implicit none
    integer, intent(in) :: n1, n2, comm
    type(HDF5_type) hdf5desc
    hdf5desc%dsetname="evc"
    hdf5desc%comm=comm
    hdf5desc%rank =2 
    hdf5desc%chunk_dim=(/n1,n2/)
    hdf5desc%size(1) = n1*2
    hdf5desc%size(2) = n2
    hdf5desc%offset(1) = 0
    hdf5desc%offset(2) = 0
   
  end subroutine initialize_hdf5_array
  
end module io_hpc

