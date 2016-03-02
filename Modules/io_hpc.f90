module io_hpc

  USE hdf5_qe
  USE mp_world, ONLY : nproc
  USE io_files, ONLY : wfc_dir, prefix,nd_nmbr, tmp_dir

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
      evc_hdf5%dsetname="evc"
      evc_hdf5%comm=comm
      evc_hdf5%rank =2 
      evc_hdf5%chunk_dim=(/npwx,nbnd/)

      !filename=trim(wfc_dir) //TRIM(prefix) //".wfchdf5"//trim(mpimestring)
      filename=trim(tmp_dir) //TRIM(prefix) //".wfchdf5"//trim(mpimestring)
      if(write.eq..true.)then
        CALL setup_file_property_hdf5(evc_hdf5, filename,.false.,.true.)
      else
        CALL setup_file_property_hdf5(evc_hdf5, filename,.false.,.false.)
      endif
      CALL prepare_index_hdf5(npwx,off_npw,npw_g,evc_hdf5%comm,nproc)
      CALL set_index_hdf5(evc_hdf5,data,off_npw,npw_g,2)
      evc_hdf5%size(1) = npwx*2
      evc_hdf5%size(2) = nbnd
      evc_hdf5%offset(1) = 0
      evc_hdf5%offset(2) = 0
      if(write.eq..true.)CALL define_dataset_hdf5(evc_hdf5,.false.)
      !call errore('','',2)
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

  
end module io_hpc

