module io_hpc

  USE hdf5_pw
  USE mp_world, ONLY : nproc
  USE wvfct,    ONLY : npwx
  USE wavefunctions_module,ONLY : evc

  contains

  subroutine initialize_io_hpc(which, comm)

    implicit none
    integer, intent(in) :: which, comm

    if(which.eq.1)then
      call initialize_hdf5()
      evc_hdf5%dsetname="evc"
      evc_hdf5%comm=comm
      evc_hdf5%rank = 2 
      CALL setup_file_property_hdf5(evc_hdf5, "evc.hdf5")
      CALL prepare_index_hdf5(npwx,off_npw,npw_g,evc_hdf5%comm,nproc)
      CALL set_index_hdf5(evc_hdf5,evc,off_npw,npw_g,2)
      CALL define_dataset_hdf5(evc_hdf5)
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

