module io_hpc

#if defined __HDF5
  USE hdf5_pw
#elif defined __ADIOS
#endif


  contains

  subroutine initialize_io_hpc(which)

    implicit none
    integer, intent(in) :: which

    if(which.eq.1)then
      call initialize_hdf5()
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

