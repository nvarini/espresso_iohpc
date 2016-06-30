# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_HDF5], [

  AC_MSG_CHECKING([hdf5])
 
AC_ARG_WITH(hdf5,
   [AS_HELP_STRING([--with-hdf5],
       [(yes|no|<path>) Use HDF5. Self-compile or a <path> can be specified (default: no)])],
   [if  test "$withval" = "yes" ; then
      with_hdf5=1
   elif  test "$withval" = "no" ; then
      with_hdf5=0
   else
      with_hdf5=2
      with_hdf5_path="$withval"
   fi],
   [with_hdf5=0])

hdf5_libs=""
 
hdf5_libs_switch="disabled"
if test "$with_hdf5" -eq 2 ; then
     hdf5_libs="-L$with_hdf5_path/lib -lhdf5_fortran -lhdf5"
     try_iflags="$try_iflags -I$with_hdf5_path/include"
     try_dflags="$try_dflags -D__HDF5"
     FFLAGS="\$(FFLAGS) -I$with_hdf5_path/include"
fi

if test "$with_hdf5" -eq 1 ; then
    hdf5_libs="-L\$(TOPDIR)/HDF5/lib -lhdf5_fortran -lhdf5"
    try_dflags="$try_dflags -D__HDF5"
    hdf5_libs_switch="enabled"
    FFLAGS="\$(FFLAGS) -I$with_hdf5_path/include"
fi

  AC_MSG_RESULT(${hdf5_libs})
  AC_SUBST(hdf5_libs_switch)
  
  ]
)
