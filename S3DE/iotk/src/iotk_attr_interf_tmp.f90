












! Input/Output Tool Kit (IOTK)
! Copyright (C) 2004-2006 Giovanni Bussi
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

!------------------------------------------------------------------------------!
! Inclusion of configuration file
! Input/Output Tool Kit (IOTK)
! Copyright (C) 2004 Giovanni Bussi
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!------------------------------------------------------------------------------!
! CONFIGURATION FILE FOR IOTK 1.1.0 for Quantum-Espresso
!------------------------------------------------------------------------------!
! The following lines map some commonly defined system macro to the internal
! iotk macros.
! Iotk macros which are not defined take their default values.
! See the manual for a list of iotk macros.


! Generic options valid for quantum-espresso
! QE uses ranks up to four and default integer/logicals only


! some compilers do not like the following
!    #define __IOTK_REAL1 selected_real_kind(6,30)
!    #define __IOTK_REAL2 selected_real_kind(14,200)
! so we use explicit kinds

! Something from an ancient past that probably is not valid anymore...
! #if defined(__NAG)
! #   define __IOTK_REAL1 1
! #   define __IOTK_REAL2 2
! #else

! Machine-dependent options
! Only for compilers that require some special tricks

    !
    ! proceed with a machine dependent def where available
    !



!------------------------------------------------------------------------------!



! The macros are defined with -D option or inside iotk_config.h
! The default values are set here

! Maximum rank of an array

! Minimum value used in iotk_free_unit

! Maximum value used in iotk_free_unit

! Unit for errors

! Unit for output

! Kind for header in binary files

! Maximum number of arguments to the iotk tool

! Character (or eventually string) for newline
! It may be adjusted for particular systems
! It is used only in binary files, surrounding the tags so that they can
! be easily isolated with grep
! Unix    achar(10)
! Mac-OS  achar(13)
! Windows ? (now it should be a single byte)

! Character for EOS

! These are the default kinds, which depend on the options used
! during the library compilation
! Only default characters are implemented
! For logical, integer and real types, the c precompiler
! looks for defined kinds. If no kind is found, the default
! is used as __IOTK_type1

! Complex are treated indentically to reals
! These lines map the definitions.

! Some useful check follow



module iotk_attr_interf
implicit none
private

public :: iotk_read
public :: iotk_write
public :: iotk_write_attr
public :: iotk_scan_attr


interface iotk_read
subroutine iotk_read_LOGICAL1(val,string,index,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  LOGICAL(kind=this_kind), intent(inout) :: val(:)
  character(len=*),        intent(in)    :: string
  integer,                 intent(inout) :: index
  integer,                 intent(out)   :: ierr
end subroutine iotk_read_LOGICAL1
subroutine iotk_read_INTEGER1(val,string,index,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  INTEGER(kind=this_kind), intent(inout) :: val(:)
  character(len=*),        intent(in)    :: string
  integer,                 intent(inout) :: index
  integer,                 intent(out)   :: ierr
end subroutine iotk_read_INTEGER1
subroutine iotk_read_REAL1(val,string,index,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  REAL(kind=this_kind), intent(inout) :: val(:)
  character(len=*),        intent(in)    :: string
  integer,                 intent(inout) :: index
  integer,                 intent(out)   :: ierr
end subroutine iotk_read_REAL1
subroutine iotk_read_REAL2(val,string,index,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  REAL(kind=this_kind), intent(inout) :: val(:)
  character(len=*),        intent(in)    :: string
  integer,                 intent(inout) :: index
  integer,                 intent(out)   :: ierr
end subroutine iotk_read_REAL2
subroutine iotk_read_COMPLEX1(val,string,index,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  COMPLEX(kind=this_kind), intent(inout) :: val(:)
  character(len=*),        intent(in)    :: string
  integer,                 intent(inout) :: index
  integer,                 intent(out)   :: ierr
end subroutine iotk_read_COMPLEX1
subroutine iotk_read_COMPLEX2(val,string,index,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  COMPLEX(kind=this_kind), intent(inout) :: val(:)
  character(len=*),        intent(in)    :: string
  integer,                 intent(inout) :: index
  integer,                 intent(out)   :: ierr
end subroutine iotk_read_COMPLEX2
end interface

interface iotk_write
subroutine iotk_write_LOGICAL1(val,string,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  LOGICAL(kind=this_kind), intent(in)  :: val(:)
  character(len=*)                     :: string
  character(len=*), intent(in)         :: fmt
  integer,                 intent(out) :: ierr
end subroutine iotk_write_LOGICAL1
subroutine iotk_write_INTEGER1(val,string,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  INTEGER(kind=this_kind), intent(in)  :: val(:)
  character(len=*)                     :: string
  character(len=*), intent(in)         :: fmt
  integer,                 intent(out) :: ierr
end subroutine iotk_write_INTEGER1
subroutine iotk_write_REAL1(val,string,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  REAL(kind=this_kind), intent(in)  :: val(:)
  character(len=*)                     :: string
  character(len=*), intent(in)         :: fmt
  integer,                 intent(out) :: ierr
end subroutine iotk_write_REAL1
subroutine iotk_write_REAL2(val,string,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  REAL(kind=this_kind), intent(in)  :: val(:)
  character(len=*)                     :: string
  character(len=*), intent(in)         :: fmt
  integer,                 intent(out) :: ierr
end subroutine iotk_write_REAL2
subroutine iotk_write_COMPLEX1(val,string,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  COMPLEX(kind=this_kind), intent(in)  :: val(:)
  character(len=*)                     :: string
  character(len=*), intent(in)         :: fmt
  integer,                 intent(out) :: ierr
end subroutine iotk_write_COMPLEX1
subroutine iotk_write_COMPLEX2(val,string,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  COMPLEX(kind=this_kind), intent(in)  :: val(:)
  character(len=*)                     :: string
  character(len=*), intent(in)         :: fmt
  integer,                 intent(out) :: ierr
end subroutine iotk_write_COMPLEX2
end interface

interface iotk_write_attr
subroutine iotk_write_attr_LOGICAL1_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL1_0
subroutine iotk_write_attr_LOGICAL1_1(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL1_1
subroutine iotk_write_attr_LOGICAL1_2(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL1_2
subroutine iotk_write_attr_LOGICAL1_3(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL1_3
subroutine iotk_write_attr_LOGICAL1_4(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  LOGICAL(kind=this_kind),           intent(in)  :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_LOGICAL1_4
subroutine iotk_write_attr_INTEGER1_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER1_0
subroutine iotk_write_attr_INTEGER1_1(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER1_1
subroutine iotk_write_attr_INTEGER1_2(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER1_2
subroutine iotk_write_attr_INTEGER1_3(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER1_3
subroutine iotk_write_attr_INTEGER1_4(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  INTEGER(kind=this_kind),           intent(in)  :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_INTEGER1_4
subroutine iotk_write_attr_REAL1_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL1_0
subroutine iotk_write_attr_REAL1_1(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL1_1
subroutine iotk_write_attr_REAL1_2(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL1_2
subroutine iotk_write_attr_REAL1_3(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL1_3
subroutine iotk_write_attr_REAL1_4(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL1_4
subroutine iotk_write_attr_REAL2_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL2_0
subroutine iotk_write_attr_REAL2_1(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL2_1
subroutine iotk_write_attr_REAL2_2(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL2_2
subroutine iotk_write_attr_REAL2_3(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL2_3
subroutine iotk_write_attr_REAL2_4(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  REAL(kind=this_kind),           intent(in)  :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_REAL2_4
subroutine iotk_write_attr_COMPLEX1_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX1_0
subroutine iotk_write_attr_COMPLEX1_1(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX1_1
subroutine iotk_write_attr_COMPLEX1_2(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX1_2
subroutine iotk_write_attr_COMPLEX1_3(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX1_3
subroutine iotk_write_attr_COMPLEX1_4(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX1_4
subroutine iotk_write_attr_COMPLEX2_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX2_0
subroutine iotk_write_attr_COMPLEX2_1(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX2_1
subroutine iotk_write_attr_COMPLEX2_2(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX2_2
subroutine iotk_write_attr_COMPLEX2_3(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX2_3
subroutine iotk_write_attr_COMPLEX2_4(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  COMPLEX(kind=this_kind),           intent(in)  :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_COMPLEX2_4
subroutine iotk_write_attr_CHARACTER1_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  character(len=*), intent(inout) :: attr
  character(len=*), intent(in)    :: name
  CHARACTER(kind=this_kind,len=*),           intent(in)  :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(in)  :: first
  logical,                         optional, intent(in)  :: newline
  character(len=*), optional, intent(in)                 :: fmt
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_write_attr_CHARACTER1_0
end interface

interface iotk_scan_attr
subroutine iotk_scan_attr_LOGICAL1_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind)                        :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default 
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL1_0
subroutine iotk_scan_attr_LOGICAL1_1(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind)                        :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL1_1
subroutine iotk_scan_attr_LOGICAL1_2(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind)                        :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL1_2
subroutine iotk_scan_attr_LOGICAL1_3(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind)                        :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL1_3
subroutine iotk_scan_attr_LOGICAL1_4(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_LOGICAL1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  LOGICAL(kind=this_kind)                        :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  LOGICAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_LOGICAL1_4
subroutine iotk_scan_attr_INTEGER1_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind)                        :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default 
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER1_0
subroutine iotk_scan_attr_INTEGER1_1(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind)                        :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER1_1
subroutine iotk_scan_attr_INTEGER1_2(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind)                        :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER1_2
subroutine iotk_scan_attr_INTEGER1_3(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind)                        :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER1_3
subroutine iotk_scan_attr_INTEGER1_4(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_INTEGER1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  INTEGER(kind=this_kind)                        :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  INTEGER(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_INTEGER1_4
subroutine iotk_scan_attr_REAL1_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind)                        :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default 
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL1_0
subroutine iotk_scan_attr_REAL1_1(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind)                        :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL1_1
subroutine iotk_scan_attr_REAL1_2(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind)                        :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL1_2
subroutine iotk_scan_attr_REAL1_3(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind)                        :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL1_3
subroutine iotk_scan_attr_REAL1_4(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind)                        :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL1_4
subroutine iotk_scan_attr_REAL2_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind)                        :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default 
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL2_0
subroutine iotk_scan_attr_REAL2_1(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind)                        :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL2_1
subroutine iotk_scan_attr_REAL2_2(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind)                        :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL2_2
subroutine iotk_scan_attr_REAL2_3(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind)                        :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL2_3
subroutine iotk_scan_attr_REAL2_4(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_REAL2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  REAL(kind=this_kind)                        :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  REAL(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_REAL2_4
subroutine iotk_scan_attr_COMPLEX1_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind)                        :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default 
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX1_0
subroutine iotk_scan_attr_COMPLEX1_1(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind)                        :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX1_1
subroutine iotk_scan_attr_COMPLEX1_2(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind)                        :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX1_2
subroutine iotk_scan_attr_COMPLEX1_3(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind)                        :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX1_3
subroutine iotk_scan_attr_COMPLEX1_4(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind)                        :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX1_4
subroutine iotk_scan_attr_COMPLEX2_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind)                        :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default 
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX2_0
subroutine iotk_scan_attr_COMPLEX2_1(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind)                        :: val (:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX2_1
subroutine iotk_scan_attr_COMPLEX2_2(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind)                        :: val (:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX2_2
subroutine iotk_scan_attr_COMPLEX2_3(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind)                        :: val (:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX2_3
subroutine iotk_scan_attr_COMPLEX2_4(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_COMPLEX2
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  COMPLEX(kind=this_kind)                        :: val (:,:,:,:)
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  COMPLEX(kind=this_kind), optional, intent(in)  :: default (:,:,:,:)
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_COMPLEX2_4
subroutine iotk_scan_attr_CHARACTER1_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_CHARACTER1
  character(len=*),                          intent(in)  :: attr
  character(len=*),                          intent(in)  :: name
  CHARACTER(kind=this_kind,len=*)                        :: val 
  type(iotk_dummytype),            optional              :: dummy
  logical,                         optional, intent(out) :: found
  CHARACTER(kind=this_kind,len=*), optional, intent(in)  :: default 
  logical,                         optional, intent(in)  :: eos
  integer,                         optional, intent(out) :: ierr
end subroutine iotk_scan_attr_CHARACTER1_0
end interface

end module iotk_attr_interf
