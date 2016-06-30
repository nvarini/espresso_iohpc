












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



module iotk_xtox_interf
use iotk_base
implicit none
private

public :: iotk_atol
public :: iotk_ltoa
public :: iotk_atoi
public :: iotk_itoa

interface iotk_atol
function iotk_atol_x(a,check)
  character(len=*),           intent(in)  :: a
  logical,          optional, intent(out) :: check
  logical                                 :: iotk_atol_x
end function iotk_atol_x
end interface

interface iotk_atoi
subroutine iotk_atoi1(i,a,check)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_integer1
  integer(kind=this_kind),           intent(out) :: i
  character(len=*),                  intent(in)  :: a
  logical,                 optional, intent(out) :: check
end subroutine iotk_atoi1
end interface

interface iotk_itoa
function iotk_itoa1(i,length)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_integer1
  integer(kind=this_kind),             intent(in)  :: i
  integer,                   optional, intent(out) :: length
  character(len=range(i)+2)                        :: iotk_itoa1
end function iotk_itoa1
end interface

interface iotk_ltoa
function iotk_ltoa1(l)
  use iotk_base
  implicit none
  integer, parameter :: this_kind = iotk_logical1
  logical(kind=this_kind), intent(in) :: l
  character                           :: iotk_ltoa1
end function iotk_ltoa1
end interface

end module iotk_xtox_interf
