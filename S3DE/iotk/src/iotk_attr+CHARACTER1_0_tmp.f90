












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





! This is needed as a workaround for bugged pack 
subroutine iotk_private_pack_CHARACTER1(out,in,n,l)
    use iotk_base
    implicit none
    integer,                                    intent(in)  :: n,l
    CHARACTER (kind=iotk_CHARACTER1,len=l), intent(out) :: out(n)
    CHARACTER (kind=iotk_CHARACTER1,len=l), intent(in)  :: in(n)
    out = in
end subroutine iotk_private_pack_CHARACTER1



subroutine iotk_write_attr_CHARACTER1_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  CHARACTER(kind=iotk_CHARACTER1,len=*), intent(in)  :: val 
  type(iotk_dummytype), optional :: dummy
  logical, optional, intent(in)  :: first
  logical, optional, intent(in)  :: newline
  character(*), optional, intent(in) :: fmt
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  character :: delim
  character(len=300) :: usefmt
  logical :: lquot,lapos
  character(iotk_vallenx) :: tmpval
  logical :: nl
  if(present(newline)) then
    nl = newline
  else
    nl = .false.
  endif
!-<
  if (present(fmt)) then
    usefmt = fmt
  else
    usefmt = "!"
  end if
!->
  ierrl = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  attlen = iotk_strlen_trim(attr)
  namlen = iotk_strlen_trim(name)
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr","iotk_attr+CHARACTER1_0.f90",85)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'Wrong tag name')
call iotk_error_write(ierrl,"name",name(1:namlen))
    goto 1
  end if
  lquot=iotk_strscan(val,'"')>0
  lapos=iotk_strscan(val,"'")>0
  if(.not.lquot) then
    delim='"'
    call iotk_deescape(tmpval,val)
  else if(.not.lapos) then
    delim="'"
    call iotk_deescape(tmpval,val)
  else
    delim='"'
    call iotk_deescape(tmpval,val,quot=.true.,apos=.true.)
  end if
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr","iotk_attr+CHARACTER1_0.f90",105)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'Attribute dummy argument is too short')
    goto 1
  end if
  if(.not. nl) then
    attr(attlen+1:attlen+vallen+namlen+5) = " "//name(1:namlen)//"="//delim//tmpval(1:vallen)//delim//iotk_eos
  else
    attr(attlen+1:attlen+vallen+namlen+len(iotk_newline)+5) &
       = iotk_newline//" "//name(1:namlen)//"="//delim//tmpval(1:vallen)//delim//iotk_eos
  endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_attr_CHARACTER1_0

subroutine iotk_scan_attr_CHARACTER1_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  CHARACTER(kind=iotk_CHARACTER1,len=*)                        :: val 
  type(iotk_dummytype), optional :: dummy
  logical,        optional, intent(out) :: found
  CHARACTER(kind=iotk_CHARACTER1,len=*), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal,namlen
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
  character(iotk_vallenx) :: valctmp
  integer :: vallen,defaultlen
  logical :: leos
  leos=.false.
  if(present(eos)) leos=eos
  ierrl = 0
  attlen=iotk_strlen(attr)
  namlen=iotk_strlen_trim(name)
  foundl = .false.
  equal = 0
  do
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) exit
    equal = equal + pos
    pos = scan(attr(equal+1:attlen),"=")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+CHARACTER1_0.f90",164)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'')
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==name(1:namlen)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+CHARACTER1_0.f90",174)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+CHARACTER1_0.f90",181)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+CHARACTER1_0.f90",187)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+CHARACTER1_0.f90",197)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
  else
    goto 1
  end if
  call iotk_escape(valctmp,valc)
  vallen = iotk_strlen(valctmp)
  if(len(val) < vallen) then
    call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+CHARACTER1_0.f90",207)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
    goto 1
  end if
  val(1:vallen) = valctmp(1:vallen)
  if(len(val) > vallen) then
    val(vallen+1:vallen+1) = iotk_eos
    if(.not.leos) then
      val(vallen+1:)=" "
    end if
  end if
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+CHARACTER1_0.f90",222)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'Attribute not found')
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
    if(leos) then
      defaultlen = min(iotk_strlen(default),len(val))
      val(1:defaultlen) = default(1:defaultlen)
      if(defaultlen<len(val)) val(defaultlen+1:defaultlen+1)=iotk_eos
    else
      val = default
    end if
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_CHARACTER1_0


subroutine iotk_attr_dummy_CHARACTER1_0
  write(0,*)
end subroutine iotk_attr_dummy_CHARACTER1_0

