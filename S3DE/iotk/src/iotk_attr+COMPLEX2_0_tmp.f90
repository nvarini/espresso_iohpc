












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
subroutine iotk_private_pack_COMPLEX2(out,in,n,l)
    use iotk_base
    implicit none
    integer,                                    intent(in)  :: n,l
    COMPLEX (kind=iotk_COMPLEX2), intent(out) :: out(n)
    COMPLEX (kind=iotk_COMPLEX2), intent(in)  :: in(n)
    out = in
end subroutine iotk_private_pack_COMPLEX2

subroutine iotk_write_COMPLEX2(val,string,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_xtox_interf
  use iotk_fmt_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  COMPLEX(kind=iotk_COMPLEX2), intent(in) :: val(:)
  character(len=*)              :: string
  character(len=*), intent(in)  :: fmt
  integer, intent(out) :: ierr
  character(len=100) :: tmpval
  integer :: index,iostat
  ierr = 0
  iostat = 0 
  string(1:1) = iotk_eos
  if(size(val)==0) return
  if(len(string)==0) then
    call iotk_error_issue(ierr,"iotk_write","iotk_attr+COMPLEX2_0.f90",63)
call iotk_error_msg(ierr,"CVS Revision: 1.21 ")
    return
  end if
  do index=1,size(val)
    if (trim(fmt)=="!") then
       write(tmpval,trim(iotk_wfmt("COMPLEX",kind(val),1,-1," ")),iostat=iostat) val(index)
    else
       write(tmpval,fmt=fmt,iostat=iostat) val(index)
    end if
    if(iostat/=0) then
      call iotk_error_issue(ierr,"iotk_write","iotk_attr+COMPLEX2_0.f90",74)
call iotk_error_msg(ierr,"CVS Revision: 1.21 ")
call iotk_error_msg(ierr,' ')
call iotk_error_write(ierr,"iostat",iostat)
      return
    end if
    call iotk_strcat(string,trim(adjustl(tmpval))//" ",ierr)
    if(ierr/=0) then
      call iotk_error_issue(ierr,"iotk_write","iotk_attr+COMPLEX2_0.f90",82)
call iotk_error_msg(ierr,"CVS Revision: 1.21 ")
      return
    end if
  end do
! the last blank is deleted
  string(iotk_strlen(string):iotk_strlen(string)) = iotk_eos
end subroutine iotk_write_COMPLEX2

subroutine iotk_read_COMPLEX2(val,string,index,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_xtox_interf
  use iotk_misc_interf
  implicit none
  COMPLEX(kind=iotk_COMPLEX2), intent(inout) :: val(:)
  character(len=*), intent(in) :: string
  integer, intent(inout) :: index
  integer, intent(out) :: ierr
  integer :: pos,pos1,iostat
  integer :: maxindex
  real(kind=iotk_COMPLEX2) :: tmpreal
  complex(kind=iotk_COMPLEX2) :: tmpcomplex
  pos = 0
  pos1= 0
  ierr = 0
  iostat = 0
   maxindex = 2 * size(val)
! for the moment, commas are considered as blanks
  do
    pos = verify(string(pos1+1:)," ,")+pos1
    if(pos==pos1) exit
    pos = pos - 1
    pos1 = scan(string(pos+1:)," ,")+pos
    if(pos1==pos) pos1 = len(string) + 1
!READ string(pos+1:pos1-1)
    index = index+1
    if(index>maxindex) then
      call iotk_error_issue(ierr,"iotk_read","iotk_attr+COMPLEX2_0.f90",123)
call iotk_error_msg(ierr,"CVS Revision: 1.21 ")
call iotk_error_msg(ierr,'Too many data')
    end if
    read(string(pos+1:pos1-1),"(G100.95)",iostat=iostat) tmpreal
    if(modulo(index,2)==1) then
      tmpcomplex = cmplx(tmpreal,aimag((val((index+1)/2))),kind=iotk_COMPLEX2)
    else
      tmpcomplex = cmplx(real(val((index+1)/2)),tmpreal,kind=iotk_COMPLEX2)
    end if
    val((index+1)/2) = tmpcomplex
    if(iostat/=0) then
      call iotk_error_issue(ierr,"iotk_read","iotk_attr+COMPLEX2_0.f90",140)
call iotk_error_msg(ierr,"CVS Revision: 1.21 ")
call iotk_error_msg(ierr,'Error reading a COMPLEX number from string')
call iotk_error_write(ierr,"string",string(pos+1:pos1-1))
call iotk_error_write(ierr,"iostat",iostat)
      return
    end if
    if(pos1>=len(string)) exit
  end do
end subroutine iotk_read_COMPLEX2

subroutine iotk_write_attr_COMPLEX2_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=iotk_COMPLEX2), intent(in)  :: val 
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
    call iotk_error_issue(ierrl,"iotk_write_attr","iotk_attr+COMPLEX2_0.f90",193)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'Wrong tag name')
call iotk_error_write(ierrl,"name",name(1:namlen))
    goto 1
  end if
  delim = '"'
  call iotk_write((/val/),tmpval,usefmt,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr","iotk_attr+COMPLEX2_0.f90",202)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
    goto 1
  end if
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr","iotk_attr+COMPLEX2_0.f90",208)
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
end subroutine iotk_write_attr_COMPLEX2_0

subroutine iotk_scan_attr_COMPLEX2_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=iotk_COMPLEX2)                        :: val 
  type(iotk_dummytype), optional :: dummy
  logical,        optional, intent(out) :: found
  COMPLEX(kind=iotk_COMPLEX2), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal,namlen
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
  integer :: index
  COMPLEX(kind=iotk_COMPLEX2), allocatable :: tmpval (:)
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
      call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+COMPLEX2_0.f90",264)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'')
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==name(1:namlen)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+COMPLEX2_0.f90",274)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+COMPLEX2_0.f90",281)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+COMPLEX2_0.f90",287)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+COMPLEX2_0.f90",297)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
  else
    goto 1
  end if
  allocate(tmpval(1))
  index = 0
  call iotk_str_clean(valc)
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+COMPLEX2_0.f90",309)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
    goto 1
  end if
  if(index/=2*1) then
    call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+COMPLEX2_0.f90",314)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'Attribute size does not match')
call iotk_error_write(ierrl,"attr",valc)
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
  val = tmpval(1)
  deallocate(tmpval)
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+COMPLEX2_0.f90",327)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'Attribute not found')
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
    val = default
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX2_0


subroutine iotk_attr_dummy_COMPLEX2_0
  write(0,*)
end subroutine iotk_attr_dummy_COMPLEX2_0




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

!------------------------------------------------------------------------------!







subroutine iotk_write_attr_COMPLEX2_1(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=iotk_COMPLEX2), intent(in)  :: val (:)
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
    call iotk_error_issue(ierrl,"iotk_write_attr","iotk_attr+COMPLEX2_0.f90",408)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'Wrong tag name')
call iotk_error_write(ierrl,"name",name(1:namlen))
    goto 1
  end if
  delim = '"'
  call iotk_write(pack(val,mask=.true.),tmpval,usefmt,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr","iotk_attr+COMPLEX2_0.f90",417)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
    goto 1
  end if
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr","iotk_attr+COMPLEX2_0.f90",423)
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
end subroutine iotk_write_attr_COMPLEX2_1

subroutine iotk_scan_attr_COMPLEX2_1(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=iotk_COMPLEX2)                        :: val (:)
  type(iotk_dummytype), optional :: dummy
  logical,        optional, intent(out) :: found
  COMPLEX(kind=iotk_COMPLEX2), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal,namlen
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
  integer :: index
  COMPLEX(kind=iotk_COMPLEX2), allocatable :: tmpval (:)
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
      call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+COMPLEX2_0.f90",479)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'')
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==name(1:namlen)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+COMPLEX2_0.f90",489)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+COMPLEX2_0.f90",496)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+COMPLEX2_0.f90",502)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+COMPLEX2_0.f90",512)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
  else
    goto 1
  end if
  allocate(tmpval(size(val)))
  index = 0
  call iotk_str_clean(valc)
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+COMPLEX2_0.f90",524)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
    goto 1
  end if
  if(index/=2*size(val)) then
    call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+COMPLEX2_0.f90",529)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'Attribute size does not match')
call iotk_error_write(ierrl,"attr",valc)
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
  val = reshape (source=tmpval,shape=shape(val))
  deallocate(tmpval)
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+COMPLEX2_0.f90",542)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'Attribute not found')
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
    val = default
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX2_1


subroutine iotk_attr_dummy_COMPLEX2_1
  write(0,*)
end subroutine iotk_attr_dummy_COMPLEX2_1




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

!------------------------------------------------------------------------------!







subroutine iotk_write_attr_COMPLEX2_2(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=iotk_COMPLEX2), intent(in)  :: val (:,:)
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
    call iotk_error_issue(ierrl,"iotk_write_attr","iotk_attr+COMPLEX2_0.f90",623)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'Wrong tag name')
call iotk_error_write(ierrl,"name",name(1:namlen))
    goto 1
  end if
  delim = '"'
  call iotk_write(pack(val,mask=.true.),tmpval,usefmt,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr","iotk_attr+COMPLEX2_0.f90",632)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
    goto 1
  end if
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr","iotk_attr+COMPLEX2_0.f90",638)
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
end subroutine iotk_write_attr_COMPLEX2_2

subroutine iotk_scan_attr_COMPLEX2_2(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
  COMPLEX(kind=iotk_COMPLEX2)                        :: val (:,:)
  type(iotk_dummytype), optional :: dummy
  logical,        optional, intent(out) :: found
  COMPLEX(kind=iotk_COMPLEX2), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal,namlen
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
  integer :: index
  COMPLEX(kind=iotk_COMPLEX2), allocatable :: tmpval (:)
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
      call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+COMPLEX2_0.f90",694)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'')
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==name(1:namlen)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+COMPLEX2_0.f90",704)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+COMPLEX2_0.f90",711)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+COMPLEX2_0.f90",717)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+COMPLEX2_0.f90",727)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
  else
    goto 1
  end if
  allocate(tmpval(size(val)))
  index = 0
  call iotk_str_clean(valc)
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+COMPLEX2_0.f90",739)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
    goto 1
  end if
  if(index/=2*size(val)) then
    call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+COMPLEX2_0.f90",744)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'Attribute size does not match')
call iotk_error_write(ierrl,"attr",valc)
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
  val = reshape (source=tmpval,shape=shape(val))
  deallocate(tmpval)
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr","iotk_attr+COMPLEX2_0.f90",757)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'Attribute not found')
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
    val = default
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX2_2


subroutine iotk_attr_dummy_COMPLEX2_2
  write(0,*)
end subroutine iotk_attr_dummy_COMPLEX2_2


