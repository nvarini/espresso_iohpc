# Makefile for IO_kernel

include ../../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../Modules \
          $(MOD_FLAG)../Tools \
          $(MOD_FLAG)../FFT_kernel \
          $(MOD_FLAG)../Para_kernel \
          $(MOD_FLAG).
IFLAGS=

LIBOBJS        = ../../flib/ptools.a ../../flib/flib.a ../../clib/clib.a ../../iotk/src/libiotk.a

IO_KERNEL_OBJS = \
mpiio.o \
my_mkdir.o \
fetch_namelist.o \
pdep_io.o \
pdep_db.o \
mod_west_io.o \
wfreq_io.o \
wfreq_restart.o  

PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a

TLDEPS= bindir mods libs pw

all : title tldeps lib_io_kernel.a

tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

lib_io_kernel.a : $(IO_KERNEL_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.d *.i *~ *.F90 *.mod *.L
	- /bin/rm -f lib_io_kernel.a

title :
	@echo "                  "
	@echo "##################"
	@echo "### IO_kernel ####"
	@echo "##################"
	@echo "                  "

include make.depend
