# Makefile for Modules

include ../../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules $(MOD_FLAG)../../LAXlib \
          $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG).
IFLAGS=

LIBOBJS        = ../../clib/clib.a ../../iotk/src/libiotk.a

MODULES_OBJS = \
west_version.o \
westcom.o

PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a

TLDEPS= bindir mods libs pw

all : title tldeps libmodules.a

tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

libmodules.a : $(MODULES_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.d *.i *~ *.F90 *.mod *.L 
	- /bin/rm -f libmodules.a
	- /bin/rm -f west_version.f90

title :
	@echo "                "
	@echo "################"
	@echo "### Modules ####"
	@echo "################"
	@echo "                "

include make.depend
