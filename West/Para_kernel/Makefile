# Makefile for Para_kernel

include ../../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules $(MOD_FLAG)../../LAXlib \
          $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../Modules \
          $(MOD_FLAG)../Tools \
          $(MOD_FLAG).
IFLAGS=

LIBOBJS        = ../../clib/clib.a ../../iotk/src/libiotk.a

PARA_KERNEL_OBJS = \
distribution_center.o \
parallel_distributed_diago.o  

PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a

TLDEPS= bindir mods libs pw

all : title tldeps lib_para_kernel.a

tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

lib_para_kernel.a : $(PARA_KERNEL_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.d *.i *~ *.F90 *.mod *.L
	- /bin/rm -f lib_para_kernel.a

title :
	@echo "                    "
	@echo "####################"
	@echo "### Para_kernel ####"
	@echo "####################"
	@echo "                    "

include make.depend
