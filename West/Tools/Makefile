# Makefile for Tools

include ../../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../Modules \
          $(MOD_FLAG).
IFLAGS=

LIBOBJS        = ../../flib/ptools.a ../../flib/flib.a ../../clib/clib.a ../../iotk/src/libiotk.a

TOOLS_OBJS = \
exx_go.o \
exx_ungo.o \
report_dynamical_memory.o \
io_push.o \
do_setup.o \
linear_algebra_kernel.o \
human_readable_time.o \
destroy_pw_arrays.o \
init_pw_arrays.o \
get_alpha_pv.o \
clean_scratchfiles.o \
set_nbndocc.o \
set_iks_l2g.o \
set_isz.o \
west_environment.o \
bar.o  

PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a

TLDEPS= bindir mods libs pw

all : title tldeps libtools.a

tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

libtools.a : $(TOOLS_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.d *.i *~ *.F90 *.mod *.L
	- /bin/rm -f libtools.a

title :
	@echo "              "
	@echo "##############"
	@echo "### Tools ####"
	@echo "##############"
	@echo "              "

include make.depend
