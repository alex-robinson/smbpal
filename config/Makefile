.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
srcdir = src
objdir = libsmbpal/include
bindir = libsmbpal/bin
libdir = libs

# Command-line options at make call
debug  ?= 0
openmp ?= 0 

## COMPILER CONFIGURATION ##
# (should be loaded from config directory)

<COMPILER_CONFIGURATION>

# Determine whether to use normal flags or debugging flags
DFLAGS   = $(DFLAGS_NODEBUG)
ifeq ($(debug), 1)
	DFLAGS   = $(DFLAGS_DEBUG)
endif

# Debugging flags with profiling output enabled
ifeq ($(debug), 2)
	DFLAGS   = $(DFLAGS_PROFILE)
endif

###############################################
##
## List of rules and source files
##
###############################################

## Individual libraries or modules ##
$(objdir)/nml.o: $(libdir)/nml.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/ncio.o: $(libdir)/ncio.f90
	$(FC) $(DFLAGS) $(FFLAGS) $(INC_NC) -c -o $@ $<

$(objdir)/interp1D.o: $(libdir)/insol/interp1D.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/interp_time.o: $(libdir)/interp_time.f90 $(objdir)/interp1D.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/insolation.o: $(libdir)/insol/insolation.f90 $(objdir)/interp_time.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/smbpal_precision.o: $(srcdir)/smbpal_precision.f90
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/smb_itm.o: $(srcdir)/smb_itm.f90 $(objdir)/smbpal_precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/smb_pdd.o: $(srcdir)/smb_pdd.f90 $(objdir)/smbpal_precision.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

$(objdir)/smbpal.o: $(srcdir)/smbpal.f90 $(objdir)/smbpal_precision.o $(objdir)/nml.o $(objdir)/insolation.o  \
					$(objdir)/interp1D.o $(objdir)/interp_time.o $(objdir)/ncio.o \
					$(objdir)/smb_pdd.o $(objdir)/smb_itm.o
	$(FC) $(DFLAGS) $(FFLAGS) -c -o $@ $<

# Group the object files together
obj_smbpal =     $(objdir)/insolation.o \
				 $(objdir)/interp1D.o \
				 $(objdir)/interp_time.o \
				 $(objdir)/nml.o \
				 $(objdir)/ncio.o \
				 $(objdir)/smb_itm.o \
				 $(objdir)/smb_pdd.o \
				 $(objdir)/smbpal.o \
				 $(objdir)/smbpal_precision.o

###############################################
##
## Compilation of complete programs
##
###############################################

$(objdir)/libsmbpal.a: $(obj_smbpal)
	ar rcs $@ $^
	@echo " "
	@echo "    $(objdir)/libsmbpal.a is ready."
	@echo " "

.PHONY: lib
lib: $(objdir)/libsmbpal.a

test: $(objdir)/libsmbpal.a
	$(FC) $(DFLAGS) $(FFLAGS) -o $(bindir)/test_smbpal.x test/test_smbpal.f90 \
		-L$(objdir) -lsmbpal -I$(objdir) $(LFLAGS)
	@echo " "
	@echo "    $(bindir)/test_smbpal.x is ready."
	@echo " "

clean:
	rm -f $(bindir)/*.x
	rm -f $(libdir)/libsmbpal.a
	rm -f  *.x gmon.out $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.so
	rm -rf *.x.dSYM