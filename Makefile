.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
objdir = .obj
libdir = libs
srcdir = .

# Command-line options at make call
env   ?= None      # options: manto,eolo,airaki,iplex
debug ?= 0 

ifeq ($(env),manto) ## env=manto

    ## IFORT OPTIONS ##
    FC  = ifort
    INC_NC  = -I/home/jalvarez/work/librairies/netcdflib/include
    LIB_NC  = -L/home/jalvarez/work/librairies/netcdflib/lib -lnetcdf
    LIB_MKL = -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

    FLAGS    = -module $(objdir) -L$(objdir) $(INC_NC)
    LFLAGS   = $(LIB_NC) $(LIB_MKL)

    DFLAGS   = -vec-report0 -O2 -fp-model precise -i_dynamic 
    ifeq ($(debug), 1)
        DFLAGS   = -C -traceback -ftrapuv -fpe0 -check all -vec-report0 -fp-model precise -i_dynamic 
    endif

else ifeq ($(env),eolo) ## env=eolo

#    ## IFORT OPTIONS ##
#    FC  = ifort
#    INC_NC  = -I/home/fispalma22/work/librairies/netcdflib/include
#    LIB_NC  = -L/home/fispalma22/work/librairies/netcdflib/lib -lnetcdf
#    INC_COORD = -I/home/fispalma25/robinson/models/EURICE/coord/.obj
#    LIB_COORD = /home/fispalma25/robinson/models/EURICE/coord/libcoordinates.a
#    LIB_MKL = -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
#
#    FLAGS    = -module $(objdir) -L$(objdir) $(INC_COORD) $(INC_NC)
#    LFLAGS   = $(LIB_COORD) $(LIB_NC) $(LIB_MKL)
#
#    DFLAGS   = -vec-report0 -O2 -fp-model precise
#    ifeq ($(debug), 1)
#        DFLAGS   = -C -traceback -ftrapuv -fpe0 -check all -vec-report0 -fp-model precise
#    endif

    ## GFORTRAN OPTIONS ##
    FC  = gfortran
    INC_NC  = -I/home/fispalma25/apps/netcdf/netcdf/include
    LIB_NC  = -L/home/fispalma25/apps/netcdf/netcdf/lib -lnetcdff -lnetcdf
    INC_COORD = -I/home/fispalma25/apps/coordinates/.obj
    LIB_COORD = /home/fispalma25/apps/coordinates/libcoordinates.a

    FLAGS  = -I$(objdir) -J$(objdir) $(INC_COORD) $(INC_NC)
    LFLAGS = $(LIB_COORD) $(LIB_NC)

    DFLAGS = -O3
    ifeq ($(debug), 1)  # ,underflow
        DFLAGS   = -w -g -p -ggdb -ffpe-trap=invalid,zero,overflow -fbacktrace -fcheck=all
    endif

else ifeq ($(env),airaki) ## env=airaki

    ## GFORTRAN OPTIONS ##
    FC  = gfortran
    INC_NC  = -I/opt/local/include
    LIB_NC  = -L/opt/local/lib -lnetcdff -lnetcdf
    INC_COORD = -I/Users/robinson/models/EURICE/coord/.obj
    LIB_COORD = /Users/robinson/models/EURICE/coord/libcoordinates.a

    FLAGS  = -I$(objdir) -J$(objdir) $(INC_COORD) $(INC_NC) 
    LFLAGS = $(LIB_COORD) $(LIB_NC)

    DFLAGS = -O3
    ifeq ($(debug), 1)  # ,underflow
        DFLAGS   = -w -g -p -ggdb -ffpe-trap=invalid,zero,overflow -fbacktrace -fcheck=all
    endif

else ifeq ($(env),pik) ## env=pik

    ## IFORT OPTIONS ##
    FC  = ifort
    INC_NC  = -I${NETCDF_FORTRANROOT}/include
    LIB_NC  = -L${NETCDF_FORTRANROOT}/lib -lnetcdff -L${NETCDF_CROOT}/lib -lnetcdf
    LIB_MKL = -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
    INC_COORD = -I/p/projects/tumble/robinson/EURICE/coord/.obj
	LIB_COORD = /p/projects/tumble/robinson/EURICE/coord/libcoordinates.a

    FLAGS    = -module $(objdir) -L$(objdir) $(INC_COORD) $(INC_NC) 
    LFLAGS   = $(LIB_COORD) $(LIB_NC)

    DFLAGS   = -O3
    ifeq ($(debug), 1)
        DFLAGS   = -C -g -traceback -ftrapuv -fpe0 -check all
    endif

else 
    
    ## None ##
    FC = $(error "Define env")

endif

## Individual libraries or modules ##
$(objdir)/nml.o: $(libdir)/nml.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/ncio.o: $(libdir)/ncio.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/interp1D.o: $(libdir)/insol/interp1D.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/interp_time.o: $(libdir)/coordinates/interp_time.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/insolation.o: $(libdir)/insol/insolation.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/smbpal_precision.o: $(srcdir)/smbpal_precision.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/smb_itm.o: $(srcdir)/smb_itm.f90 $(objdir)/smbpal_precision.o
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/smb_pdd.o: $(srcdir)/smb_pdd.f90 $(objdir)/smbpal_precision.o
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/smbpal.o: $(srcdir)/smbpal.f90 $(objdir)/smbpal_precision.o $(objdir)/nml.o $(objdir)/insolation.o  \
					$(objdir)/interp1D.o $(objdir)/interp_time.o $(objdir)/ncio.o \
					$(objdir)/smb_pdd.o $(objdir)/smb_itm.o
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

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

## Complete programs

test: $(obj_smbpal) 
	$(FC) $(DFLAGS) $(FLAGS) -o test_smbpal.x $^ test_smbpal.f90 $(LFLAGS)
	@echo " "
	@echo "    test_smbpal.x is ready."
	@echo " "

clean:
	rm -f *.x $(objdir)/*.o $(objdir)/*.mod
