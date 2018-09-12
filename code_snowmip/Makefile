.SUFFIXES:
.SUFFIXES: .o .F

include user_build_options


OBJS =	module_sf_ruclsm.o \
	module_model_constants.o \
	kwm_date_utilities.o \
	module_netcdf_io.o \
	module_io.o \
	module_ascii_io.o \
	module_sfcdif_wrf.o \
	module_Noahlsm_utility.o

CMD = driver.exe 

FFLAGS 	= -DEM_CORE=1 $(FREESOURCE) $(F90FLAGS)

all:	$(CMD)

driver.exe:	$(OBJS)
	$(COMPILERF90) -o $(@) $(OBJS) simple_driver.o $(NETCDFLIB)

.F.o:
	$(RM) $(*).f
	$(CPP) $(CPPMACROS) $(NETCDFINC) $(*).F > $(*).f
	$(COMPILERF90) -c $(FFLAGS) $(NETCDFINC) $(*).f
	$(RM) $(*).f

clean:
	$(RM) $(OBJS) simple_driver.o $(CMD) *.mod *~

driver.exe:		simple_driver.o
module_sf_ruclsm.o:	module_model_constants.o
simple_driver.o:	module_io.o
simple_driver.o:	kwm_date_utilities.o
simple_driver.o:	module_Noahlsm_utility.o
module_io.o:		module_ascii_io.o
module_io.o:		module_netcdf_io.o
