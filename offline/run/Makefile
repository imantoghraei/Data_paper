obj_launch = acama_gwd_rando.o flott_gwd_rando.o orografi_strato.o laun_gwd_healpix.o

all: laun_gwd_healpix

laun_gwd_healpix: ${obj_launch}
	$(FC) $(LDFLAGS) $^ $(LDLIBS) -o $@


LDLIBS = -lnetcdff -lnetcdf -lhealpix -lcfitsio -L/home/itoghraei/Healpix_3.82/lib  -L/home/itoghraei/cfitsio-4.3.0/lib -Wl,-rpath,/home/itoghraei/cfitsio-4.3.0/lib
#FFLAGS = -I/usr/include -mcmodel=large 
FFLAGS = -I/lib64/gfortran/modules -mcmodel=large -I/home/itoghraei/Healpix_3.82/include 
FC = gfortran

%.o: %.f90
	$(COMPILE.f90) $(OUTPUT_OPTION) $<
COMPILE.f90 = $(FC) $(F90FLAGS) $(TARGET_ARCH) -c

F90FLAGS= -ffree-form ${FFLAGS}

laun_gwd_healpix.o: acama_gwd_rando.o flott_gwd_rando.o

clean:
	\rm  *.h *.f90  *.nc *.o *.mod input_data laun_gwd_healpix data *ctl ../netsto/gwd_healpix_temp.nc


