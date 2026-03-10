FC = mpif90

subroutines = name:calc_tau_straight,calc_tau_cross,flux4

OUTPUT_OPTIN = -cuda -acc -fast -Mchkptr -Minline=reshape,$(subroutines) -Minfo=accel,inline -Mchkptr -gpu=cc89,ptxinfo,rdc,lto -lnvhpcwrapnvtx 

.SUFFIXES : .f90

%.o: %.f90
		$(COMPILE.f) $(OUTPUT_OPTIN) $<

%.mod: %.f90 %.o
		@:

TARGET = TGV
vpath %f90 ../../src

OBJ =	mod_globals.o \
      mod_constant.o \
      cpu_gpu_mpi.o \
      set_coordinate.o \
      calc_physical_quantities.o \
      calc_keep_kernel.o \
      calc_flux_base.o \
      calc_steps.o \
      set_bc_common.o \
      set.o \
      print.o \
      preprocess.o \
      calc_time_dev.o \
      main.o

a.out: $(OBJ)
		$(FC) -cuda -acc -fast -Minfo -gpu=cc89,rdc,lto -lnvhpcwrapnvtx -I./ $^

$(TARGET): $(OBJ)
		$(LINK.f) $^ $(LOADLIBES) $(LDLIBS) -o $@

mod_constant.o: mod_globals.mod
set_bc_common.o: mod_globals.mod
set.o: mod_globals.mod mod_constant.mod set_bc_common.mod set_coordinate.mod
calc_physical_quantities.o: mod_globals.mod mod_constant.mod
calc_keep_kernel.o: mod_globals.mod calc_keep_3d.f90
calc_flux_base.o: mod_globals.mod calc_physical_quantities.mod set.mod calc_keep_kernel.mod 
calc_steps.o: mod_globals.mod
preprocess.o: mod_globals.mod print.mod
calc_time_dev.o: mod_globals.mod preprocess.mod set.mod calc_flux_base.mod calc_steps.mod print.mod
print.o: mod_globals.mod
main.o: main.f90 mod_globals.mod set_coordinate.mod set.mod calc_time_dev.mod

clean:
		rm -f *.mod *.o *.out

