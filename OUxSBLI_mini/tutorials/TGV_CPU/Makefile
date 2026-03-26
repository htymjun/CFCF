FC = mpif90

OUTPUT_OPTIN = -O3 -fopenmp -ffree-line-length-none

.SUFFIXES : .f90

%.o: %.f90
		$(COMPILE.f) $(OUTPUT_OPTIN) $<

%.mod: %.f90 %.o
		@:

TARGET = TGV
vpath %f90 ../../src_cpu

OBJ =	mod_globals.o \
      mod_constant.o \
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
		$(FC) -O3 -fopenmp -I./ $^

$(TARGET): $(OBJ)
		$(LINK.f) $^ $(LOADLIBES) $(LDLIBS) -o $@

mod_constant.o: mod_globals.mod
set_bc_common.o: mod_globals.mod
set.o: mod_globals.mod mod_constant.mod set_bc_common.mod set_coordinate.mod
calc_physical_quantities.o: mod_globals.mod mod_constant.mod
calc_keep_kernel.o: mod_globals.mod calc_keep_3d.f90
calc_flux_base.o: mod_globals.mod calc_physical_quantities.mod calc_keep_kernel.mod
calc_steps.o: mod_globals.mod
preprocess.o: mod_globals.mod print.mod
calc_time_dev.o: mod_globals.mod preprocess.mod set.mod calc_flux_base.mod calc_steps.mod print.mod
print.o: mod_globals.mod
main.o: main.f90 mod_globals.mod set_coordinate.mod set.mod calc_time_dev.mod

clean:
		rm -f *.mod *.o *.out
