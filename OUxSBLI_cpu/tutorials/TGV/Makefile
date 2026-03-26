FC = gfortran
FFLAGS = -O3 -fdefault-real-8 -fdefault-double-8 -Wall -cpp

.SUFFIXES : .f90

%.o: %.f90
	$(FC) -c $(FFLAGS) $<

%.mod: %.f90 %.o
	@:

TARGET = a_cpu.out
vpath %f90 ../../src

OBJ = mod_globals.o \
      mod_constant.o \
      set_coordinate.o \
      calc_physical_quantities.o \
      calc_keep_kernel.o \
      calc_flux_base.o \
      calc_steps.o \
      set_bc_common.o \
      set.o \
      print.o \
      calc_time_dev.o \
      main.o

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) $^ -o $@

mod_constant.o: mod_globals.mod
set_bc_common.o: mod_globals.mod
set.o: mod_globals.mod mod_constant.mod set_bc_common.mod set_coordinate.mod
calc_physical_quantities.o: mod_globals.mod mod_constant.mod
calc_keep_kernel.o: mod_constant.mod
calc_flux_base.o: mod_globals.mod calc_physical_quantities.mod calc_keep_kernel.mod
calc_steps.o: mod_globals.mod
set_coordinate.o: mod_globals.mod
calc_time_dev.o: mod_globals.mod set.mod calc_flux_base.mod calc_steps.mod print.mod
print.o: mod_globals.mod
main.o: main.f90 mod_globals.mod set_coordinate.mod set.mod calc_time_dev.mod

clean:
	rm -f *.mod *.o *.out

.PHONY: clean
