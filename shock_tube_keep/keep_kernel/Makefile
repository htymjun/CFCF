# compiler
FC = nvfortran

# options
FFLAGS = -cuda -acc -fast -mfma -Mchkptr -Minfo -gpu=ptxinfo,rdc,lto

TARGET = sod

SRC = parameters.f90 solver.f90 main.f90

OBJ = $(SRC:.f90=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod $(TARGET)

