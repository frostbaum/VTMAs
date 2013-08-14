FC=gfortran
FCFLAGS+= -O2
p_NAME := VTMAs

.PHONY: all, clean, distclean

all: $(p_NAME)

$(p_NAME): v3d_func_rep.o calc_vtma.o
	$(FC) $(FCFLAGS) v3d_func_rep.o calc_vtma.o -o $(p_NAME)

v3d_func_rep.o v3d_func_rep.mod: 3d-vectorOPs/v3d_func_rep.f95
	$(FC) $(FCFLAGS) -c 3d-vectorOPs/v3d_func_rep.f95
	@ touch v3d_func_rep.mod

calc_vtma.o: calc_vtma.f95 v3d_func_rep.mod
	$(FC) $(FCFLAGS) -c calc_vtma.f95

clean:
	@- $(RM) $(p_NAME)
	@- $(RM) *.o *.mod

distclean: clean
