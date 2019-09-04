FCFLAGS = -g -O3
LDFLAGS = -lblas
F90 = gfortran

colint: spmat_types.o task_list_types.o task_list_methods.o grid_base_ref.o colint.o
	$(F90) -o colint spmat_types.o task_list_types.o task_list_methods.o \
                  grid_base_ref.o colint.o \
                  $(FCFLAGS) $(LDFLAGS)

.o: %.f90
	$(F90) $(FCFLAGS) -o $@ -c $<

spmat_types.o: spmat_types.f90
	$(F90) $(FCFLAGS) -c spmat_types.f90

task_list_types.o: task_list_types.f90
	$(F90) $(FCFLAGS) -c task_list_types.f90

task_list_methods.o: task_list_methods.f90
	$(F90) $(FCFLAGS) -c task_list_methods.f90

grid_base_ref.o: grid_base_ref.f90
	$(F90) $(FCFLAGS) -c grid_base_ref.f90

colint.o: colint.f90
	$(F90) $(FCFLAGS) -c colint.f90

.PHONY: clean
clean:
	rm -f *.o *.mod colint


