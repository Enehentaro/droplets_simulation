PROGRAM = droplet

FC = gfortran
FCFLAGS = -O -fbacktrace -g
# FCFLAGS = -Wall -fbounds-check -O -Wuninitialized -fbacktrace -g

FCCOMPILE = ${FC} ${FCFLAGS}

OBJS = filename_mod.o csv_reader.o caseList_mod.o path_operator.o stl_reader.o adhesion_onSTL.o \
	fld_reader.o plot3d_operator.o CUBE_mod.o vtkMesh_operator.o\
    unstructured_grid.o adjacency_solver.o flow_field.o equation_mod.o drop_motion.o \
	main.o

${PROGRAM}: ${OBJS}
	${FC} -o ${PROGRAM} ${OBJS}

%.o: %.f90
	${FCCOMPILE} -o $@ -c $<

clean:
	- del *.o *.mod *.exe
# - rm -f *.o *~ *.mod
