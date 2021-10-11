PROGRAM = droplet

FC = gfortran
FCFLAGS = -O -fbacktrace -g
# FCFLAGS = -Wall -fbounds-check -O -Wuninitialized -fbacktrace -g

FCCOMPILE = ${FC} ${FCFLAGS}

OBJS = csv_reader.o cases_reader.o stl_reader.o fld_reader.o plot3d_operator.o CUBE_mod.o \
    unstructured_grid.o adjacency_solver.o flow_field.o equation_mod.o drop_motion.o \
	main.o

${PROGRAM}: ${OBJS}
	${FCCOMPILE} -o ${PROGRAM} ${OBJS}

%.o:%.f90
	${FCCOMPILE} -c $<

clean:
	- del *.o *.mod
# - rm -f *.o *~ *.mod
