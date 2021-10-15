PROGRAM = droplet

FC = ifort
FCFLAGS = -traceback -CB -g -O0
# FCFLAGS = -qopenmp

# FC = gfortran
# FCFLAGS = -O -fbacktrace -g
# FCFLAGS = -Wall -fbounds-check -O -Wuninitialized -fbacktrace -g

FCCOMPILE = ${FC} ${FCFLAGS}

OBJS = filename_mod.o csv_reader.o caseList_mod.o path_operator.o stl_reader.o adhesion_onSTL.o \
	fld_reader.o plot3d_operator.o CUBE_mod.o vtkMesh_operator.o\
    unstructured_grid.o adjacency_solver.o flow_field.o equation_mod.o drop_motion.o \
	main.o

SRCDIR    = src
OBJDIR    = obj
OBJECTS   = $(addprefix $(OBJDIR)/, $(OBJS))
MODDIR = ${SRCDIR}

$(PROGRAM): $(OBJECTS)
	$(FC) -o $@ $^

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	@if [ ! -d $(OBJDIR) ]; then \
		mkdir -p $(OBJDIR); \
	fi
	$(FCCOMPILE) -o $@ -c $< -module $(MODDIR)
# $(FCCOMPILE) -o $@ -c $<

clean:
	- rm -f -r $(OBJDIR) $(SRCDIR)/*.mod $(PROGRAM)
# - del /Q ${OBJDIR}\*.o *.mod *.exe
# - del *.o *.mod *.exe
