PROGRAM = droplet

FC = ifort
FCFLAGS = -traceback -CB -g -O0 -fpe0
# FCFLAGS = -qopenmp

# FC = gfortran
# FCFLAGS = -O -fbacktrace -g
# FCFLAGS = -Wall -fbounds-check -O -Wuninitialized -fbacktrace -g

FCCOMPILE = ${FC} ${FCFLAGS}

OBJS = filename_mod.o csv_reader.o caseList_mod.o path_operator.o vector.o\
	fld_reader.o  vtkMesh_operator.o unstructured_grid.o adjacency_solver.o \
	stl_reader.o adhesion_onSTL.o plot3d_operator.o CUBE_mod.o \
    flow_field.o equation_mod.o virusDroplet_mod.o drop_motion.o \
	dropletManager.o main.o

SRCDIR    = src
OBJDIR    = obj
OBJECTS   = $(addprefix $(OBJDIR)/, $(OBJS))
MODDIR = ${OBJDIR}

$(PROGRAM): $(OBJECTS)
	$(FC) -o $@ $^

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	@if [ ! -d $(OBJDIR) ]; then \
		mkdir -p $(OBJDIR); \
	fi
	$(FCCOMPILE) -o $@ -c $< -module $(MODDIR)
# $(FCCOMPILE) -o $@ -c $<

clean:
	- rm -f $(PROGRAM) -r $(OBJDIR)
# - del /Q ${OBJDIR}\*.o *.mod *.exe
# - del *.o *.mod *.exe
