# Intel Fortran on Linux

PROGRAM = droplet

FC = ifort
FCFLAGS = -traceback -CB -g -O0 -fpe0
# FCFLAGS = -qopenmp

OBJS = filename_mod.o csv_reader.o caseList_mod.o path_operator.o vector.o\
	SCTfile_reader.o  vtkMesh_operator.o unstructured_grid.o adjacency_solver.o \
	stl_reader.o adhesion_onSTL.o plot3d_operator.o CUBE_mod.o \
    flow_field.o equation_mod.o virusDroplet_mod.o drop_motion.o \
	dropletManager.o main.o

SRCDIR    = src
OBJDIR    = obj
OBJECTS   = $(addprefix $(OBJDIR)/, $(OBJS))
MODDIR = ${OBJDIR}

$(PROGRAM): $(OBJECTS)
	$(FC) $^ -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	@if [ ! -d $(OBJDIR) ]; then \
		mkdir -p $(OBJDIR); \
	fi
	$(FC) $< -o $@ -c -module $(MODDIR) $(FCFLAGS)

clean:
	$(RM) $(PROGRAM) -r $(OBJDIR)
