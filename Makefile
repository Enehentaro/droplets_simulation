# GNU Fortran on Windows
PROGRAM = droplet_Kishi

FC = gfortran
FCFLAGS = -O0 -fbacktrace -g
# -Wall -fbounds-check -Wuninitialized

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
	$(FC) $^ -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	@if not exist $(OBJDIR) ( \
		md $(OBJDIR) \
	)
	$(FC) $< -o $@ -c -J$(MODDIR) ${FCFLAGS}

clean:
	del /Q ${OBJDIR} $(PROGRAM).exe
