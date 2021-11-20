# GNU Fortran on Windows

TARGET = droplet_Kishi

FC = gfortran
FCFLAGS = -O0 -fbacktrace -g
# -Wall -fbounds-check -Wuninitialized

TARGET1 = CUBE2USG
TARGET2 = droplet2CSV

OBJS = filename_mod.o csv_reader.o caseList_mod.o path_operator.o vector.o terminalControler.o\
	SCTfile_reader.o  vtkMesh_operator.o unstructured_grid.o adjacency_solver.o \
	stl_reader.o adhesion_onSTL.o plot3d_operator.o CUBE_mod.o \
    flow_field.o equation_mod.o virusDroplet_mod.o drop_motion.o
	
MAINOBJS = dropletManager.o main.o

SRCDIR    = src
OBJDIR    = obj
OBJECTS   = $(addprefix $(OBJDIR)/, $(OBJS))
MAINOBJECTS   = $(addprefix $(OBJDIR)/, $(MAINOBJS))
MODDIR = ${OBJDIR}

$(TARGET): $(OBJECTS) $(MAINOBJECTS)
	$(FC) $^ -o $@

$(TARGET1): $(OBJECTS) $(OBJDIR)/$(TARGET1).o
	$(FC) $^ -o $@

$(TARGET2): $(OBJECTS) $(OBJDIR)/$(TARGET2).o
	$(FC) $^ -o $@		

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	@if not exist $(OBJDIR) ( \
		md $(OBJDIR) \
	)
	$(FC) $< -o $@ -c -J$(MODDIR) ${FCFLAGS}

all: $(TARGET) $(TARGET1) $(TARGET2)

clean:
	del /Q ${OBJDIR} $(TARGET).exe $(TARGET1).exe $(TARGET2).exe
