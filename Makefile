# GNU Fortran on Windows

TARGET = droplet

FC = gfortran
FCFLAGS = #-O0 -fbacktrace -g

TARGET1 = CUBE2USG
TARGET2 = droplet2CSV
TARGET3 = dropletCount

OBJS = filename_mod.o simpleFile_reader.o caseName.o path_operator.o vector.o terminalControler.o \
	SCTfile_reader.o  vtkMesh_operator.o unstructured_grid.o adjacency_solver.o \
    flow_field.o dropletEquation.o virusDroplet.o dropletGroup.o dropletMotionSimulation.o
	
MAINOBJS = dropletManager.o main.o

TARGET1OBJS = simpleFile_reader.o vtkMesh_operator.o plot3d_operator.o CUBE_mod.o CUBE2USG.o
TARGET3OBJS = boxCounter.o dropletCount.o

SRCDIR    = src
OBJDIR    = obj
OBJECTS   = $(addprefix $(OBJDIR)/, $(OBJS))
MAINOBJECTS   = $(addprefix $(OBJDIR)/, $(MAINOBJS))
TARGET1OBJECTS   = $(addprefix $(OBJDIR)/, $(TARGET1OBJS))
TARGET3OBJECTS   = $(addprefix $(OBJDIR)/, $(TARGET3OBJS))
MODDIR = ${OBJDIR}

$(TARGET): $(OBJECTS) $(MAINOBJECTS)
	$(FC) $^ $(FCFLAGS) -o $@

$(TARGET1): $(TARGET1OBJECTS)
	$(FC) $^ $(FCFLAGS) -o $@

$(TARGET2): $(OBJECTS) $(OBJDIR)/$(TARGET2).o
	$(FC) $^ $(FCFLAGS) -o $@

$(TARGET3): $(OBJECTS) $(TARGET3OBJECTS)
	$(FC) $^ $(FCFLAGS) -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	@if not exist $(OBJDIR) ( \
		md $(OBJDIR) \
	)
	$(FC) $< -o $@ -c -J$(MODDIR) ${FCFLAGS}

all: $(TARGET) $(TARGET1) $(TARGET2) $(TARGET3)

clean:
	del /Q ${OBJDIR} *.exe
