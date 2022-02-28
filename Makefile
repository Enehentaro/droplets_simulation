# Intel Fortran on Linux

MAINTARGET = droplet

FC = ifort

FCFLAGS = -traceback -CB -g -O0 -fpe0 -mcmodel=large
# FCFLAGS = -qopenmp

TARGET1 = CUBE2USG
TARGET2 = droplet2CSV
TARGET3 = dropletCount
TARGET4 = initialTranslate
TARGET5 = boxFlow

COMMONOBJ = filename_mod.o simpleFile_reader.o path_operator.o vector.o terminalControler.o caseName.o conditionValue.o \
    	dropletEquation.o virusDroplet.o
	
MAINOBJ = $(COMMONOBJ) SCTfile_reader.o vtkMesh_operator.o adjacency_solver.o unstructured_grid.o flow_field.o \
			dropletGenerator.o dropletMotionSimulation.o main.o

TARGET1OBJ = simpleFile_reader.o vtkMesh_operator.o plot3d_operator.o CUBE_mod.o CUBE2USG.o
TARGET2OBJ = $(COMMONOBJ) droplet2CSV.o
TARGET3OBJ = $(COMMONOBJ) vtkMesh_operator.o boxCounter.o dropletCount.o
TARGET4OBJ = $(COMMONOBJ) initial_translate.o

TARGET5OBJ = $(COMMONOBJ) SCTfile_reader.o vtkMesh_operator.o adjacency_solver.o unstructured_grid.o flow_field.o \
				boxCounter.o boxFlowField.o

SRCDIR    = src
OBJDIR    = obj
MAINOBJECTS   = $(addprefix $(OBJDIR)/, $(MAINOBJ))
TARGET1OBJECTS   = $(addprefix $(OBJDIR)/, $(TARGET1OBJ))
TARGET2OBJECTS   = $(addprefix $(OBJDIR)/, $(TARGET2OBJ))
TARGET3OBJECTS   = $(addprefix $(OBJDIR)/, $(TARGET3OBJ))
TARGET4OBJECTS   = $(addprefix $(OBJDIR)/, $(TARGET4OBJ))
TARGET5OBJECTS   = $(addprefix $(OBJDIR)/, $(TARGET5OBJ))
MODDIR = ${OBJDIR}

$(MAINTARGET): $(MAINOBJECTS)
	$(FC) $^ $(FCFLAGS) -o $@

$(TARGET1): $(TARGET1OBJECTS)
	$(FC) $^ $(FCFLAGS) -o $@

$(TARGET2): $(TARGET2OBJECTS)
	$(FC) $^ $(FCFLAGS) -o $@

$(TARGET3): $(TARGET3OBJECTS)
	$(FC) $^ $(FCFLAGS) -o $@

$(TARGET4): $(TARGET4OBJECTS)
	$(FC) $^ $(FCFLAGS) -o $@

$(TARGET5): $(TARGET5OBJECTS)
	$(FC) $^ $(FCFLAGS) -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	@if [ ! -d $(OBJDIR) ]; then \
		mkdir -p $(OBJDIR); \
	fi
	$(FC) $< -o $@ -c -module $(MODDIR) $(FCFLAGS)

all: $(MAINTARGET) $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5)

clean:
	$(RM) $(MAINTARGET) $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5) -r $(OBJDIR)
