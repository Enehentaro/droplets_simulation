# GNU Fortran on Windows

<<<<<<< HEAD
TARGET = droplet_Kishi

FC = gfortran
FCFLAGS = -O0 -fbacktrace -g
=======
MAINTARGET = droplet

FC = ifort

FCFLAGS = -traceback -CB -g -O0 -fpe0
# FCFLAGS = -qopenmp
>>>>>>> origin

TARGET1 = CUBE2USG
TARGET2 = droplet2CSV
TARGET3 = dropletCount
TARGET4 = initialTranslate

COMMONOBJ = filename_mod.o simpleFile_reader.o path_operator.o vector.o terminalControler.o caseName.o conditionValue.o \
    	dropletEquation.o virusDroplet.o
	
MAINOBJ = $(COMMONOBJ) SCTfile_reader.o vtkMesh_operator.o adjacency_solver.o unstructured_grid.o flow_field.o \
			dropletGenerator.o dropletMotionSimulation.o main.o

TARGET1OBJ = simpleFile_reader.o vtkMesh_operator.o plot3d_operator.o CUBE_mod.o CUBE2USG.o
TARGET2OBJ = $(COMMONOBJ) droplet2CSV.o
TARGET3OBJ = $(COMMONOBJ) vtkMesh_operator.o boxCounter.o dropletCount.o
TARGET4OBJ = $(COMMONOBJ) initial_translate.o

SRCDIR    = src
OBJDIR    = obj
MAINOBJECTS   = $(addprefix $(OBJDIR)/, $(MAINOBJ))
TARGET1OBJECTS   = $(addprefix $(OBJDIR)/, $(TARGET1OBJ))
TARGET2OBJECTS   = $(addprefix $(OBJDIR)/, $(TARGET2OBJ))
TARGET3OBJECTS   = $(addprefix $(OBJDIR)/, $(TARGET3OBJ))
TARGET4OBJECTS   = $(addprefix $(OBJDIR)/, $(TARGET4OBJ))
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

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	@if not exist $(OBJDIR) ( \
		md $(OBJDIR) \
	)
	$(FC) $< -o $@ -c -J$(MODDIR) ${FCFLAGS}

all: $(MAINTARGET) $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4)

clean:
<<<<<<< HEAD
	del /Q ${OBJDIR} $(TARGET).exe $(TARGET1).exe $(TARGET2).exe
=======
	$(RM) $(MAINTARGET) $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) -r $(OBJDIR)
>>>>>>> origin
