# Intel Fortran on Linux

TARGET = droplet

FC = ifort
FCFLAGS = -traceback -CB -g -O0 -fpe0
#FCFLAGS = -qopenmp

TARGET1 = CUBE2USG
TARGET2 = droplet2CSV
TARGET3 = dropletCount
TARGET4 = initialTranslate

OBJS = filename_mod.o simpleFile_reader.o path_operator.o vector.o terminalControler.o caseName.o  conditionValue.o\
	SCTfile_reader.o  vtkMesh_operator.o unstructured_grid.o adjacency_solver.o \
    flow_field.o dropletEquation.o virusDroplet.o dropletGroup.o dropletMotionSimulation.o
	
MAINOBJS = dropletManager.o main.o

TARGET1OBJS = simpleFile_reader.o vtkMesh_operator.o plot3d_operator.o CUBE_mod.o CUBE2USG.o
TARGET3OBJS = boxCounter.o dropletCount.o
TARGET4OBJS = initial_translate.o

SRCDIR    = src
OBJDIR    = obj
OBJECTS   = $(addprefix $(OBJDIR)/, $(OBJS))
MAINOBJECTS   = $(addprefix $(OBJDIR)/, $(MAINOBJS))
TARGET1OBJECTS   = $(addprefix $(OBJDIR)/, $(TARGET1OBJS))
TARGET3OBJECTS   = $(addprefix $(OBJDIR)/, $(TARGET3OBJS))
TARGET4OBJECTS   = $(addprefix $(OBJDIR)/, $(TARGET4OBJS))
MODDIR = ${OBJDIR}

$(TARGET): $(OBJECTS) $(MAINOBJECTS)
	$(FC) $^ $(FCFLAGS) -o $@

$(TARGET1): $(TARGET1OBJECTS)
	$(FC) $^ $(FCFLAGS) -o $@

$(TARGET2): $(OBJECTS) $(OBJDIR)/$(TARGET2).o
	$(FC) $^ $(FCFLAGS) -o $@

$(TARGET3): $(OBJECTS) $(TARGET3OBJECTS)
	$(FC) $^ $(FCFLAGS) -o $@

$(TARGET4): $(OBJECTS) $(TARGET4OBJECTS)
	$(FC) $^ $(FCFLAGS) -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	@if [ ! -d $(OBJDIR) ]; then \
		mkdir -p $(OBJDIR); \
	fi
	$(FC) $< -o $@ -c -module $(MODDIR) $(FCFLAGS)

all: $(TARGET) $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4)

clean:
	$(RM) $(TARGET) $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) -r $(OBJDIR)
