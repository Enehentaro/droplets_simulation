# Intel Fortran on Linux

MAINTARGET = droplet

FC = ifort

FCFLAGS = -traceback -CB -g -O0 -fpe0
# FCFLAGS = -qopenmp

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

SRCDIR = src
OBJDIR = obj
BINDIR = bin
MODDIR = ${OBJDIR}
MAINOBJECTS = $(addprefix $(OBJDIR)/, $(MAINOBJ))
TARGET1OBJECTS = $(addprefix $(OBJDIR)/, $(TARGET1OBJ))
TARGET2OBJECTS = $(addprefix $(OBJDIR)/, $(TARGET2OBJ))
TARGET3OBJECTS = $(addprefix $(OBJDIR)/, $(TARGET3OBJ))
TARGET4OBJECTS = $(addprefix $(OBJDIR)/, $(TARGET4OBJ))

$(MAINTARGET): $(MAINOBJECTS)
	$(call linkCompile,$^,$@)

$(TARGET1): $(TARGET1OBJECTS)
	$(call linkCompile,$^,$@)

$(TARGET2): $(TARGET2OBJECTS)
	$(call linkCompile,$^,$@)

$(TARGET3): $(TARGET3OBJECTS)
	$(call linkCompile,$^,$@)

$(TARGET4): $(TARGET4OBJECTS)
	$(call linkCompile,$^,$@)

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(call mkdirIfNotExist,$(OBJDIR))
	$(FC) $< -o $@ -c -module $(MODDIR) $(FCFLAGS)

all: $(MAINTARGET) $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4)

clean:
	$(RM) -r $(OBJDIR) $(BINDIR)

#===========以下、自作関数=============================================

#ディレクトリが存在しなければ作成する
define mkdirIfNotExist
	@if [ ! -d ${1} ]; then \
		mkdir -p ${1}; \
	fi
endef

#分割コンパイルされたオブジェクトをリンクし、実行ファイルをbinディレクトリに出力
define linkCompile
	$(call mkdirIfNotExist,$(BINDIR))
	$(FC) ${1} $(FCFLAGS) -o $(BINDIR)/${2}
endef
