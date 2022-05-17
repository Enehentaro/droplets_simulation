# Intel Fortran on Linux
FC = ifort
FCFLAGS = -traceback -CB -g -O0 -fpe0
# FCFLAGS = -qopenmp

#プログラム名配列
PROGRAMS = droplet CUBE2USG droplet2CSV dropletCount initialTranslate boxFlow

#頻繁に使うファイル
COMMONOBJ = filename_mod simpleFile_reader path_operator vector terminalControler caseName conditionValue \
    	dropletEquation virusDroplet timeKeeper

#プログラムコンパイルに必要な依存ファイル配列
#配列の名前は、プログラム名＋"_OBJ"
#左から順にコンパイルするので、記述順に注意
droplet_OBJ = $(COMMONOBJ) SCTfile_reader vtkMesh_operator adjacency_solver unstructured_grid flow_field \
			dropletGenerator dropletMotionSimulation main
CUBE2USG_OBJ = $(COMMONOBJ) vtkMesh_operator plot3d_operator CUBE2USG
droplet2CSV_OBJ = $(COMMONOBJ) droplet2CSV
dropletCount_OBJ = $(COMMONOBJ) vtkMesh_operator boxCounter dropletCount
initialTranslate_OBJ = $(COMMONOBJ) initial_translate
boxFlow_OBJ = $(COMMONOBJ) SCTfile_reader vtkMesh_operator adjacency_solver unstructured_grid flow_field \
				boxCounter boxFlowField
        
#ディレクトリ指定
SRCDIR = src
OBJDIR = obj
BINDIR = bin
MODDIR = $(OBJDIR)


#ディレクトリが存在しなければ作成する関数
define mkdirIfNotExist
	@if [ ! -d ${1} ]; then \
		mkdir -p ${1}; \
	fi
endef

#ルール（レシピ？）を記述するための関数
#分割コンパイルされたオブジェクトをリンクし、実行ファイルをbinディレクトリに出力
define makeRule
${1}: ${2}
	$(call mkdirIfNotExist,$(BINDIR))
	$(FC) ${2} $(FCFLAGS) -o $(BINDIR)/${1}
endef

#これが一番重要
#プログラム名配列が展開され、各プログラムのルールを記述する関数が呼ばれる
$(foreach prg,${PROGRAMS},$(eval $(call makeRule,$(prg),$(addprefix $(OBJDIR)/,$(addsuffix .o,$($(prg)_OBJ))))))

#オブジェクトファイルのレシピ
#ソースコードを分割コンパイルし、オブジェクトファイルをobjディレクトリに作成
$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(call mkdirIfNotExist,$(OBJDIR))
	$(FC) $< -o $@ -c -module $(MODDIR) $(FCFLAGS)


all: $(PROGRAMS)

clean:
	$(RM) -r $(OBJDIR) $(BINDIR)
