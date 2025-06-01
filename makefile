UNAME_S := $(shell uname -s)
# Detect Windows-like environments (MinGW or Cygwin)
ifneq (,$(filter MINGW% CYGWIN%,$(UNAME_S)))
    MIMP_ON_WINDOWS := 1
else
    MIMP_ON_LINUX := 1
endif

MIMP_SOURCE_PATH = ../src

CSCOPE_EXE = c:/msys64/usr/bin/cscope.exe
GTAGS_EXE = c:/msys64/mingw64/bin/gtags

export GTAGSFORCECPP = 1
export TMPDIR = c:/tmp
export TMP = c:/tmp

vpath %.c $(MIMP_SOURCE_PATH)
vpath %.cpp $(MIMP_SOURCE_PATH)
vpath %.h $(MIMP_SOURCE_PATH)
vpath %.ispc $(MIMP_SOURCE_PATH)

MFILE = makefile

MAKELIB = $(AR) rv

CFLAGS_MT = -D_REENTRANT 
LFLAGS_MT = -lpthread

CFLAGS_OPT = -O3
#CFLAGS_OPT = -g

ifdef MIMP_ON_LINUX
CFLAGS_PLATFORM = -DMIMP_ON_LINUX
else
CFLAGS_PLATFORM = -DMIMP_ON_WINDOWS
endif

#-------------------------------------------------------------------------

CFLAGS = $(CFLAGS_OPT) $(CFLAGS_PROF) $(CFLAGS_PLATFORM)\
	 $(CFLAGS_MIMP_ENV) \
	 $(CFLAGS_TBB) -DTBB_SUPPRESS_DEPRECATED_MESSAGES\
	 -Wall $(CFLAGS_EXT) \
	 $(CFLAGS_WARN) \
	 -Wstrict-aliasing=2 \
	 -Wno-unused-label \
	 -Wno-unused-function \
	 -Wno-unused-but-set-variable \
	 -DANN_WITH_FLOAT $(CFLAGS_MT) \
	 -ffast-math -fno-finite-math-only \
	 $(CFLAGS_OMP) \
	 -I. -I$(MIMP_SOURCE_PATH)/../external

# misc. additional flags 
# CFLAGS3 = -Wa,-mbig-obj
#CFLAGS3 = -E
#CFLAGS3 = -S -fverbose-asm 
#CFLAGS3 = -Wshadow
#CFLAGS3 = -save-temps
#CFLAGS3 = -ftree-vectorizer-verbose=2
#CFLAGS3 = -ftree-vectorizer-verbose=10
#CFLAGS3 = -fno-tree-vectorize -ftree-vectorizer-verbose=2

CPP_FLAGS_ALL = $(CFLAGS) $(CFLAGS2) $(CFLAGS_BW) $(CFLAGS3) $(CPP_FLAGS)

LDFLAGS =

CTAGS = $(MIMP_SOURCE_PATH)/tags
GTAGS =  $(MIMP_SOURCE_PATH)/GTAGS
CSCOPE = $(MIMP_SOURCE_PATH)/cscope.out
#
# source files
#
#

SOURCE = main.cpp msbg_demo.cpp gwx.c mtool.c util.c util2.cpp plot.c rand.c panel.c \
	 readpng.c sbg.cpp thread.cpp fastmath.cpp blockpool.cpp grid.c \
	 bitmap.c bitmap2.c pnoise.c pnoise2.cpp pnoise3.cpp halo.cpp msbg.cpp \
	 msbg2.cpp msbg3.cpp msbg4.cpp visualizeSlices.cpp render.cpp \
	 msbgaux.cpp mm.c mm2.c 

HEADER = MGWin.h VEC3.h bitmap.h blockpool.h common_headers.h fastmath.h globdef.h \
	grid.h gwx.h halo.h kernels_ispc.h log.h mm.h msbg.h msbg_demo.h msbgcell.h \
	mtool.h	mtool2.h noise.h panel.h plot.h pnoise.h rand.h readpds.h render.h \
	sbg.h stdmath.h thread.h util.h vectorclass_util.h vectorclass_util2.h \
	windows_linux.h wturbglob.h internal.h

SOURCE_GTAGS = $(SOURCE) $(HEADER)

LIBNAME = msbg

STATIC_LIB = lib$(LIBNAME).a
SHARED_LIB = lib$(LIBNAME).so

%.$(OBJE) : %.cpp 
	$(CC) -c $(CPP_FLAGS_ALL) \
	      -o $@ $<

%.$(OBJE) : %.c 
	$(CC) -c $(CFLAGS) $(CFLAGS2) $(CFLAGS_BW) $(CFLAGS3) \
	      -o $@ $<

all: msbg_demo$(EXE) 

#test: 
#	@echo "MIMP_ON_LINUX="$(MIMP_ON_LINUX)
#	@echo "MIMP_ON_WINDOWS="$(MIMP_ON_WINDOWS)


ifdef MIMP_ON_LINUX
  LDFLAGS_GCC:=
else  

ifeq ($(CC),clang)
  LDFLAGS_GCC:=-mwindows
else
  LDFLAGS_GCC:=-Wl,--subsystem,console -mwindows -lgomp
endif

endif


ifdef MIMP_ON_LINUX
  LDFLAGS_PLATFORM:=
else
  LDFLAGS_PLATFORM:=-lwinmm
endif

# 
# library
#
OBJS_MSBG_LIB = gwx.$(OBJE) mtool.$(OBJE) util.$(OBJE) util2.$(OBJE) \
		 plot.$(OBJE) rand.$(OBJE) panel.$(OBJE) readpng.$(OBJE) sbg.$(OBJE) \
		 thread.$(OBJE) fastmath.$(OBJE) blockpool.$(OBJE) grid.$(OBJE) bitmap.$(OBJE) bitmap2.$(OBJE) \
		 pnoise.$(OBJE) pnoise2.$(OBJE) pnoise3.$(OBJE) halo.$(OBJE) msbg.$(OBJE) msbg2.$(OBJE) msbg3.$(OBJE) msbg4.$(OBJE) \
		 visualizeSlices.$(OBJE) render.$(OBJE) msbgaux.$(OBJE) \
		 mm.$(OBJE) mm2.$(OBJE) 

$(STATIC_LIB): $(OBJS_MSBG_LIB)
	ar rcs $@ $^

$(SHARED_LIB): $(OBJS_MSBG_LIB)
	$(CXX) -shared -o $@ $^

# 
# msbg_demo
#

OBJS_MSBG_DEMO = main.$(OBJE) msbg_demo.$(OBJE) 

LD_LIBS_FOR_MSBG_DEMO = \
	    -lpng \
	    -ljpeg \
	    -ltbbmalloc -ltbb -ltbbmalloc_proxy \
	    -lz \
	    -lm 

ifdef MIMP_ON_LINUX

msbg_demo$(EXE): $(OBJS_MSBG_DEMO) $(STATIC_LIB)
	$(LD) $(LDFLAGS) $@ \
	  $(OBJS_MSBG_DEMO) \
	  -L. -l$(LIBNAME) \
	  $(LD_LIBS_FOR_MSBG_DEMO) \
		$(LDFLAGS_BW) \
		$(LDFLAGS_PROF) \
		-lgomp \
		$(LFLAGS_MT)
else
# 
# native windows application via MSYS2/MinGw64
#
msbg_demo$(EXE): $(OBJS_MSBG_DEMO) $(STATIC_LIB)
	$(LD) $(LDFLAGS) $@ \
	  -static-libgcc -static-libstdc++ \
	  $(OBJS_MSBG_DEMO) \
	  -L. -l$(LIBNAME) \
	  $(LD_LIBS_FOR_MSBG_DEMO) \
	  $(LDFLAGS_BW) \
	  $(LDFLAGS_PROF) \
	  -lgomp \
	  -lwinmm \
	  -Wl,--enable-auto-import \
	  -Wl,--subsystem,console -mwindows \
	  $(LFLAGS_MT)
endif



# 
# clean
#
OBJS_ALL = $(OBJS_MIMP) $(OBJS_TSTMKL) $(OBJS_TSTENV) \
	   $(OBJS_ISPC) $(OBJS_MSBG_DEMO)

clean:
	rm $(OBJS_ALL)
	rm msbg_demo$(EXE)

gtags_files:
	@rm -f $(MIMP_SOURCE_PATH)/gtags.files
	for F in $(SOURCE_GTAGS); do \
	  echo $$F >> $(MIMP_SOURCE_PATH)/gtags.files; \
	done  

gtags_make: gtags_files
	$(GTAGS_EXE) -C $(MIMP_SOURCE_PATH) -f $(MIMP_SOURCE_PATH)/gtags.files

gtags_update: 
	$(GTAGS_EXE) -C $(MIMP_SOURCE_PATH) -i

gtags_test: $(GTAGS)


$(CTAGS): $(SOURCE) $(ISPC_SRC) $(HEADER)
	$(CTAGS_EXE) -f $(CTAGS) \
	  $(addprefix $(shell cygpath -m $(MIMP_SOURCE_PATH))/, $(SOURCE) $(ISPC_SRC) $(HEADER))\


$(CSCOPE): $(SOURCE) $(ISPC_SRC) $(HEADER)
	cd $(MIMP_SOURCE_PATH) && $(CSCOPE_EXE) -bk -f  $(shell cygpath -m $(CSCOPE)) $(SOURCE) $(ISPC_SRC) $(HEADER)

#
#
#
#

asm_mark: all
	objdump -d -M intel mimp.exe | sed -n '/1111111111111111/,/2222222222222222/p'	

#  
#  Dependencies (automatically generated via "make depend")
#

depend_test:
	$(CC) $(CPPFLAGS) -std=gnu++17 -MG -MM $(MIMP_SOURCE_PATH)/msbgaux.cpp

depend: gtags_make
	rm -f $(MIMP_SOURCE_PATH)/dependencies.mk
	for F in $(SOURCE); do \
	  D=`dirname $$F | sed "s/^\.\///"`; \
	  echo -n "$$D/" >> $(MIMP_SOURCE_PATH)/dependencies.mk; \
	  $(CC) $(CPPFLAGS) -std=gnu++17 -MG -MM $$F \
	  >> $(MIMP_SOURCE_PATH)/dependencies.mk; \
	done
	for F in $(PCH_CPP_HEADER); do \
	  D=`dirname $$F | sed "s/^\.\///"`; \
	  echo -n "$$D/" >> $(MIMP_SOURCE_PATH)/dependencies.mk; \
	  $(CC) $(CPPFLAGS) -std=gnu++17 -MG -MM -MT $(PCH_CPP_HEADER_0).h -MT $(PCH_CPP_OUT) -c $$F \
	  >> $(MIMP_SOURCE_PATH)/dependencies.mk; \
	done

#include $(MIMP_SOURCE_PATH)/dependencies.mk

