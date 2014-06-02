#!/usr/bin/make -f  
_VERSION_=
#(64 bit) 1.1.0.0

INCPATH=-I./ -I/usr/include 
LIBPATH=-L./ 
# c++0x supports array initialization in constructor
FLAGS=-DLINUX -Wall -O3 
DEBUG=0
CORONA=0

LIBS = $(LIBPATH) \
	   -lglut -lGLU -lGL -lXmu -lXext -lX11 -lm \
	   -lpthread -lrt 

CC=g++

LIBSRCS= GLoader.cpp

ifneq ($(DEBUG), 0)
	FLAGS=-Wno-deprecated -DDEBUG -g _VERSION_ += (debug version)
endif

LIBOBJ=$(LIBSRCS:.cpp=.o)

%.o : %.cpp
	$(CC) -c $< -o $@ $(INCPATH) $(FLAGS) -DVERSION=\""$(_VERSION_)\""

###############################################################################

EXE  = $(t)
SRCS = $(EXE).cpp $(LIBSRCS)
OBJS = $(SRCS:.cpp=.o) 

$(EXE) : $(OBJS)
	$(CC) -o $@ $< $(LIBOBJ) $(LIBS)
	-@rm -rf $(OBJS)

.PHONY : single 
single : $(EXE)
	@echo "Single target found got built."; 

.PHONY : clean
clean : 
	-@rm -rf $(t) $(OBJS) *~

.PHONY : rebuild
rebuild : clean single

.PHONY : release
release : rebuild

###############################################################################

EXES = stereoTesting
ALLSRCS = $(foreach x,$(EXES), ${x}.cpp) $(LIBSRCS)
ALLOBJS = $(ALLSRCS:.cpp=.o) 

$(EXES) : $(foreach x,$(EXES),${x}.o) $(LIBOBJ)
	$(CC) -o $@ $@.o $(LIBOBJ) $(LIBS)
	-@rm -rf $@.o

.PHONY : all
all : $(EXES)
	@echo "All found got built."; 

.PHONY : cleanall 
cleanall : 
	-@rm -rf $(EXES) $(ALLOBJS)

.PHONY : rebuildall
rebuild : cleanall all

.PHONY : releaseall
release : rebuildall

###############################################################################

help : usage

.PHONY : usage
usage : 
	@echo " [ Usage ] "
	@echo " build single :       make t=<target> | clean | release | rebuild"
	@echo " build all    :       make all | cleanall | releaseall | rebuildall" 
	@echo

