
#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.
HOME = ./..
WRITERS = $(HOME)/utils/Writers
PROTO_DIR = $(HOME)/proto/include
FFT = $(HOME)/src/fftTools
CONV = $(HOME)/src/Hockney
PARTICLES = $(HOME)/src/Particles
TIMER = $(HOME)/utils/timer
VPATH= . $(HOME) $(PARTICLES) $(PROTO_DIR) $(FFT) $(CONV) $(TIMER) $(WRITERS)
ANIMATION = 1

DIM=2
CXX=g++
# CXX=clang++
FFTWDIR = /usr/local
#CFLAGS = -g -Wall
CXXFLAGS = -O3
CXXFLAGS += -std=c++11 -I./util.h -I./w.h -I$(PARTICLES) -I$(PROTO_DIR) -I$(FFT) -I$(CONV) -I$(TIMER) -I$(WRITERS) -I$(FFTWDIR)/include

odir = ./o.$(DIM)d
ddir = ./d.$(DIM)d

CXXFLAGS += -DDIM=$(DIM) -D ANIMATION=$(ANIMATION)
LIBS:=  -L$(FFTWDIR)/lib -lfftw3

SRCFILES:= $(notdir $(wildcard $(TIMER)/*.cpp $(PROTO_DIR)/*.cpp $(WRITERS)/*.cpp $(FFT)/*.cpp ./*.cpp $(CONV)/*.cpp $(PARTICLES)/*.cpp))

OBJS:=$(patsubst %.cpp,$(odir)/%.o,$(SRCFILES))


DEPS:=$(patsubst $(odir)/%.o,$(ddir)/%.d, $(OBJS))

$(odir)/%.o:%.cpp GNUmakefile
	mkdir -p $(odir);$(CXX) -c $(CXXFLAGS) $< -o $@
	mkdir -p $(ddir);$(CXX) -MM $(CXXFLAGS) $< | sed '1s/^/o.$(DIM)d\//' > $*.d;mv $*.d $(ddir)

vortex2D: GNUmakefile util.cpp w.cpp $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LIBS) -o execPlot

-include $(DEPS)

clean:
	rm -r *.exe $(odir) $(ddir) *.vtk

listsrc:
	@echo $(SRCFILES)
listobj:
	@echo $(OBJS)
listdep:
	@echo $(DEPS)

