#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.
DIM=2
PROTO_DIR = /Users/colella/Desktop/gitstuff/resources/proto/include
VPATH = ../timer $(PROTO_DIR) ../Writers
system := $(shell uname)
FFTWDIR = /usr/local
FFTLIBLOC = ../lib
CXX=clang++
#CXX=g++
FFTINCFLAGS = -I$(FFTWDIR)/include
#CXXFLAGS = -g -Wall -DDIM=$(DIM)
CXXFLAGS = -O3 -DDIM=$(DIM)
CXXFLAGS += -std=c++11 -I. -I../Writers -I../timer -I$(PROTO_DIR) $(FFTINCFLAGS)
ifeq ($(CXX) , clang++)
CXXFLAGS += -stdlib=libc++
endif
LIBS =  -L$(FFTWDIR)/lib -lfftw3 

LIBTOOL = ar cru 
ifeq ($(system),Darwin)
  LIBTOOL = libtool -static -s -o
endif

osuffix:=$(DIM)D.o
dsuffix:=$(DIM)D.d
%.o: %.cpp GNUmakefile
	$(CXX) -c $(CXXFLAGS) $< -o $@
	$(CXX) -MM $(CXXFLAGS) $< > $*.d

%$(osuffix): %.cpp GNUmakefile
	$(CXX) -c $(CXXFLAGS) $< -o $@
	$(CXX) -MM $(CXXFLAGS) $< > $*$(dsuffix)

1DFFTOBJS = FFT1DBRI.o FFTCTBRI.o FFTW1D.o FFT1DRecursive.o 
#add these files when they are implemented.

ifneq ($(DIM),1)
MDFFTOBJS = FFTMD$(osuffix)
libfft$(DIM)D.a: GNUmakefile $(MDFFTOBJS)
	$(LIBTOOL) libfft$(DIM)D.a $(MDFFTOBJS) 
	mkdir -p ../lib;mv libfft$(DIM)D.a ../lib
endif

libfft1D.a: GNUmakefile $(1DFFTOBJS)
	$(LIBTOOL) libfft1D.a $(1DFFTOBJS)
	mkdir -p ../lib;mv libfft1D.a ../lib

clean:
	rm *.o *.d *.a
-include $(MDFFTOBJS:.o=.d)