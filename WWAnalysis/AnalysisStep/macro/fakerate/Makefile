ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs) -lTMVA

CXX           = g++
CXXFLAGS      = -g -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE -g -O2
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared


ARCH         := $(shell root-config --arch)
PLATFORM     := $(shell root-config --platform)


CXXFLAGS      += $(ROOTCFLAGS)
#CXX           += -I./
LIBS           = $(ROOTLIBS) 

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit
GLIBS          = $(filter-out -lNew, $(NGLIBS))

INCLUDEDIR       = ./ 
CXX	         += -I$(INCLUDEDIR)
OUTLIB	         = lib/

.SUFFIXES: .cc,.C,.hh,.h
.PREFIXES: ./lib/

$(OUTLIB)EfficiencyEvaluator.o: src/EfficiencyEvaluator.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)EfficiencyEvaluator.o $<

fakerate: src/fakerateHzz4l.cc \
	src/estimateMuonFakeRateHzz4lTree.C \
	src/estimateElectronFakeRateHzz4lTree.C \
	$(OUTLIB)EfficiencyEvaluator.o
	$(CXX) $(CXXFLAGS) -o fakerate $(OUTLIB)*.o $(GLIBS) $ $<

clean:
	rm -f $(OUTLIB)*.o
	rm -f fakerate

all: fakerate

