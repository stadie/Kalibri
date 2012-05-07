CC=g++
LD=g++
GCC_VERSION=$(shell g++ --version | head -1 | cut -f 3 -d' ' | cut -f 1 -d'.')
ifeq ($(GCC_VERSION), 3)
 F77=g77
 F77LDFLAGS=-lg2c
else
 F77=gfortran
 F77LDFLAGS=-lgfortran
endif


#O2 for optimization, g for debugging, pg for profiling
SPECIALFLAGS= -fpic -g -Wall -O2 #-O1 #-pg# -O2
ROOTAUXCFLAGS=$(shell root-config --auxcflags)
ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs) -lMinuit2
LFLAGS= $(SPECIALFLAGS) -lz $(F77LDFLAGS) -lgsl -lgslcblas -lm
BOOSTLINKFLAGS=-lboost_thread -lpthread
# change path for MacPort or fink@MacOS
 ifeq (exists, $(shell [ -d /opt/local/include/boost ] && echo exists)) 
 BOOSTFLAGS=-I/opt/local/include/boost -I/opt/local/include
 LFLAGS += -L/opt/local/lib
 BOOSTLINKFLAGS=-lboost_thread-mt -lpthread
else ifeq (exists, $(shell [ -d /sw/include/boost ] && echo exists)) 
 BOOSTFLAGS=-I/sw/include/boost -I/sw/include
 SPECIALFLAGS += -arch i386
 LFLAGS += -L/sw/lib
else
 BOOSTFLAGS=-I/usr/include/boost
endif 
RCXX=$(SPECIALFLAGS) $(BOOSTFLAGS) -I. $(ROOTCFLAGS) -DSTANDALONE 
RLXX=$(LFLAGS) $(ROOTLIBS) $(BOOSTLINKFLAGS)  #-lrt -lpthread # -lposix4
ROOTSYS=$(shell root-config --prefix)
KALIBRIDIR=$(PWD)

#other .cc files
OTHERSRCS=ControlPlotsComparison.cc caliber.cc compareControlPlots.cc JetMETCorFactorsFactory.cc toy.cc

#all files that end up in the Kalibri library:
SRCS=$(filter-out $(OTHERSRCS),$(wildcard *.cc))

OBJS = $(SRCS:.cc=.o)

.PHONY: clean bins libs plugins all

all: libs bins

bin:
	@mkdir -p bin
lib:
	@mkdir -p lib
tmp:
	@mkdir -p tmp

clean:
	@rm -rf ti_files tmp lib bin share
	@rm -f *~ *.o *# *.d *.bkp junk caliber libKalibri.so *.cfi fort.* .#*

libs: lib include/lbfgs.h lib/libKalibri.so lib/liblbfgs.so

bins: bin bin/junk bin/caliber

plugins: lib PUReweighting lib/libJetMETCor.so 

lbfgs.o: lbfgs.F
	$(F77) $(RCXX) -fno-automatic -fno-backslash -O -c lbfgs.F

lib/libKalibri.so: $(OBJS) lbfgs.o
	$(LD) $(RCXX) -shared $^ $(RLXX) -o lib/libKalibri.so
	@echo '-> Kalibri library created.'

include/lbfgs.h lib/liblbfgs.a lib/liblbfgs.so: liblbfgs-1.10
	cd liblbfgs-1.10 && $(MAKE) && $(MAKE) install
	@echo '-> shared library lib/liblbfgs-1.10.so created.'

liblbfgs-1.10:
	wget --no-check-certificate https://github.com/downloads/chokkan/liblbfgs/liblbfgs-1.10.tar.gz
	tar zxvf liblbfgs-1.10.tar.gz
	cd liblbfgs-1.10 && ./configure --prefix=$(KALIBRIDIR)

bin/junk: $(OBJS) lbfgs.o caliber.o lib/liblbfgs.a
	$(LD) $^ $(RLXX) lib/liblbfgs.a -o bin/junk
	@ln -s -f bin/junk
	@echo '-> static Kalibri executable created.'

bin/caliber: caliber.o lib/libKalibri.so lib/liblbfgs.so
	$(LD) caliber.o $(RLXX) -Llib -lKalibri -llbfgs -o bin/caliber
	@ln -f -s bin/caliber
	@echo '-> shared Kalibri executable created.'

#special targets:
toy:	ToyMC.o toy.o ConfigFile.o
	$(LD) ToyMC.o toy.o ConfigFile.o $(RLXX) -o toy
	@echo '-> toy MC executable created.'

comp: 	ControlPlotsComparison.o compareControlPlots.o
	$(LD) ControlPlotsComparison.o compareControlPlots.o $(RLXX) -o compControlPlots
	@echo '-> Comparison executable created. Type "compControlPlots" to compare control plots.'

#target for plugins:
lib/libJetMETCor.so: lib/libJetMETObjects.so lib/libKalibri.so JetMETCorFactorsFactory.o
	$(LD) $(RCXX) -shared JetMETCorFactorsFactory.o lib/libJetMETObjects.so lib/libKalibri.so $(RLXX) -o lib/libJetMETCor.so
	@echo '-> JetMETCor plugin created.'

lib/libJetMETObjects.so: bin lib tmp JetMETObjects 
	cd JetMETObjects && $(MAKE) STANDALONE_DIR=${PWD} ROOTSYS=${ROOTSYS}  CXXFLAGS='${RCXX}' lib

JetMETObjects:
	@cvs -d :gserver:cmssw.cvs.cern.ch:/local/reps/CMSSW co -r V03-02-02 -d JetMETObjects CMSSW/CondFormats/JetMETObjects
	patch -p0 < JetMETObjects.patch
	rm -f JetMETObjects/CondFormats; ln -sf ../ JetMETObjects/CondFormats

PUReweighting:
	@cvs -d :gserver:cmssw.cvs.cern.ch:/local/reps/CMSSW co -d PUReweighting CMSSW/PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h
	cd PUReweighting && patch LumiReweightingStandAlone.h ../LumiReweightingStandAlone.patch

#rules
.cc.o:
	$(CC) $(RCXX) -MMD -c -o $@ $<
	@sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
             -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.d

-include $(SRCS:%.cc=%.d) $(OTHERSRCS:%.cc=%.d)

