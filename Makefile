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
else ifneq ($(wildcard /usr/lib64/libboost_thread-mt.so),)
 BOOSTLINKFLAGS=-lboost_thread-mt 
#-lpthread
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

.PHONY: clean bins libs plugins all PUreweighting lbfgs

all: include libs plugins bins

bin:
	@mkdir -p bin
lib:
	@mkdir -p lib
tmp:
	@mkdir -p tmp

include:
	@mkdir -p include

clean:
	@rm -rf ti_files tmp lib bin share  
	@rm -f *~ *.o *# *.d *.bkp junk caliber libKalibri.so *.cfi fort.* .#* 


libs: include lib 

bins: plugins bin bin/junk bin/caliber

plugins: PUreweighting lbfgs lib lib/libJetMETCor.so 

PUreweighting: PUReweighting/LumiReweightingStandAlone.h


lbfgs:  include/lbfgs.h lib/liblbfgs.so

lbfgs.o: lbfgs.F
	$(F77) $(RCXX) -fno-automatic -fno-backslash -O -c lbfgs.F

lib/libKalibri.so:  include/lbfgs.h $(OBJS) lbfgs.o
	$(LD) $(RCXX) -shared $^ $(RLXX) -o lib/libKalibri.so
	@echo '-> Kalibri library created.'

liblbfgs/configure: liblbfgs/configure.in
	@cd liblbfgs &&  autoconf

liblbfgs/Makefile: liblbfgs/configure
	@cd liblbfgs && ./configure --enable-sse2 --prefix=$(KALIBRIDIR)
	@echo '-> liblbfgs configured.'

liblbfgs/configure.in:
	@git submodule init
	@git submodule update
	@cd liblbfgs && libtoolize --force
	@cd liblbfgs && aclocal
	@cd liblbfgs && autoheader
	@cd liblbfgs && automake --force-missing --add-missing


include/lbfgs.h lib/liblbfgs.so lib/liblbfgs.a: liblbfgs/Makefile
	@cd liblbfgs && $(MAKE) && $(MAKE) install
	@echo '-> shared library lib/liblbfgs.so created.'


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
	cp -r /afs/cern.ch/cms/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_20/src/CondFormats/JetMETObjects .
	cd JetMETObjects && patch -p0 < ../JetMETObjects.patch
	rm -f JetMETObjects/CondFormats; ln -sf ../ JetMETObjects/CondFormats


PUReweighting/LumiReweightingStandAlone.h:
	mkdir PUReweighting
	cp /afs/cern.ch/cms/slc5_amd64_gcc462/cms/cmssw/CMSSW_5_3_20/src/PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h PUReweighting/. 
	cd PUReweighting && patch LumiReweightingStandAlone.h ../LumiReweightingStandAlone.patch

#rules
.cc.o:  include/lbfgs.h 
	$(CC) $(RCXX) -MMD -c -o $@ $<
	@sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
             -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.d

-include $(SRCS:%.cc=%.d) $(OTHERSRCS:%.cc=%.d)

