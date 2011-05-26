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
ROOTLIBS=$(shell root-config --libs) -lMinuit
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


#all files that end up in the Kalibri library:
SRCS=Kalibri.cc GammaJetSel.cc ZJetSel.cc NJetSel.cc TopSel.cc ConfigFile.cc CalibData.cc Parametrization.cc Parameters.cc ControlPlots.cc ControlPlotsProfile.cc ControlPlotsFunction.cc ControlPlotsConfig.cc ControlPlotsResolution.cc ToyMC.cc EventReader.cc PhotonJetReader.cc DiJetReader.cc ThreadedDiJetReader.cc TriJetReader.cc ZJetReader.cc TopReader.cc ParameterLimitsReader.cc EventProcessor.cc EventWeightProcessor.cc Jet.cc JetTruthEvent.cc JetWithTowers.cc TwoJetsInvMassEvent.cc TwoJetsPtBalanceEvent.cc JetWithTracks.cc DiJetResolutionEvent.cc JetConstraintEvent.cc CorFactorsFactory.cc JetBin.cc Binning.cc Function.cc ResolutionParametrization.cc ResolutionFunction.cc ParameterLimit.cc JetWidthEvent.cc EventBinning.cc DiJetEventWeighting.cc 

OBJS = $(SRCS:.cc=.o)

.PHONY: clean

all: dirs lib bin

dirs:
	@mkdir -p bin
	@mkdir -p lib
	@mkdir -p tmp

clean:
	@rm -rf ti_files
	@rm -f *~
	@rm -f *.o 
	@rm -f *#
	@rm -f *.d
	@rm -f .nfs*
	@rm -f *.bkp
	@rm -f junk
	@rm -f caliber
	@rm -f libKalibri.so
	@rm -rf tmp
	@rm -rf lib
	@rm -rf bin
	@rm -f *.cfi
	@rm -f fort.*
	@rm -f .#*

lib: dirs lib/libKalibri.so

bin: dirs junk caliber

plugins: dirs lib/libJetMETCor.so

lbfgs.o: lbfgs.F
	$(F77) $(RCXX) -fno-automatic -fno-backslash -O -c lbfgs.F


lib/libKalibri.so: $(OBJS) lbfgs.o
	$(LD) $(RCXX) -shared $^ $(RLXX) -o lib/libKalibri.so
	@echo '-> Kalibri library created.'


junk: $(OBJS) lbfgs.o caliber.o
	$(LD) $^ $(RLXX) -o bin/junk
	@ln -s -f bin/junk
	@echo '-> static Kalibri executable created.'

caliber: caliber.o lib/libKalibri.so
	$(LD) caliber.o $(RLXX) -Llib -lKalibri -o bin/caliber
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

lib/libJetMETObjects.so: dirs JetMETObjects
	@env STANDALONE_DIR=${PWD} ROOTSYS=${ROOTSYS}  CXXFLAGS='${RCXX} -I.'  /bin/sh -c 'make -e -C JetMETObjects'

JetMETObjects:
	@cvs -d :pserver:anonymous@cmscvs.cern.ch:/cvs_server/repositories/CMSSW co -r V03-02-02 -d JetMETObjects CMSSW/CondFormats/JetMETObjects
	patch -p0 < JetMETObjects.patch



.cc.o:
	$(CC) $(RCXX) -MMD -c -o $@ $<
	@sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
             -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.d

-include $(SRCS:%.cc=%.d) ControlPlotsComparison.d caliber.d compareControlPlots.d JetMETCorFactorsFactory.d

