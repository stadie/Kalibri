C=g++
LD=g++
F77=g77
F77LDFLAGS=-lg2c
F77EXE=$(shell which $(F77) )
ifeq ($(F77EXE),)
  F77=gfortran
  F77LDFLAGS=-lgfortran
endif

#O2 for optimization, g for debugging, pg for profiling
SPECIALFLAGS= -g #-g -Wall -pg#-O2
ROOTAUXCFLAGS=$(shell root-config --auxcflags)
ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs) -lMinuit
#-I. -I./include -I$(SRT_PUBLIC_CONTEXT)/include 
CFLAGS = $(SPECIALFLAGS) -Wall $(ROOTAUXCFLAGS)
#-L../../lib/$(SRT_SUBDIR)/
LFLAGS = $(SPECIALFLAGS) -lz $(F77LDFLAGS)

RCXX=$(SPECIALFLAGS) -Wno-deprecated -Wall $(ROOTCFLAGS)
RLXX=$(LFLAGS) $(ROOTLIBS) -lboost_thread -lpthread  #-lrt -lpthread # -lposix4


SRC=caliber.cc GammaJetSel.cc ZJetSel.cc TrackClusterSel.cc NJetSel.cc TopSel.cc ConfigFile.cc CalibData.cc Parameters.cc ControlPlots.cc ToyMC.cc EventReader.cc PhotonJetReader.cc DiJetReader.cc TriJetReader.cc ZJetReader.cc TopReader.cc ParameterLimitsReader.cc TowerConstraintsReader.cc TrackClusterReader.cc EventProcessor.cc Jet.cc JetTruthEvent.cc JetWithTowers.cc TwoJetsInvMassEvent.cc JetWithTracks.cc

%.o: %.cc
		$(C) $(RCXX) -c $<

all: runjunk

lbfgs.o: lbfgs.F
	$(F77) -fno-automatic -fno-backslash -O -c lbfgs.F

ConfigFile.o: ConfigFile.cc ConfigFile.h
	$(C) $(CFLAGS) -c ConfigFile.cc

GammaJetSel.o: GammaJetSel.cc GammaJetSel.h
	$(C) $(RCXX) -c GammaJetSel.cc

ZJetSel.o: ZJetSel.cc ZJetSel.h
	$(C) $(RCXX) -c ZJetSel.cc

TopSel.o: TopSel.cc TopSel.h
	$(C) $(RCXX) -c TopSel.cc

TrackTowerSel.o: TrackTowerSel.cc TrackTowerSel.h
	$(C) $(RCXX) -c TrackTowerSel.cc

TrackClusterSel.o: TrackClusterSel.cc TrackClusterSel.h
	$(C) $(RCXX) -c TrackClusterSel.cc

NJetSel.o: NJetSel.cc NJetSel.h
	$(C) $(RCXX) -c NJetSel.cc

CalibData.o: CalibData.cc CalibData.h Parametrization.h Parameters.h
	$(C) $(RCXX) -c CalibData.cc

Parameters.o: Parameters.cc Parameters.h Parametrization.h Function.h ConfigFile.h
	$(C) $(RCXX) -c Parameters.cc

ControlPlots.o: ControlPlots.cc ControlPlots.h CalibData.h CalibMath.h ConfigFile.h TwoJetsInvMassEvent.h
	$(C) $(RCXX) -c ControlPlots.cc

EventReader.o: EventReader.h EventReader.cc Parameters.h ConfigFile.h 
	$(C) $(RCXX) -c EventReader.cc

PhotonJetReader.o: EventReader.h PhotonJetReader.h PhotonJetReader.cc  GammaJetSel.h ToyMC.h Parameters.h ConfigFile.h Jet.h JetTruthEvent.h JetWithTowers.h Function.h
	$(C) $(RCXX) -c PhotonJetReader.cc

DiJetReader.o: EventReader.h DiJetReader.h DiJetReader.cc NJetSel.h ToyMC.h Parameters.h ConfigFile.h Jet.h JetTruthEvent.h JetWithTowers.h Function.h
	$(C) $(RCXX) -c DiJetReader.cc

TriJetReader.o: EventReader.h TriJetReader.h TriJetReader.cc NJetSel.h Parameters.h ConfigFile.h
	$(C) $(RCXX) -c TriJetReader.cc

ZJetReader.o: EventReader.h ZJetReader.h ZJetReader.cc ZJetSel.h Parameters.h ConfigFile.h Jet.h JetTruthEvent.h JetWithTowers.h Function.h
	$(C) $(RCXX) -c ZJetReader.cc

TopReader.o: EventReader.h TopReader.h TopReader.cc TopSel.h Parameters.h ConfigFile.h
	$(C) $(RCXX) -c TopReader.cc

ParameterLimitsReader.o: EventReader.h ParameterLimitsReader.h ParameterLimitsReader.cc Parameters.h ConfigFile.h
	$(C) $(RCXX) -c ParameterLimitsReader.cc

TowerConstraintsReader.o:  EventReader.h TowerConstraintsReader.h TowerConstraintsReader.cc Parameters.h ConfigFile.h
	$(C) $(RCXX) -c TowerConstraintsReader.cc

TrackClusterReader.o: EventReader.h TrackClusterReader.h TrackClusterReader.cc TrackClusterSel.h Parameters.h ConfigFile.h
	$(C) $(RCXX) -c TrackClusterReader.cc

EventProcessor.o: ConfigFile.h Parameters.h EventProcessor.h EventProcessor.cc
	$(C) $(RCXX) -c EventProcessor.cc

Jet.o: CalibData.h Jet.h Jet.cc Parametrization.h Function.h
	$(C) $(RCXX) -c Jet.cc	

JetTruthEvent.o: CalibData.h Jet.h JetTruthEvent.h JetTruthEvent.cc Jet.h
	$(C) $(CFLAGS) -c JetTruthEvent.cc

TwoJetsInvMassEvent.o: CalibData.h Jet.h TwoJetsInvMassEvent.h TwoJetsInvMassEvent.cc Jet.h 
	$(C) $(RCXX) -c TwoJetsInvMassEvent.cc

JetWithTowers.o: CalibData.h Jet.h JetWithTowers.h Function.h JetWithTowers.cc Parametrization.h
	$(C) $(RCXX) -c JetWithTowers.cc

JetWithTracks.o: CalibData.h Jet.h JetWithTracks.h Function.h JetWithTracks.cc Parametrization.h
	$(C) $(RCXX) -c JetWithTracks.cc

caliber.o: caliber.cc caliber.h CalibMath.h external.h ConfigFile.h CalibData.h Parameters.h ControlPlots.h EventReader.h DiJetReader.h TriJetReader.h ZJetReader.h TopReader.h ParameterLimitsReader.h TowerConstraintsReader.h TrackClusterReader.h EventProcessor.h Jet.h TwoJetsInvMassEvent.h
	$(C) $(RCXX)  -I/usr/include/boost -c caliber.cc 

runjunk: $(SRC:.cc=.o) lbfgs.o
	$(LD) $(SRC:.cc=.o) lbfgs.o $(RLXX) $(JCORR) -o junk
	@echo '-> Calibration executable created.'

clean:
	@rm -rf ti_files
	@rm -f *~
	@rm -f *.o 
	@rm -f *#
	@rm -f .nfs*
	@rm -f *.bkp
	@rm -f junk
	@rm -f *.ps
	@rm -f *.eps
	@rm -f *.cfi
	@rm -f fort.*
	@rm -f .#*



ToyMC.o: ToyMC.h ToyMC.cc ConfigFile.h
	$(C) $(RCXX) -c ToyMC.cc

toy:	ToyMC.o toy.o ConfigFile.o
	$(LD) ToyMC.o toy.o ConfigFile.o $(RLXX) -o toy
	@echo '-> toy MC executable created.'

ControlPlotsComparison.o: ControlPlotsComparison.cc ControlPlotsComparison.h
	$(C) $(RCXX) -c ControlPlotsComparison.cc

comp: 	ControlPlotsComparison.o compareControlPlots.o
	$(LD) ControlPlotsComparison.o compareControlPlots.o $(RLXX) -o compControlPlots
	@echo '-> Comparison executable created. Type "compControlPlots" to compare control plots.'
