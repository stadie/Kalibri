C=g++
LD=g++
F77=g77
#O2 for optimization, g for debugging
SPECIALFLAGS=-O3 #-g -Wall -pg#-O2
ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs) -lMinuit

CFLAGS = $(SPECIALFLAGS) -I. -I./include -I$(SRT_PUBLIC_CONTEXT)/include -I$(ROOTSYS)/include -Wno-deprecated -Wall
LFLAGS = $(SPECIALFLAGS) -L../../lib/$(SRT_SUBDIR)/ -lz -lg2c


RCXX=$(CFLAGS) $(ROOTCFLAGS) -I/usr/include/boost
RLXX=$(LFLAGS) $(ROOTLIBS)  -I/usr/include/boost -lboost_thread -lpthread  #-lrt -lpthread # -lposix4

SRC=caliber.cc GammaJetSel.cc ZJetSel.cc TrackTowerSel.cc TrackClusterSel.cc NJetSel.cc ConfigFile.cc CalibData.cc Parameters.cc ControlPlots.cc

%.o: %.cc
		$(C) $(RCXX) -c $<

all: runjunk

lbfgs.o: lbfgs.F
		$(F77) -fno-automatic -fno-backslash -O -c lbfgs.F

ConfigFile.o: ConfigFile.cc ConfigFile.h
		$(C) $(RCXX) -c ConfigFile.cc

GammaJetSel.o: GammaJetSel.cc GammaJetSel.h
		$(C) $(RCXX) -c GammaJetSel.cc

ZJetSel.o: ZJetSel.cc ZJetSel.h
		$(C) $(RCXX) -c ZJetSel.cc

TrackTowerSel.o: TrackTowerSel.cc TrackTowerSel.h
		$(C) $(RCXX) -c TrackTowerSel.cc

TrackClusterSel.o: TrackClusterSel.cc TrackClusterSel.h
		$(C) $(RCXX) -c TrackClusterSel.cc

NJetSel.o: NJetSel.cc NJetSel.h
		$(C) $(RCXX) -c NJetSel.cc

CalibData.o: CalibData.cc CalibData.h
		$(C) $(RCXX) -c CalibData.cc

Parameters.o: Parameters.cc Parameters.h Parametrization.h
		$(C) $(RCXX) -c Parameters.cc

ControlPlots.o: ControlPlots.cc ControlPlots.h
		$(C) $(RCXX) -c ControlPlots.cc


caliber.o: caliber.cc caliber.h CalibMath.h external.h GammaJetSel.h TrackTowerSel.h ZJetSel.h NJetSel.h ConfigFile.h CalibData.h Parameters.h ControlPlots.h
		$(C) $(RCXX) -c caliber.cc 

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
