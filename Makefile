C=g++
LD=g++
F77=g77
#O2 for optimization, g for debugging
<<<<<<< Makefile
SPECIALFLAGS=-O2 #-g -Wall#-O2
=======
SPECIALFLAGS=-Wall -O4
>>>>>>> 1.4
ROOTCFLAGS=$(shell root-config --cflags)
ROOTLIBS=$(shell root-config --libs) -lMinuit

CFLAGS = $(SPECIALFLAGS) -I. -I./include -I$(SRT_PUBLIC_CONTEXT)/include -I$(ROOTSYS)/include -Wno-deprecated -m32 -Wa,--32
LFLAGS = $(SPECIALFLAGS) -L../../lib/$(SRT_SUBDIR)/ -lz -m32 -lg2c


RCXX=$(CFLAGS) $(ROOTCFLAGS)
RLXX=$(LFLAGS) $(ROOTLIBS)

SRC=caliber.C GammaJetSel.C TrackTowerSel.C TrackClusterSel.C JetJetSel.C ConfigFile.C CalibData.C Parameters.C ControlPlots.C

%.o: %.C
		$(C) $(RCXX) -c $<

all: runjunk

lbfgs.o: lbfgs.F
		$(F77) -fno-automatic -fno-backslash -m32 -O -c lbfgs.F

ConfigFile.o: ConfigFile.C ConfigFile.h
		$(C) $(RCXX) -c ConfigFile.C

GammaJetSel.o: GammaJetSel.C GammaJetSel.h
		$(C) $(RCXX) -c GammaJetSel.C

TrackTowerSel.o: TrackTowerSel.C TrackTowerSel.h
		$(C) $(RCXX) -c TrackTowerSel.C

TrackClusterSel.o: TrackClusterSel.C TrackClusterSel.h
		$(C) $(RCXX) -c TrackClusterSel.C

JetJetSel.o: JetJetSel.C JetJetSel.h
		$(C) $(RCXX) -c JetJetSel.C

CalibData.o: CalibData.C CalibData.h
		$(C) $(RCXX) -c CalibData.C

Parameters.o: Parameters.C Parameters.h
		$(C) $(RCXX) -c Parameters.C

ControlPlots.o: ControlPlots.C ControlPlots.h
		$(C) $(RCXX) -c ControlPlots.C

caliber.o: caliber.C caliber.h CalibMath.h external.h GammaJetSel.h TrackTowerSel.h ConfigFile.h CalibData.h Parameters.h ControlPlots.h
		$(C) $(RCXX) -c caliber.C 

runjunk: $(SRC:.C=.o) lbfgs.o
		$(LD) $(SRC:.C=.o) lbfgs.o $(RLXX) $(JCORR) -o junk
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
#		@rm -f *.cfi
		@rm -f fort.*
		@rm -f .#*

