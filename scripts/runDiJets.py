#! /usr/bin/python

import os

def writeCfg(filename):
    config="""
#
# Configuration File for the Calibration Program
# Hamburg, 2007/08/15
#

#---------------------------------------------------------------------------------
#   Mode
#---------------------------------------------------------------------------------
# Running mode: Calibration (0), Jet smearing (1)
Mode = 0


#---------------------------------------------------------------------------------
#Fit method: LVMini=1, use existing calibration=2, write start values=3
Fit method = 2
#Parametrization
Parametrization Class = MeanWidthParametrization

#Error Parametrization
tower error parametrization = const
#jet error parametrization   = jet et
#jet error parametrization   = const

start values = 1.0
#jet start values        = 1. 0. 0. 0. 0. 0.
jet start values = 2.22197 -2.43273 0.161247 -1.8384 -1.12056 3.76558 -1.28309 -1.21227 4.97975 -1.06063 -3.92938 3.9911 -1.81438 1.18677 -7.45932 -0.710319 4.92905 -1.52939 -5.97873 -8.6159
global jet start values =  0.981669 4.16531 3.03254 2.0096

# Scaling of residuals: 0 - none, 1 - with Cauchy-Function, 2 - with Huber-Function
Residual Scaling Scheme    = 1111 #   221 : default
Outlier Cut on Chi2        = 1000.0 # Applied before each iteration with no scaling

BFGS derivative step     = 1e-05
BFGS mvec                = 20
BFGS niter               = 1000
BFGS eps                 = 1e-03
BFGS 1st wolfe parameter = 1.E-04
BFGS 2nd wolfe parameter = 0.9
BFGS print derivatives   = false


#---------------------------------------------------------------------------------
#   Geometry / Binning for fitting
#---------------------------------------------------------------------------------
maximum eta twr used  = 82
granularity in eta    = 1   # allowed values are: 1,3,5,11,21,41 (*2)
granularity in phi    = 1   # allowed values are: 1,2,3,6,9,18,36,72
symmetry in eta       = true

jet granularity in eta = 4 #  allowed values are: 1,3,5,11,21,41 (*2)
jet granularity in phi = 1 #   1 : default

track granularity in eta = 1
track granularity in phi = 1

jet binning variables = eta;pt
jet binning eta bins = -5.191 -4.889 -4.716 -4.538 -4.363 -4.191 -4.013 -3.839 -3.664 -3.489 -3.314 -3.139 -2.964 -2.853 -2.650 -2.500 -2.322 -2.172 -2.043 -1.930 -1.830 -1.740 -1.653 -1.566 -1.479 -1.392 -1.305 -1.218 -1.131 -1.044 -0.957 -0.879 -0.783 -0.696 -0.609 -0.522 -0.435 -0.348 -0.261 -0.174 -0.087 0.000 0.087 0.174 0.261 0.348 0.435 0.522 0.609 0.696 0.783 0.879 0.957 1.044 1.131 1.218 1.305 1.392 1.479 1.566 1.653 1.740 1.830 1.930 2.043 2.172 2.322 2.500 2.650 2.853 2.964 3.139 3.314 3.489 3.664 3.839 4.013 4.191 4.363 4.538 4.716 4.889 5.191
jet binning pt bins = 5 10 12 15 18 22 26 30 35 40 45 51 57 64 72 80 90 105 120 135 150 175 200 250 300 350 400 500 650 800 1000 1500 5000
jet binning emf bins = 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
#---------------------------------------------------------------------------------
#   Analyses (GammaJet, TrackTower, ... )
#---------------------------------------------------------------------------------
Et cut on jet              = 0.0
Et cut on gamma            = 0.0
Et cut on Z                = 20.0
Et cut on tower            = 0.0
Et cut on cluster          = 0.0
Et cut on track            = 0.0
Et cut on n+1 Jet          = 9999.0
Min Delta Phi              = 2.7
Eta cut on jet             = 5.2
Relative Rest Jet Cut      = 0.2      #NonLeadingJetsEt / PhotonEt
Min had fraction           = -0.05    #Default: 0.07
Max had fraction           = 1.05    #Default: 0.95
Et genJet min              = 0.0
Et genJet max              = 7000.0
DeltaR cut on jet matching = 9999
Max cut on relative n+1 Jet Et = 0.3
Max cut on relative Soft Jet Et = 0.1

#---------------------------------------------------------------------------------
#   Input / Output
#---------------------------------------------------------------------------------

#jet correction source = JetMETCor
#jet correction name   = Spring10_AK5TRK
#jet correction source = JetMETCor
#jet correction name   = Spring10_AK5CaloData
Default Tree Name      = CalibTree

# List of input files:
Gamma-Jet tree         = GammaJetTree
Z-Jet tree             = ZJetTree
Track-Tower tree       = TrackTowerTree
Track-Cluster tree     = TrackClusterTree
Di-Jet tree            = DiJetTree
Di-Jet Control1 tree   = DiJetTree
Tri-Jet tree           = TriJetTree
Top tree               = TopTree

#Di-Jet input file = input/dijetlistspring10Track
#Di-Jet input file =  input/toy_dijet_const_5-500_uniform.root

#Top input file = /scratch/current/cms/user/stadie/Top_Madgraph.root
#Z-Jet input file = /scratch/current/cms/user/stadie/ZJet_Track_230_300_rereco.root; /scratch/current/cms/user/stadie/ZJet_Track_300_INF_rereco.root; 

# -1: use all events, 1000: run over 1000 events
use Gamma-Jet events     = 0
use Track-Tower events   = 0
use Track-Cluster events = 0
#use Di-Jet events        = -1
use Tri-Jet events       = 0
use Z-Jet events         = 0
use Top events           = 0

Gamma-Jet data class     = 1
Z-Jet data class     = 3
#Di-Jet data class    = 21
Top data class       = 1

correct jets L2L3 = false
#Di-Jet prescale = 1000
#-----------------------------------------------------------------
#  Control plots
#-----------------------------------------------------------------
#  General parameters
create plots                     = true
plots output directory           = plots
#plots format                      = pdf

# JetTruthEvent plots
create JetTruthEvent plots    =  false
create TwoJetsPtBalanceEvent plots = true

TwoJetsPtBalanceEvent plots names =  AsymmetryVsEta;AsymmetryVsPt25;AsymmetryVsPt20;AsymmetryVsPt15;AsymmetryVsPt10;AsymmetryVsPt05

AsymmetryVsEta x variable        =  Eta
AsymmetryVsEta x edges           =  26 -5.2 5.2
AsymmetryVsEta y variable        =  Asymmetry
AsymmetryVsEta y edges           =  51 -1 1 -0.5 0.8 -0.5 0.8 0.7 1.3 0.7 1.3
AsymmetryVsEta bin variable      =  MeanPt
AsymmetryVsEta bin edges         =  20 30 50 80 120 200 360 500 900 7000
AsymmetryVsEta cut variable      =  ThirdJetFraction
AsymmetryVsEta cut edges         =  0.0 0.20
AsymmetryVsEta correction types  =  L2L3; L2L3Res
AsymmetryVsEta 1 correction types  =  L2L3
AsymmetryVsEta profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AsymmetryVsEta legend label      =  L2L3:CMS default
AsymmetryVsEta input samples     =  0:data;1:MC

AsymmetryVsPt25 x variable        =  MeanPt; log
AsymmetryVsPt25 x edges           =  30 20 2000
AsymmetryVsPt25 y variable        =  Asymmetry
AsymmetryVsPt25 y edges           =  51 -1 1 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
AsymmetryVsPt25 bin variable      =  Eta
#AsymmetryVsPt25 bin edges        =  0.0 1.3 2.6 3.0 5.2
AsymmetryVsPt25 bin edges         =  -5 -4.5 -4.2 -3.9 -3.6 -3.3 -3.0 -2.7 -2.4 -2.1 -1.8 -1.5 -1.3 -1.1 -0.9 -0.6 -0.3 0.0 0.3 0.6 0.9 1.1 1.3 1.5 1.8 2.1 2.4 2.7 3.0 3.3 3.6 3.9 4.2 4.5 5 
AsymmetryVsPt25 cut variable      =  ThirdJetFraction
AsymmetryVsPt25 cut edges         =  0.0 0.25
AsymmetryVsPt25 correction types  =  L2L3; L2L3Res
AsymmetryVsPt25 1 correction types  =  L2L3
AsymmetryVsPt25 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AsymmetryVsPt25 legend label      =  L2L3:CMS default
AsymmetryVsPt25 input samples     =  0:data;1:MC

AsymmetryVsPt20 x variable        =  MeanPt; log
AsymmetryVsPt20 x edges           =  30 20 2000
AsymmetryVsPt20 y variable        =  Asymmetry
AsymmetryVsPt20 y edges           =  51 -1 1 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
AsymmetryVsPt20 bin variable      =  Eta
#AsymmetryVsPt20 bin edges        =  0.0 1.3 2.6 3.0 5.2
AsymmetryVsPt20 bin edges         =  -5 -4.5 -4.2 -3.9 -3.6 -3.3 -3.0 -2.7 -2.4 -2.1 -1.8 -1.5 -1.3 -1.1 -0.9 -0.6 -0.3 0.0 0.3 0.6 0.9 1.1 1.3 1.5 1.8 2.1 2.4 2.7 3.0 3.3 3.6 3.9 4.2 4.5 5 
AsymmetryVsPt20 cut variable      =  ThirdJetFraction
AsymmetryVsPt20 cut edges         =  0.0 0.2
AsymmetryVsPt20 correction types  =  L2L3; L2L3Res
AsymmetryVsPt20 1 correction types  =  L2L3
AsymmetryVsPt20 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AsymmetryVsPt20 legend label      =  L2L3:CMS default
AsymmetryVsPt20 input samples     =  0:data;1:MC

AsymmetryVsPt15 x variable        =  MeanPt; log
AsymmetryVsPt15 x edges           =  30 20 2000
AsymmetryVsPt15 y variable        =  Asymmetry
AsymmetryVsPt15 y edges           =  51 -1 1 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
AsymmetryVsPt15 bin variable      =  Eta
#AsymmetryVsPt15 bin edges        =  0.0 1.3 2.6 3.0 5.2
AsymmetryVsPt15 bin edges         =  -5 -4.5 -4.2 -3.9 -3.6 -3.3 -3.0 -2.7 -2.4 -2.1 -1.8 -1.5 -1.3 -1.1 -0.9 -0.6 -0.3 0.0 0.3 0.6 0.9 1.1 1.3 1.5 1.8 2.1 2.4 2.7 3.0 3.3 3.6 3.9 4.2 4.5 5 
AsymmetryVsPt15 cut variable      =  ThirdJetFraction
AsymmetryVsPt15 cut edges         =  0.0 0.15
AsymmetryVsPt15 correction types  =  L2L3; L2L3Res
AsymmetryVsPt15 1 correction types  =  L2L3
AsymmetryVsPt15 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AsymmetryVsPt15 legend label      =  L2L3:CMS default
AsymmetryVsPt15 input samples     =  0:data;1:MC

AsymmetryVsPt10 x variable        =  MeanPt; log
AsymmetryVsPt10 x edges           =  30 20 2000
AsymmetryVsPt10 y variable        =  Asymmetry
AsymmetryVsPt10 y edges           =  51 -1 1 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
AsymmetryVsPt10 bin variable      =  Eta
#AsymmetryVsPt10 bin edges        =  0.0 1.3 2.6 3.0 5.2
AsymmetryVsPt10 bin edges         =  -5 -4.5 -4.2 -3.9 -3.6 -3.3 -3.0 -2.7 -2.4 -2.1 -1.8 -1.5 -1.3 -1.1 -0.9 -0.6 -0.3 0.0 0.3 0.6 0.9 1.1 1.3 1.5 1.8 2.1 2.4 2.7 3.0 3.3 3.6 3.9 4.2 4.5 5 
AsymmetryVsPt10 cut variable      =  ThirdJetFraction
AsymmetryVsPt10 cut edges         =  0.0 0.1
AsymmetryVsPt10 correction types  =  L2L3; L2L3Res
AsymmetryVsPt10 1 correction types  =  L2L3
AsymmetryVsPt10 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AsymmetryVsPt10 legend label      =  L2L3:CMS default
AsymmetryVsPt10 input samples     =  0:data;1:MC

AsymmetryVsPt05 x variable        =  MeanPt; log
AsymmetryVsPt05 x edges           =  30 20 2000
AsymmetryVsPt05 y variable        =  Asymmetry
AsymmetryVsPt05 y edges           =  51 -1 1 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
AsymmetryVsPt05 bin variable      =  Eta
#AsymmetryVsPt05 bin edges        =  0.0 1.3 2.6 3.0 5.2
AsymmetryVsPt05 bin edges         =  -5 -4.5 -4.2 -3.9 -3.6 -3.3 -3.0 -2.7 -2.4 -2.1 -1.8 -1.5 -1.3 -1.1 -0.9 -0.6 -0.3 0.0 0.3 0.6 0.9 1.1 1.3 1.5 1.8 2.1 2.4 2.7 3.0 3.3 3.6 3.9 4.2 4.5 5 
AsymmetryVsPt05 cut variable      =  ThirdJetFraction
AsymmetryVsPt05 cut edges         =  0.0 0.05
AsymmetryVsPt05 correction types  =  L2L3; L2L3Res
AsymmetryVsPt05 1 correction types  =  L2L3
AsymmetryVsPt05 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AsymmetryVsPt05 legend label      =  L2L3:CMS default
AsymmetryVsPt05 input samples     =  0:data;1:MC


AsymmetryVsEtaEta x variable       =  momentEtaEta
AsymmetryVsEtaEta x edges          =  30 0 0.3
AsymmetryVsEtaEta y variable        =  Asymmetry
AsymmetryVsEtaEta y edges           =  51 -1 1 -0.5 0.5 -0.5 0.5 0.5 2.0 0.5 2.5
AsymmetryVsEtaEta bin variable      =  AbsEta
AsymmetryVsEtaEta bin edges         =  0.0 1.3 2.6 3.0 5.2

AsymmetryVsEtaEta correction types  =  Uncorrected; L2L3
AsymmetryVsEtaEta profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AsymmetryVsEtaEta legend label      =  L2L3:CMS default
AsymmetryVsEtaEta input samples     =  0:data;1:MC

AsymmetryVsEmf x variable       =  EMF
AsymmetryVsEmf x edges          =  20 0 1.0
AsymmetryVsEmf y variable        =  Asymmetry
AsymmetryVsEmf y edges           =  51 -1 1 -0.5 0.5 -0.5 0.5 0.0 2.5 0.0 2.5
AsymmetryVsEmf bin variable      =  AbsEta
AsymmetryVsEmf bin edges         =  0.0 1.3 2.6 3.0 5.2
AsymmetryVsEmf correction types  =  Uncorrected; L2L3
AsymmetryVsEmf profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AsymmetryVsEmf legend label      =  L2L3:CMS default
AsymmetryVsEmf input samples     =  0:data;1:MC

AsymmetryVsMeanMoment x variable       =  meanMoment
AsymmetryVsMeanMoment x edges          =  30 0 0.3
AsymmetryVsMeanMoment y variable        =  Asymmetry
AsymmetryVsMeanMoment y edges           =  51 -1 1 -0.5 0.5 -0.5 0.5 0.0 2.5 0.0 2.5
AsymmetryVsMeanMoment bin variable      =  AbsEta
AsymmetryVsMeanMoment bin edges         =  0 1.3 2.6 3.0 5.2
AsymmetryVsMeanMoment cut variable      =  MeanPt
AsymmetryVsMeanMoment cut edges         =  20 7000
AsymmetryVsMeanMoment correction types  =  Uncorrected; L2L3
AsymmetryVsMeanMoment profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AsymmetryVsMeanMoment legend label      =  L2L3:CMS default
AsymmetryVsMeanMoment input samples     =  0:data;1:MC

AsymmetryVsMeanMomentMeanPt x variable       =  meanMoment
AsymmetryVsMeanMomentMeanPt x edges          =  30 0 0.3
AsymmetryVsMeanMomentMeanPt y variable        =  Asymmetry
AsymmetryVsMeanMomentMeanPt y edges           =  51 -1 1 -0.5 0.5 -0.5 0.5 0.0 2.5 0.0 2.5
AsymmetryVsMeanMomentMeanPt bin variable      =  MeanPt
AsymmetryVsMeanMomentMeanPt bin edges         =  5 10 20 30 50 80 120 200 360 7000
AsymmetryVsMeanMomentMeanPt cut variable      =  Eta
AsymmetryVsMeanMomentMeanPt cut edges         =  -1.3 1.3
AsymmetryVsMeanMomentMeanPt correction types  =  Uncorrected; L2L3
AsymmetryVsMeanMomentMeanPt profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AsymmetryVsMeanMomentMeanPt legend label      =  L2L3:CMS default
AsymmetryVsMeanMomentMeanPt input samples     =  0:data;1:MC

AsymmetryVsMeanMomentMeanPt2 x variable       =  meanMoment
AsymmetryVsMeanMomentMeanPt2 x edges          =  30 0 0.3
AsymmetryVsMeanMomentMeanPt2 y variable        =  Asymmetry
AsymmetryVsMeanMomentMeanPt2 y edges           =  51 -1 1 -0.5 0.5 -0.5 0.5  0.0 2.5 0.0 2.5
AsymmetryVsMeanMomentMeanPt2 bin variable      =  MeanPt
AsymmetryVsMeanMomentMeanPt2 bin edges         =  5 10 20 30 50 80 120 200 360 7000
AsymmetryVsMeanMomentMeanPt2 cut variable      =  AbsEta
AsymmetryVsMeanMomentMeanPt2 cut edges         =  1.3 2.6
AsymmetryVsMeanMomentMeanPt2 correction types  =  Uncorrected; L2L3
AsymmetryVsMeanMomentMeanPt2 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AsymmetryVsMeanMomentMeanPt2 legend label      =  L2L3:CMS default
AsymmetryVsMeanMomentMeanPt2 input samples     =  0:data;1:MC

"""
    fcfg = open(filename, "w")
    fcfg.write(config)
    fcfg.write("Di-Jet input file = dijetlist\n")
    fcfg.write("Di-Jet Control1 input file = mcdijetlist\n")
    fcfg.write("Output file       = "+output+"\n");
    fcfg.write("Number of Threads = "+str(nthreads)+"\n")
    fcfg.write("Number of IO Threads = 5\n")
    if(useconstraint):
        fcfg.write("jet constraints =  5.0 10.0 0.0 1.2 1 10.0 15.0 0.0 1.2 1 15.0 20.0 0 1.2 1 20.0 25.0 0 1.2 1 25.0 30.0 0 1.2 1 30.0 40.0 0 1.2 1 40.0 50.0 0 1.2  1 50.0 60.0 0 1.2 1 60.0 70.0 0 1.2  1 70.0 80 0 1.2 1  80.0 90.0 0 1.2 1 90.0 100.0 0 1.2 1 100.0 120. 0 1.2 1 120 150 0 1.2 1 150 200 0 1.2 1 200 280 0 1.2 1 280 350 0 1.2 1 350 500 0 1.2 1 500 800 0 1.2 1 800 1400 0 1.2 1 1400 7000 0 1.2 1\n")
    else:
        fcfg.write("fixed jet parameters = 7 1\n")
        
    if(binned):
        fcfg.write("Di-Jet data class    = 1\n")
        fcfg.write("jet error parametrization   = const\n")
    else:
        fcfg.write("Di-Jet data class    = 1\n")
        fcfg.write("jet error parametrization   = jet et\n")

    fcfg.write("use Di-Jet events = "+str(nevents)+"\n")
    fcfg.write("use Di-Jet Control1 events = "+str(nevents)+"\n")
    fcfg.write("Di-Jet trigger names = HLT_DiJetAve15U;HLT_DiJetAve30U;HLT_DiJetAve50U;HLT_DiJetAve70U;HLT_DiJetAve100U;HLT_DiJetAve140U\n")
    if(jettype == "ak5Calo"):
        fcfg.write("Di-Jet trigger thresholds = 38 59 86 111 147 196\n")
    if(jettype == "ak5PF"):
        fcfg.write("Di-Jet trigger thresholds = 43 70 100 127 168 214\n")
    if(input != ""):
        fcfg.write("input calibration = Kalibri; "+input+"\n");

    fcfg.close()
    return


#main program starts here!!!
#change these variables to steer the fit
jettype = "ak5Calo"
#datadir = "/scratch/hh/current/cms/user/stadie/QCDFlat_Pt15to3000Spring10-START3X_V26_S09-v1C"
#datadir = "/scratch/hh/current/cms/user/stadie/JetMET_Run2010A-PromptReco-v4_DCSONLY"
#datadir = "/scratch/hh/current/cms/user/stadie/Spring10QCDDiJetPt*"
datadir = "/scratch/hh/current/cms/user/stadie/DiJetNov4ReReco_v1A"
datadirmc = "/scratch/hh/current/cms/user/stadie/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Fall10-START38_V12-v1_matthias_merged"
nthreads = 4
nevents =  -1
dirname = "/afs/naf.desy.de/group/cms/scratch/stadie/dijetsFall10AK5Calo_weighted"
useconstraint = False
batch = False
doBinnedFit = False
doUnbinnedFit = True


#write configs and run the fit
print "producing validation plots for "+jettype+" using "+datadir+" and "+datadirmc+" as MC sample\n"


if os.path.exists(dirname):
    os.system("rm -rf "+dirname+"/*")
else:
    os.system("mkdir "+dirname)
    
os.system("cp junk "+dirname+"/.")
os.system("ln -s $PWD/kalibriLogoSmall.gif "+dirname+"/")
os.system("ln -s $PWD/L4JWfit.txt "+dirname+"/")
os.system("ln -s $PWD/lib "+dirname+"/")
os.system("ln -s $PWD/JetMETObjects "+dirname+"/")
os.system("ls "+datadir+"/*"+jettype+"*.root > "+dirname+"/dijetlist");
os.system("ls "+datadirmc+"/*"+jettype+"*.root > "+dirname+"/mcdijetlist");

binned= True
output="Kalibri.txt"
input ="Kalibri.txt"
writeCfg(dirname+"/L2L3.cfg")
binned = False
output="Kalibri2.txt"
input ="Kalibri.txt"
writeCfg(dirname+"/L2L3b.cfg")

if batch:
    fjob = open(dirname+"/fitL2L3.sh", "w")
    fjob.write("#! /bin/sh\n")
    fjob.write("#\n")
    fjob.write("#$ -V\n")
    fjob.write("#$ -pe  multicore "+str(nthreads)+"\n")
    fjob.write("#$ -R Y\n")
    fjob.write("#$ -l h_cpu=23:00:00\n")
    fjob.write("#$ -l h_vmem=5000M\n")
    fjob.write("cd "+os.getcwd()+"/"+dirname+"\n")
    fjob.write("date\n")
    if doBinnedFit:
        fjob.write("./junk L2L3.cfg > $TMPDIR/L2L3.log\n")
        fjob.write("date\n")
        fjob.write("mv $TMPDIR/L2L3.log .\n")

    if doUnbinnedFit:
        fjob.write("./junk L2L3b.cfg > $TMPDIR/L2L3b.log\n")
        fjob.write("date\n")
        fjob.write("mv $TMPDIR/L2L3b.log .\n")

    fjob.close()
    qsubcmd = "qsub "+dirname+"/fitL2L3.sh"
    print "running "+qsubcmd
    os.system(qsubcmd)
else:
    if doBinnedFit:
        kalibricmd = "cd "+dirname+"; ./junk L2L3.cfg; cd -";
        print "running "+kalibricmd
        os.system(kalibricmd)

    if doUnbinnedFit:
        kalibricmd = "cd "+dirname+"; ./junk L2L3b.cfg; cd -";
        print "running "+kalibricmd
        os.system(kalibricmd)


