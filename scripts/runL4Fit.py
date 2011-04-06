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
Fit method = 1
#Parametrization
#Parametrization Class = SimplePhiPhiParametrization
Parametrization Class = PhiPhiParametrization
#Error Parametrization
tower error parametrization = const
#jet error parametrization   = jet et
#jet error parametrization   = const

start values = 1.0
#jet start values = -9.0 5.4 -0.5 -1.9 22.4 -105.5 0.32 -40.3 2.13 -4.86
jet start values = 7.4 -21.4 46.0 0.95 -117.8 -410.2 -0.3333 2.97 137.5 -0.015 .0.18 26.1 1.40 275.2

# Scaling of residuals: 0 - none, 1 - with Cauchy-Function, 2 - with Huber-Function
Residual Scaling Scheme    = 111 #   221 : default
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

jet granularity in eta = 3 #  allowed values are: 1,3,5,11,21,41 (*2)
jet granularity in phi = 1 #   1 : default

track granularity in eta = 1
track granularity in phi = 1
 
jet binning variables = eta;pt;sigmaphi
jet binning eta bins = -5.191 -2.964 -1.392  1.392 2.964 5.191
jet binning pt bins = 5 10 12 15 18 22 26 30 35 40 45 51 57 64 72 80 90 105 120 135 150 175 200 250 300 350 400 500 650 800 1000 1500 5000
jet binning emf bins = 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
jet binning sigmaphi bins = 0 0.02 0.04 0.06 0.08 0.09 0.10 0.11 0.12 0.13 0.14 0.15 0.16 0.18 0.20 0.22 0.24 1.0
#---------------------------------------------------------------------------------
#   Analyses (GammaJet, TrackTower, ... )
#---------------------------------------------------------------------------------
Et cut on jet              = 0.0
Et cut on gamma            = 0.0
Et cut on Z                = 20.0
Et cut on tower            = 0.0
Et cut on cluster          = 0.0
Et cut on track            = 0.0
Et cut on n+1 Jet          = 0.0
Eta cut on jet             = 6.0
Relative Rest Jet Cut      = 0.2      #NonLeadingJetsEt / PhotonEt
Min had fraction           = -0.05    #Default: 0.07
Max had fraction           = 1.05    #Default: 0.95
Et genJet min              = 0.0
Et genJet max              = 7000.0
DeltaR cut on jet matching = 0.25

#---------------------------------------------------------------------------------
#   Input / Output
#---------------------------------------------------------------------------------

#jet correction source = JetMETCor
#jet correction name   = Spring10_AK5TRK
Default Tree Name      = CalibTree

# List of input files:
Gamma-Jet tree         = GammaJetTree
Z-Jet tree             = ZJetTree
Track-Tower tree       = TrackTowerTree
Track-Cluster tree     = TrackClusterTree
Di-Jet tree            = DiJetTree
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

correct jets L2L3 = true
#Di-Jet prescale = 1000
#-----------------------------------------------------------------
#  Control plots
#-----------------------------------------------------------------
#  General parameters
create plots                     = true
plots output directory           = plots
#plots format                      = pdf

# JetTruthEvent plots
create JetTruthEvent plots    =  true

JetTruthEvent plots names =  MCTruthResponseVsGenJetPt; MCTruthResponseVsEta; MCTruthResponseVsPhiPhi; MCTruthResponseVsEtaEta;MCTruthRespFlavorVsGenJetPt
MCTruthResponseVsGenJetPt x variable        =  GenJetPt;  log
MCTruthResponseVsGenJetPt x edges           =  30 10 2000
MCTruthResponseVsGenJetPt y variable        =  GenJetResponse
MCTruthResponseVsGenJetPt y edges           =  51 0 2.0 0.9 1.1 0.9 1.1 0.0 0.5
MCTruthResponseVsGenJetPt bin variable      =  Eta
MCTruthResponseVsGenJetPt bin edges         =  -5.0 -3.0 -1.3 1.3 3.0 5.0
MCTruthResponseVsGenJetPt correction types  =  Uncorrected; Kalibri
MCTruthResponseVsGenJetPt profile types     =  Mean; GaussFitMean; GaussFitWidth
#MCTruthResponseVsGenJetPt distributions    =  Uncorrected; Kalibri
MCTruthResponseVsGenJetPt legend label      =  Uncorrected:CMS L2L3

MCTruthResponseVsEta x variable         =  Eta
MCTruthResponseVsEta x edges            =  20 -5 5
MCTruthResponseVsEta y variable         =  GenJetResponse
MCTruthResponseVsEta y edges            =  51 0 2.0 0.9 1.1 0.9 1.1 0.0 0.5
MCTruthResponseVsEta bin variable       =  GenJetPt
MCTruthResponseVsEta bin edges          =  10 50 100 500 2000
MCTruthResponseVsEta correction types   =  Uncorrected; Kalibri
MCTruthResponseVsEta profile types      =  Mean; GaussFitMean; GaussFitWidth
#MCTruthResponseVsEta distributions      =  Uncorrected; Kalibri
MCTruthResponseVsEta legend label       =  Uncorrected:CMS L2L3

MCTruthResponseVsPhiPhi x variable         =  momentPhiPhi
MCTruthResponseVsPhiPhi x edges            =  20 0 0.5
MCTruthResponseVsPhiPhi y variable         =  GenJetResponse
MCTruthResponseVsPhiPhi y edges            =  51 0 2 0.9 1.1 0.9 1.1 0.0 0.5
MCTruthResponseVsPhiPhi bin variable       =  GenJetPt
MCTruthResponseVsPhiPhi bin edges          =  10 30 50 80 120 300 600 2000
MCTruthResponseVsPhiPhi correction types   =  Uncorrected; Kalibri
MCTruthResponseVsPhiPhi profile types      =  Mean; GaussFitMean; GaussFitWidth
#MCTruthResponseVsPhiPhi distributions     =  Uncorrected; Kalibri
MCTruthResponseVsPhiPhi legend label       =  Uncorrected:CMS L2L3

MCTruthResponseVsEtaEta x variable         =  momentEtaEta
MCTruthResponseVsEtaEta x edges            =  20 0 0.5
MCTruthResponseVsEtaEta y variable         =  GenJetResponse
MCTruthResponseVsEtaEta y edges            =  1 0 2 0.9 1.1 0.9 1.1 0.0 0.5
MCTruthResponseVsEtaEta bin variable       =  GenJetPt
MCTruthResponseVsEtaEta bin edges          =  10 30 50 80 120 300 600 2000
MCTruthResponseVsEtaEta correction types   =  Uncorrected; Kalibri
MCTruthResponseVsEtaEta profile types      =  Mean; GaussFitMean; GaussFitWidth
#MCTruthResponseVsEtaEta distributions     =  Uncorrected; Kalibri
MCTruthResponseVsEtaEta legend label       =  Uncorrected:CMS L2L3

MCTruthRespFlavorVsGenJetPt x variable        =  GenJetPt;  log
MCTruthRespFlavorVsGenJetPt x edges           =  30 10 3000
MCTruthRespFlavorVsGenJetPt y variable        =  GenJetResponse
MCTruthRespFlavorVsGenJetPt y edges           =  51 0 2 0.9 1.1 0.9 1.1
MCTruthRespFlavorVsGenJetPt bin variable      =  Flavor
MCTruthRespFlavorVsGenJetPt bin edges         =  -0.5 0.5 1.5
MCTruthRespFlavorVsGenJetPt correction types  =  Uncorrected; Kalibri
MCTruthRespFlavorVsGenJetPt profile types     =  Mean; GaussFitMean
#MCTruthResponseVsGenJetPt distributions      =  Uncorrected; Kalibri
MCTruthRespFlavorVsGenJetPt legend label      =  Uncorrected:CMS L2L3

"""
    fcfg = open(filename, "w")
    fcfg.write(config)
    fcfg.write("Di-Jet input file = dijetlist\n")
    fcfg.write("Output file       = "+output+"\n");
    fcfg.write("Number of Threads = "+str(nthreads)+"\n")
    if(useconstraint):
        fcfg.write("jet constraints =  5.0 10.0 0.0 1.2 1 10.0 15.0 0.0 1.2 1 15.0 20.0 0 1.2 1 20.0 25.0 0 1.2 1 25.0 30.0 0 1.2 1 30.0 40.0 0 1.2 1 40.0 50.0 0 1.2  1 50.0 60.0 0 1.2 1 60.0 70.0 0 1.2  1 70.0 80 0 1.2 1  80.0 90.0 0 1.2 1 90.0 100.0 0 1.2 1 100.0 120. 0 1.2 1 120 150 0 1.2 1 150 200 0 1.2 1 200 280 0 1.2 1 280 350 0 1.2 1 350 500 0 1.2 1 500 800 0 1.2 1 800 1400 0 1.2 1 1400 7000 0 1.2 1\n")
        
    if(binned):
        fcfg.write("Di-Jet data class    = 21\n")
        fcfg.write("jet error parametrization   = const\n")
    else:
        fcfg.write("Di-Jet data class    = 11\n")
        fcfg.write("jet error parametrization   = jet et\n")
            
    fcfg.write("use Di-Jet events = "+str(nevents)+"\n")
            
    if(input != ""):
        fcfg.write("input calibration = Kalibri; "+input+"\n");


    fcfg.write("correct jets L1 = true");
    fcfg.close()
    return


#main program starts here!!!
#change these variables to steer the fit
jettype = "Calo"
datadir = "/scratch/hh/current/cms/user/stadie/QCDFlat_Pt15to3000Spring10-START3X_V26_S09-v1C"
nthreads = 5
nevents =  -1
dirname = "L4small"
useconstraint = False
batch = True
doBinnedFit = True
doUnbinnedFit = True


#write configs and run the fit
print "fit L4 correction for "+jettype+" using "+datadir

if os.path.exists(dirname):
    os.system("rm -rf "+dirname+"/*")
else:
    os.system("mkdir "+dirname)
    
os.system("cp junk "+dirname+"/.")
os.system("cd "+dirname+"; ln -s ../kalibriLogoSmall.gif ; cd -")

os.system("ls "+datadir+"/*"+jettype+"*.root > "+dirname+"/dijetlist");
binned= True
output="Kalibri.txt"
input =""
writeCfg(dirname+"/L4.cfg")
binned = False
output="Kalibri2.txt"
if doBinnedFit:
    input ="Kalibri.txt"

writeCfg(dirname+"/L4b.cfg")

if batch:
    fjob = open(dirname+"/fitL4.sh", "w")
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
        fjob.write("./junk L4.cfg > $TMPDIR/L4.log\n")
        fjob.write("date\n")
        fjob.write("mv $TMPDIR/L2L3.log .\n")

    if doUnbinnedFit:
        fjob.write("./junk L4b.cfg > $TMPDIR/L4b.log\n")
        fjob.write("date\n")
        fjob.write("mv $TMPDIR/L4b.log .\n")

    fjob.close()
    qsubcmd = "qsub "+dirname+"/fitL4.sh"
    print "running "+qsubcmd
    os.system(qsubcmd)
else:
    if doBinnedFit:
        kalibricmd = "cd "+dirname+"; ./junk L4.cfg; cd -";
        print "running "+kalibricmd
        os.system(kalibricmd)

    if doUnbinnedFit:
        kalibricmd = "cd "+dirname+"; ./junk L4b.cfg; cd -";
        print "running "+kalibricmd
        os.system(kalibricmd)


