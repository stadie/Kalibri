#! /usr/bin/python

import os

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
Parametrization Class = L2L3JetParametrization

Number of IO Threads = 5

#Error Parametrization
tower error parametrization = const
jet error parametrization   = jet et
#jet error parametrization   = toy

start values = 1.0
jet start values        = 1. 0. 0. 0. 0. 0.
global jet start values = 0.998 4.997 3.084 2.048

#---------------------------------------------------------------------------------
#   Geometry / Binning for fitting
#---------------------------------------------------------------------------------
maximum eta twr used  = 82
granularity in eta    = 2   # allowed values are: 1,3,5,11,21,41 (*2)
granularity in phi    = 1   # allowed values are: 1,2,3,6,9,18,36,72
symmetry in eta       = false

jet granularity in eta = 82 #  allowed values are: 1,3,5,11,21,41 (*2)
jet granularity in phi = 1 #   1 : default

track granularity in eta = 1
track granularity in phi = 1

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
Eta max cut on jet         = 6.0
Relative Rest Jet Cut      = 0.2      #NonLeadingJetsEt / PhotonEt
Min had fraction           = -0.05    #Default: 0.07
Max had fraction           = 9999    #Default: 0.95
Et genJet min              = 0.0
Et genJet max              = 7000.0
DeltaR cut on jet matching = 0.25

#---------------------------------------------------------------------------------
#   Input / Output
#---------------------------------------------------------------------------------
Number of IO Threads = 5

#jet correction source = JetMETCor
#jet correction name   = Spring10_AK5TRK
Default Tree Name      = CalibTree

# List of input files:
Gamma-Jet tree         = GammaJetTree
Z-Jet tree             = ZJetTree
Track-Tower tree       = TrackTowerTree
Track-Cluster tree     = TrackClusterTree
Di-Jet tree            = DiJetTree
Di-Jet Control1 tree   = DiJetTree
Di-Jet Control2 tree   = DiJetTree
Tri-Jet tree           = TriJetTree
Top tree               = TopTree

#Di-Jet input file = input/dijetlistspring10Track
#Di-Jet input file =  input/toy_dijet_const_5-500_uniform.root

#Top input file = /scratch/current/cms/user/stadie/Top_Madgraph.root
#Z-Jet input file = /scratch/current/cms/user/stadie/ZJet_Track_230_300_rereco.root; /scratch/current/cms/user/stadie/ZJet_Track_300_INF_rereco.root; 

# -1: use all events, 1000: run over 1000 events
use Gamma-Jet events       = 0
use Track-Tower events     = 0
use Track-Cluster events   = 0
use Di-Jet events          = 100000
use Di-Jet Control1 events = 100000
use Di-Jet Control2 events = 100000
use Tri-Jet events         = 0
use Z-Jet events           = 0
use Top events             = 0

Gamma-Jet data class     = 1
Z-Jet data class     = 3
Di-Jet data class    = 11
Top data class       = 1


#---------------------------------------------------------------------------------
#   Event processor setup
#---------------------------------------------------------------------------------

PUTruthReweighting = false

PU TruthWeighting Reweight all eventvectors (for MC validation) = true
PU TruthWeighting = kheine/PUDistributions/Cert_2012_190456-208686_ReReco
PU TruthWeighting MC distribution = kirschen/PUDistributions/TrueDistributions/Summer12S10CMSSW53
Di-Jet trigger names = HLT_PFJet400 ## will reweight to PUTruth distribution of unprescaled PF400-trigger as determined by Kristin
Di-Jet trigger thresholds = 5 ## set threshold low so all JetTruthEvents should pass



#Di-Jet prescale = 1000
#-----------------------------------------------------------------
#  Control plots
#-----------------------------------------------------------------
#  General parameters
create plots                     = true
#plots output directory           = L2L3Plots
#plots format                      = pdf

# JetTruthEvent plots
create JetTruthEvent plots    =  true

JetTruthEvent plots names =  MCTruthResponseVsNPUTruth; MCTruthResponseVsGenJetPt; MCTruthResponseVsEta; MCTruthResponsePU; MCTruthResolPU#; MCTruthRespFlavorVsGenJetPt
MCTruthResponseVsNPUTruth x variable        =  NPUTruth
MCTruthResponseVsNPUTruth x edges           =  50 0 50
MCTruthResponseVsNPUTruth y variable        =  GenJetResponse
MCTruthResponseVsNPUTruth y edges           =  51 0 2.0 0.9 1.1 0.0 0.5
MCTruthResponseVsNPUTruth bin variable      =  Eta
MCTruthResponseVsNPUTruth bin edges         =  -5.0 5.0
MCTruthResponseVsNPUTruth correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResponseVsNPUTruth profile types     =  GaussFitMean; GaussFitWidth
MCTruthResponseVsNPUTruth legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResponseVsNPUTruth input samples     =  0:Z2star; 1:Z2; 2:Herwig




MCTruthResponseVsGenJetPt x variable        =  GenJetPt;  log
MCTruthResponseVsGenJetPt x edges           =  50 10 3000
MCTruthResponseVsGenJetPt y variable        =  GenJetResponse
MCTruthResponseVsGenJetPt y edges           =  51 0 2.0 0.9 1.1 0.0 0.5
MCTruthResponseVsGenJetPt bin variable      =  Eta
MCTruthResponseVsGenJetPt bin edges         =  -5.0 -3.0 -1.3 1.3 3.0 5.0
MCTruthResponseVsGenJetPt correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResponseVsGenJetPt profile types     =  GaussFitMean; GaussFitWidth
MCTruthResponseVsGenJetPt legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResponseVsGenJetPt input samples     =  0:Z2star; 1:Z2; 2:Herwig

MCTruthResponseVsEta x variable         =  Eta
MCTruthResponseVsEta x edges            =  40 -5 5
MCTruthResponseVsEta y variable         =  GenJetResponse
MCTruthResponseVsEta y edges            =  51 0 2.0 0.9 1.1
MCTruthResponseVsEta bin variable       =  GenJetPt
MCTruthResponseVsEta bin edges          =  20 50 100 500 2000
MCTruthResponseVsEta correction types   =  Uncorrected; L2L3
#; L2L3L4 
MCTruthResponseVsEta profile types      =  GaussFitMean
MCTruthResponseVsEta legend label       =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResponseVsEta input samples     =  0:Z2star; 1:Z2; 2:Herwig

MCTruthResponsePU x variable         =  Eta
MCTruthResponsePU x edges            =  40 -5 5
MCTruthResponsePU y variable         =  GenJetResponse
MCTruthResponsePU y edges            =  51 0 2.0 0.9 1.1
MCTruthResponsePU bin variable       =  NPU
MCTruthResponsePU bin edges          =  0 1 10 20 30
MCTruthResponsePU cut variable       =  GenJetPt
MCTruthResponsePU cut edges          =  50 100
MCTruthResponsePU correction types   =  Uncorrected; L2L3
#; L2L3L4 
MCTruthResponsePU profile types      =  GaussFitMean
MCTruthResponsePU legend label       =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResponsePU input samples     =  0:Z2star; 1:Z2; 2:Herwig

MCTruthResolPU x variable        =  GenJetPt;  log
MCTruthResolPU x edges           =  25 10 3000
MCTruthResolPU y variable        =  GenJetResponse
MCTruthResolPU y edges           =  51 0 1.0 0 0.5
MCTruthResolPU bin variable      =  NPU
MCTruthResolPU bin edges         =  0 1 10 20 30
MCTruthResolPU cut variable   =  Eta
MCTruthResolPU cut edges      =  -1.3 1.3
MCTruthResolPU correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolPU profile types     =  GaussFitWidth
MCTruthResolPU legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPU input samples     =  0:Z2star; 1:Z2; 2:Herwig

MCTruthResolVsGenJetPt x variable        =  GenJetPt;  log
MCTruthResolVsGenJetPt x edges           =  50 10 3000
MCTruthResolVsGenJetPt y variable        =  GenJetResponse
MCTruthResolVsGenJetPt y edges           =  51 0 1.0 0 0.5
MCTruthResolVsGenJetPt bin variable      =  Eta
MCTruthResolVsGenJetPt bin edges         =  -5.0 -3.0 -1.3 1.3 3.0 5.0
MCTruthResolVsGenJetPt correction types  =  Uncorrected; L2L3
#; L2L3L4
MCTruthResolVsGenJetPt profile types     =  GaussFitWidth
MCTruthResolVsGenJetPt legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolVsGenJetPt input samples     =  0:Z2star; 1:Z2; 2:Herwig

MCTruthRespFlavorVsGenJetPt x variable        =  GenJetPt;  log
MCTruthRespFlavorVsGenJetPt x edges           =  50 10 3000
MCTruthRespFlavorVsGenJetPt y variable        =  GenJetResponse
MCTruthRespFlavorVsGenJetPt y edges           =  51 0 2 0.9 1.1 0.9 1.1
MCTruthRespFlavorVsGenJetPt bin variable      =  Flavor
MCTruthRespFlavorVsGenJetPt bin edges         =  -0.5 0.5 1.5 2.5 3.5
MCTruthRespFlavorVsGenJetPt correction types  =  Uncorrected; L2L3
MCTruthRespFlavorVsGenJetPt profile types     =  Mean; GaussFitMean
MCTruthRespFlavorVsGenJetPt legend label      =  L2L3:CMS L1L2L3
MCTruthRespFlavorVsGenJetPt input samples     =  0:Z2star; 1:Z2; 2:Herwig

MCTruthResponseVsMeanWidth x variable         =  meanMoment
MCTruthResponseVsMeanWidth x edges            =  30 0 0.5
MCTruthResponseVsMeanWidth y variable         =  GenJetResponse
MCTruthResponseVsMeanWidth y edges            =  51 0 2 0.9 1.1 0.9 1.1 0.0 0.5
MCTruthResponseVsMeanWidth bin variable       =  GenJetPt
MCTruthResponseVsMeanWidth bin edges          =  10 30 50 80 120 300 600 2000
MCTruthResponseVsMeanWidth correction types   =  Uncorrected; L2L3
MCTruthResponseVsMeanWidth profile types      =  Mean; GaussFitMean; GaussFitWidth
MCTruthResponseVsMeanWidth legend label       =  L2L3:CMS L1L2L3
MCTruthResponseVsMeanWidth input samples     =  0:Z2star; 1:Z2; 2:Herwig

"""

jettypes = ["ak5FastPF"]#, "ak5PFCHS"]
datadir = "/scratch/hh/dust/naf/cms/user/kriheine/CalibNTupel/MC/Z2star_pythia_v3" #Z2star
datadirmc = "/scratch/hh/dust/naf/cms/user/kriheine/CalibNTupel/MC/Z2_pythia_v3" #Z2
datadir2mc = "/scratch/hh/dust/naf/cms/user/kriheine/CalibNTupel/MC/EE3C_herwigpp_v3" #Herwig
jecname = "Dec12_Z2star_Z2_Hpp"
datasetname="Z2star_Herwig_Comparison_TEST"

correctJets=False
CutAwayPU=False


for jettype in jettypes:
    print "make plots for jettype "+jettype
    if os.path.exists("tempdijetlist"):
        os.remove("tempdijetlist")
        os.remove("mcdijetlist")
##        os.remove("mc2dijetlist")
    if os.path.exists("tempplots"):
        os.system("rm tempplots/*")

    os.system("ls "+datadir+"/*_"+jettype+"*.root > tempdijetlist");
    os.system("ls "+datadirmc+"/*_"+jettype+"*.root > mcdijetlist");
    os.system("ls "+datadir2mc+"/*_"+jettype+"*.root > mc2dijetlist");
#    os.system("ls "+datadir+"*_"+jettype+"*.root > tempdijetlist");
#    os.system("ls "+datadirmc+"*_"+jettype+"*.root > mcdijetlist");
#    os.system("ls "+datadir2mc+"*_"+jettype+"*.root > mc2dijetlist");
    fcfg = open("valid.cfg", "w")
    #change labels
    if CutAwayPU:
        config = config.replace('CMS L1L2L3','CMS L2L3')

    fcfg.write(config)
    fcfg.write("plots output directory = tempplots\n")
    fcfg.write("Di-Jet input file = tempdijetlist\n")
    fcfg.write("Di-Jet Control1 input file = mcdijetlist\n")
    fcfg.write("Di-Jet Control2 input file = mc2dijetlist\n")
    jetalgo = jettype[0:3]
    
    if correctJets:
        fcfg.write("jet correction source = JetMETCor\n");
        fcfg.write("jet correction name   = "+jecname+"_"+jetalgo.upper()+jettype[3:len(jettype)]+"\n");
        #fcfg.write("jet correction name   = "+jecname+"_"+jetalgo.upper()+jettype[3:len(jettype)]+"NoOffset\n");
    if CutAwayPU:
        #fcfg.write("MAX n PU from MC = 8\n");
        fcfg.write("correct jets L1 = false\n")
        jec = jecname+'NoPU'
    else:
        fcfg.write("correct jets L1 = true\n")
        jec = jecname
    fcfg.close()
    
    kalibricmd = "./junk valid.cfg"
    print "running "+kalibricmd
    os.system(kalibricmd)
    tarball=jec+"plots"+jettype+".tar"
    tarcmd = "cd tempplots; tar cf ../"+tarball+" *Eta[0-9].eps *Pt[0-9].eps *Eta[0-9]_zoom.eps *Pt[0-9]_zoom.eps *Flavor[0-9].eps *Flavor[0-9]_zoom.eps *NPU[0-9].eps *NPU[0-9]_zoom.eps *.root; cd -"
    print "running "+tarcmd
    os.system(tarcmd)

print "please run:"
for jettype in jettypes:
    tarball=os.getcwd()+"/"+jec+"plots"+jettype+".tar"
    
    webcmd = "./scripts/createJECValidationHtmlPage.sh jetmet "+tarball+" \""+jec+"\" \""+jecname+"\" "+datasetname+" "+jettype[0:3] +" "+jettype[3:len(jettype)].lower()
    print webcmd
    #os.system(webcmd)
    

