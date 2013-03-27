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
Number of Threads = 5

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
use Di-Jet events          = -1
use Di-Jet Control1 events = -1
use Di-Jet Control2 events = -1
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

PUTruthReweighting = false #true
Di-Jet trigger names = HLT_PFJet40
Di-Jet trigger thresholds = 68
PU TruthWeighting = kirschen/PUDistributions/TTJetDists
PU TruthWeighting MC distribution = kirschen/PUDistributions/TTJetDists/FullSimTTJets



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

JetTruthEvent plots names =  MCTruthResponseVsGenJetPt; MCTruthResponseVsEta; MCTruthResponsePU; MCTruthResolPU; MCTruthUnmatchedAndGluonsRespFlavorVsGenJetPt;MCTruthRespFlavorVsGenJetPt;MCTruthEndCapRespFlavorVsGenJetPt;MCTruthHFRespFlavorVsGenJetPt; MCTruthResponseVsMeanWidth; MCTruthResponseVsClosestJetdR
MCTruthResponseVsGenJetPt x variable        =  GenJetPt;  log
MCTruthResponseVsGenJetPt x edges           =  50 10 3000
MCTruthResponseVsGenJetPt y variable        =  GenJetResponse
MCTruthResponseVsGenJetPt y edges           =  51 0 2.0 0.9 1.1 0.9 1.1 0.0 0.5
MCTruthResponseVsGenJetPt bin variable      =  Eta
MCTruthResponseVsGenJetPt bin edges         =  -5.0 -3.0 -1.3 1.3 3.0 5.0
MCTruthResponseVsGenJetPt correction types  =  L2L3 ### Uncorrected; L2L3
#; L2L3L4
MCTruthResponseVsGenJetPt profile types     =  Mean; GaussFitMean; GaussFitWidth
MCTruthResponseVsGenJetPt legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResponseVsGenJetPt input samples     =  0:Z2star; 1:Z2; 2:Herwig

MCTruthResponseVsEta x variable         =  Eta
MCTruthResponseVsEta x edges            =  40 -5 5
MCTruthResponseVsEta y variable         =  GenJetResponse
MCTruthResponseVsEta y edges            =  51 0 2.0 0.9 1.1 0.9 1.1
MCTruthResponseVsEta bin variable       =  GenJetPt
MCTruthResponseVsEta bin edges          =  20 50 100 500 2000
MCTruthResponseVsEta correction types   =  L2L3 ### Uncorrected; L2L3
#; L2L3L4 
MCTruthResponseVsEta profile types      =  Mean; GaussFitMean
MCTruthResponseVsEta legend label       =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResponseVsEta input samples     =  0:Z2star; 1:Z2; 2:Herwig

MCTruthResponsePU x variable         =  Eta
MCTruthResponsePU x edges            =  40 -5 5
MCTruthResponsePU y variable         =  GenJetResponse
MCTruthResponsePU y edges            =  51 0 2.0 0.9 1.1 0.9 1.1
MCTruthResponsePU bin variable       =  NPU
MCTruthResponsePU bin edges          =  0 1 10 20 30
MCTruthResponsePU cut variable       =  GenJetPt
MCTruthResponsePU cut edges          =  50 100
MCTruthResponsePU correction types   =  L2L3 ### Uncorrected; L2L3
#; L2L3L4 
MCTruthResponsePU profile types      =  Mean; GaussFitMean
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
MCTruthResolPU correction types  =  L2L3 ### Uncorrected; L2L3
#; L2L3L4
MCTruthResolPU profile types     =  GaussFitWidth
MCTruthResolPU legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolPU input samples     =  0:Z2star; 1:Z2; 2:Herwig

MCTruthResolVsGenJetPt x variable        =  GenJetPt;  log
MCTruthResolVsGenJetPt x edges           =  50 10 3000
MCTruthResolVsGenJetPt y variable        =  GenJetResponse
MCTruthResolVsGenJetPt y edges           =  51 0 1.0 0 0.5
#MCTruthResolVsGenJetPt bin variable      =  Eta
#MCTruthResolVsGenJetPt bin edges         =  -5.0 -3.0 -1.3 1.3 3.0 5.0
MCTruthResolVsGenJetPt bin variable      =  AbsEta
MCTruthResolVsGenJetPt bin edges         =  0 1.3 3.0 5.0
MCTruthResolVsGenJetPt correction types  =  L2L3 ### Uncorrected; L2L3
#; L2L3L4
MCTruthResolVsGenJetPt profile types     =  GaussFitWidth
MCTruthResolVsGenJetPt legend label      =  L2L3:CMS L1L2L3
#; L2L3L4:CMS L2L3 + L4JW
MCTruthResolVsGenJetPt input samples     =  0:Z2star; 1:Z2; 2:Herwig

MCTruthUnmatchedAndGluonsRespFlavorVsGenJetPt x variable        =  GenJetPt;  log
MCTruthUnmatchedAndGluonsRespFlavorVsGenJetPt x edges           =  50 10 3000
MCTruthUnmatchedAndGluonsRespFlavorVsGenJetPt y variable        =  GenJetResponse
MCTruthUnmatchedAndGluonsRespFlavorVsGenJetPt y edges           =  51 0 2 0.9 1.1 0.9 1.1
MCTruthUnmatchedAndGluonsRespFlavorVsGenJetPt bin variable      =  Flavor
MCTruthUnmatchedAndGluonsRespFlavorVsGenJetPt bin edges         =  -1.5 0.5
MCTruthUnmatchedAndGluonsRespFlavorVsGenJetPt cut variable      =  AbsEta
MCTruthUnmatchedAndGluonsRespFlavorVsGenJetPt cut edges         =  0 1.3
MCTruthUnmatchedAndGluonsRespFlavorVsGenJetPt correction types  =  L2L3 ### Uncorrected; L2L3
MCTruthUnmatchedAndGluonsRespFlavorVsGenJetPt profile types     =  Mean; GaussFitMean
MCTruthUnmatchedAndGluonsRespFlavorVsGenJetPt legend label      =  L2L3:CMS L1L2L3
MCTruthUnmatchedAndGluonsRespFlavorVsGenJetPt input samples     =  0:Z2star; 1:Z2; 2:Herwig


MCTruthRespFlavorVsGenJetPt x variable        =  GenJetPt;  log
MCTruthRespFlavorVsGenJetPt x edges           =  50 10 3000
MCTruthRespFlavorVsGenJetPt y variable        =  GenJetResponse
MCTruthRespFlavorVsGenJetPt y edges           =  51 0 2 0.9 1.1 0.9 1.1
MCTruthRespFlavorVsGenJetPt bin variable      =  Flavor
MCTruthRespFlavorVsGenJetPt bin edges         =  -1.5 -0.5 0.5 1.5 2.5 3.5
MCTruthRespFlavorVsGenJetPt cut variable      =  AbsEta
MCTruthRespFlavorVsGenJetPt cut edges         =  0 1.3
MCTruthRespFlavorVsGenJetPt correction types  =  L2L3 ### Uncorrected; L2L3
MCTruthRespFlavorVsGenJetPt profile types     =  Mean; GaussFitMean
MCTruthRespFlavorVsGenJetPt legend label      =  L2L3:CMS L1L2L3
MCTruthRespFlavorVsGenJetPt input samples     =  0:Z2star; 1:Z2; 2:Herwig

MCTruthEndCapRespFlavorVsGenJetPt x variable        =  GenJetPt;  log
MCTruthEndCapRespFlavorVsGenJetPt x edges           =  50 10 3000
MCTruthEndCapRespFlavorVsGenJetPt y variable        =  GenJetResponse
MCTruthEndCapRespFlavorVsGenJetPt y edges           =  51 0 2 0.9 1.1 0.9 1.1
MCTruthEndCapRespFlavorVsGenJetPt bin variable      =  Flavor
MCTruthEndCapRespFlavorVsGenJetPt bin edges         =  -1.5 -0.5 0.5 1.5 2.5 3.5
MCTruthEndCapRespFlavorVsGenJetPt cut variable      =  AbsEta
MCTruthEndCapRespFlavorVsGenJetPt cut edges         =  1.3 3.0 
MCTruthEndCapRespFlavorVsGenJetPt correction types  =  L2L3 ### Uncorrected; L2L3
MCTruthEndCapRespFlavorVsGenJetPt profile types     =  Mean; GaussFitMean
MCTruthEndCapRespFlavorVsGenJetPt legend label      =  L2L3:CMS L1L2L3
MCTruthEndCapRespFlavorVsGenJetPt input samples     =  0:Z2star; 1:Z2; 2:Herwig


MCTruthHFRespFlavorVsGenJetPt x variable        =  GenJetPt;  log
MCTruthHFRespFlavorVsGenJetPt x edges           =  50 10 3000
MCTruthHFRespFlavorVsGenJetPt y variable        =  GenJetResponse
MCTruthHFRespFlavorVsGenJetPt y edges           =  51 0 2 0.9 1.1 0.9 1.1
MCTruthHFRespFlavorVsGenJetPt bin variable      =  Flavor
MCTruthHFRespFlavorVsGenJetPt bin edges         =  -1.5 -0.5 0.5 1.5 2.5 3.5
MCTruthHFRespFlavorVsGenJetPt cut variable      =  AbsEta
MCTruthHFRespFlavorVsGenJetPt cut edges         =  3.0 5.0 
MCTruthHFRespFlavorVsGenJetPt correction types  =  L2L3 ### Uncorrected; L2L3
MCTruthHFRespFlavorVsGenJetPt profile types     =  Mean; GaussFitMean
MCTruthHFRespFlavorVsGenJetPt legend label      =  L2L3:CMS L1L2L3
MCTruthHFRespFlavorVsGenJetPt input samples     =  0:Z2star; 1:Z2; 2:Herwig

MCTruthResponseVsMeanWidth x variable         =  meanMoment
MCTruthResponseVsMeanWidth x edges            =  30 0 0.5
MCTruthResponseVsMeanWidth y variable         =  GenJetResponse
MCTruthResponseVsMeanWidth y edges            =  51 0 2 0.9 1.1 0.9 1.1 0.0 0.5
MCTruthResponseVsMeanWidth bin variable       =  GenJetPt
MCTruthResponseVsMeanWidth bin edges          =  10 30 50 80 120 300 600 2000
MCTruthResponseVsMeanWidth correction types   =  L2L3 ### Uncorrected; L2L3
MCTruthResponseVsMeanWidth profile types      =  Mean; GaussFitMean; GaussFitWidth
MCTruthResponseVsMeanWidth legend label       =  L2L3:CMS L1L2L3
MCTruthResponseVsMeanWidth input samples     =  0:Z2star; 1:Z2; 2:Herwig


MCTruthResponseVsClosestJetdR x variable         =  ClosestJetdR
MCTruthResponseVsClosestJetdR x edges            =  20 0 2
MCTruthResponseVsClosestJetdR y variable         =  GenJetResponse
MCTruthResponseVsClosestJetdR y edges            =  51 0 2.0 0.9 1.1 0.9 1.1
MCTruthResponseVsClosestJetdR bin variable       =  GenJetPt
MCTruthResponseVsClosestJetdR bin edges          =  20 50 100 500 2000
MCTruthResponseVsClosestJetdR cut variable       =  AbsEta
MCTruthResponseVsClosestJetdR cut edges          =  0.0 1.3 
MCTruthResponseVsClosestJetdR correction types   =  L2L3 ### Uncorrected; L2L3
MCTruthResponseVsClosestJetdR profile types      =  Mean; GaussFitMean
MCTruthResponseVsClosestJetdR legend label       =  L2L3:CMS L1L2L3
MCTruthResponseVsClosestJetdR input samples      =  0:Z2star; 1:Z2; 2:Herwig




"""

jettypes = ["ak5PFCHS" ,"ak5withNuPFCHS","ak5FastPF", "ak5FastCalo"]

datadir = "/scratch/hh/dust/naf/cms/user/kriheine/CalibNTupel/MC/Z2star_pythia_v3/" #full?!
datadirmc = "/scratch/hh/dust/naf/cms/user/kriheine/CalibNTupel/MC/Z2_pythia_v3/" #full?!
datadir2mc = "/scratch/hh/dust/naf/cms/user/kriheine/CalibNTupel/MC/EE3C_herwigpp_v3/" #fast?!
jecname = "2012FallV5"
datasetname="Z2star_Z2_Herwig_Comparison"



correctJets=False
CutAwayPU=False

flavorDefs = ["algo", "phys"]

for flavorDef in flavorDefs:
    for jettype in jettypes:
        print "make plots for jettype "+jettype
        if os.path.exists("tempdijetlist"):
            os.remove("tempdijetlist")
            os.remove("mcdijetlist")
    ##        os.remove("mc2dijetlist")
        if os.path.exists("tempplots"):
            os.system("rm tempplots/*")
    
    #    os.system("ls "+datadir+"/*_"+jettype+"*.root > tempdijetlist");
    #    os.system("ls "+datadirmc+"/*_"+jettype+"*.root > mcdijetlist");
    #    os.system("ls "+datadir2mc+"/*_"+jettype+"*.root > mc2dijetlist");
        os.system("ls "+datadir+"*_"+jettype+".root > tempdijetlist");
        os.system("ls "+datadirmc+"*_"+jettype+".root > mcdijetlist");
        os.system("ls "+datadir2mc+"*_"+jettype+".root > mc2dijetlist");
        fcfg = open("valid.cfg", "w")
        #change labels
        if CutAwayPU:
            config = config.replace('CMS L1L2L3','CMS L2L3')
    
        fcfg.write(config)
        if(flavorDef=="phys"):
            fcfg.write("Use physical flavor definition instead of algorithmic = true\n")
        else:
            fcfg.write("Use physical flavor definition instead of algorithmic = false\n")
        fcfg.write("plots output directory = tempplots\n")
        fcfg.write("Di-Jet input file = tempdijetlist\n")
        fcfg.write("Di-Jet Control1 input file = mcdijetlist\n")
        fcfg.write("Di-Jet Control2 input file = mc2dijetlist\n")
        jetalgo = jettype[0:3]
        
        if correctJets:
            fcfg.write("jet correction source = JetMETCor\n");
            fcfg.write("jet correction name   = "+jecname+"_"+jetalgo.upper()+jettype[3:len(jettype)]+"_MC \n");
            #fcfg.write("jet correction name   = "+jecname+"_"+jetalgo.upper()+jettype[3:len(jettype)]+"NoOffset\n");
        if CutAwayPU:
            #fcfg.write("MAX n PU from MC = 8\n");
            fcfg.write("correct jets L1 = false\n")
            jec = jecname+'NoPU'
        else:
            fcfg.write("correct jets L1 = true\n")
            jec = jecname
        fcfg.close()
        
        jec=jec+"_"+flavorDef
        kalibricmd = "./junk valid.cfg"
        print "running "+kalibricmd
        os.system(kalibricmd)
        tarball=jec+"plots"+jettype+".tar"
        tarcmd = "cd tempplots; tar cf ../"+tarball+" *Eta[0-9].eps *Pt[0-9].eps *Eta[0-9]_zoom.eps *Pt[0-9]_zoom.eps *Flavor[0-9].eps *Flavor[0-9]_zoom.eps *NPU[0-9].eps *NPU[0-9]_zoom.eps *.root; cd -"
        print "running "+tarcmd
        os.system(tarcmd)

print "please run:"
for flavorDef in flavorDefs:
    for jettype in jettypes:
        if CutAwayPU:
            jec = jecname+'NoPU'
        else:
            jec = jecname
        jec=jec+"_"+flavorDef
        tarball=os.getcwd()+"/"+jec+"plots"+jettype+".tar"
        
        webcmd = "./scripts/createJECValidationHtmlPage.sh jetmet "+tarball+" \""+jec+"\" \""+jecname+"\" "+datasetname+" "+jettype[0:3] +" "+jettype[3:len(jettype)].lower()
        print webcmd
        #os.system(webcmd)
    

