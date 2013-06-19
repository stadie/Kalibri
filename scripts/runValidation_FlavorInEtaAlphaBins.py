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
Max cut on relative n+1 Jet Et = 1.0 ## to get more dijet-like events

#---------------------------------------------------------------------------------
#   Input / Output
#---------------------------------------------------------------------------------
Number of IO Threads = 5

#deprecated#jet correction source  = JetMETCor
#deprecated#jet correction name    = 2013SummerV1_AK5FastPFCHS_MC
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



"""

jettypes = ["ak5FastPFCHS","ak5FastPF"]#"ak5PFCHS","ak5withNuPFCHS","ak5FastPF", "ak5FastCalo"]

datadir = "/afs/naf.desy.de/user/d/draeger/public/PYTHIAZ2Star/" #full?!
datadirmc = "/afs/naf.desy.de/user/d/draeger/public/PYTHIAZ2/" #full?!
datadir2mc = "/afs/naf.desy.de/user/d/draeger/public/Herwig/" #fast?!
jecname = "2013SummerV1"
datasetname="Z2star_Z2_Herwig_Comparison_various_alpha_cuts"



correctJets=True
CutAwayPU=False

flavorDefs = ["phys", "algo"] #["algo", "phys"]

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


        labelforcorrections = "CMS L1L2L3"
        #change labels
        if CutAwayPU:
            labelforcorrections = "CMS L2L3"
    
        fcfg.write(config)


        eta_min_list = ['0.0','1.3','2.5','3.2']
        eta_max_list = ['1.3','2.5','3.2','5.2']
        eta_suffix_list = ['0013','1325','2532','3252']

        cut_list=['40','30','20','10']
        cut_no_list=['.40','.30','.20','.10']

        fcfg.write("JetTruthEvent plots cut_list =   ")
        for index_cut, cut in enumerate(cut_list):
            fcfg.write(cut + " ")
        fcfg.write("\n")
        
        fcfg.write("JetTruthEvent plots cut_no_list =   ")
        for index_cut, cut in enumerate(cut_no_list):
            fcfg.write(cut + " ")
            
        fcfg.write("\n")
            
        
        
        fcfg.write("JetTruthEvent plots names =   ")
        for index_cut, cut in enumerate(cut_list):
            fcfg.write("MCTruthResponseVsGenJetPt_" +cut + ";")
            for index_eta, eta_suffix in enumerate(eta_suffix_list):
                fcfg.write("MCTruthRespFlavorVsGenJetPt"+eta_suffix+ "_" +cut + ";")
        fcfg.write("\n")



        for index_cut, cut in enumerate(cut_list):
            fcfg.write("MCTruthResponseVsGenJetPt_" +cut + " x variable        =  GenJetPt;  log\n")
            fcfg.write("MCTruthResponseVsGenJetPt_" +cut + " x edges           =  50 10 3000\n")
            fcfg.write("MCTruthResponseVsGenJetPt_" +cut + " y variable        =  GenJetResponse\n")
            fcfg.write("MCTruthResponseVsGenJetPt_" +cut + " y edges           =  51 0 2.0 0.9 1.1 0.9 1.1 0.0 0.5\n")
            fcfg.write("MCTruthResponseVsGenJetPt_" +cut + " bin variable      =  AbsEta\n")
            fcfg.write("MCTruthResponseVsGenJetPt_" +cut + " bin edges         =  0 1.3 2.5 3.0 3.2 5.2\n")
            fcfg.write("MCTruthResponseVsGenJetPt_" +cut + " cut2 variable     =  ThirdJetFractionPlain\n")
            fcfg.write("MCTruthResponseVsGenJetPt_" +cut + " cut2 edges        =  0.0 " + cut_no_list[index_cut]+ "\n")
            fcfg.write("MCTruthResponseVsGenJetPt_" +cut + " correction types  =  L2L3 ### Uncorrected; L2L3\n")
            fcfg.write("MCTruthResponseVsGenJetPt_" +cut + " profile types     =  Mean; GaussFitMean; GaussFitWidth\n")
            fcfg.write("MCTruthResponseVsGenJetPt_" +cut + " legend label      =  L2L3:" + labelforcorrections + "\n")
            fcfg.write("MCTruthResponseVsGenJetPt_" +cut + " input samples     =  0:Z2star; 1:Z2; 2:Herwig\n")
            fcfg.write("\n")
            fcfg.write("\n")

            for index_eta, eta_suffix in enumerate(eta_suffix_list):
                fcfg.write("MCTruthRespFlavorVsGenJetPt"+eta_suffix+ "_" +cut + " x variable        =  GenJetPt;  log\n")
                fcfg.write("MCTruthRespFlavorVsGenJetPt"+eta_suffix+ "_" +cut + " x edges           =  50 10 3000\n")
                fcfg.write("MCTruthRespFlavorVsGenJetPt"+eta_suffix+ "_" +cut + " y variable        =  GenJetResponse\n")
                fcfg.write("MCTruthRespFlavorVsGenJetPt"+eta_suffix+ "_" +cut + " y edges           =  51 0 2 0.9 1.1 0.9 1.1\n")
                fcfg.write("MCTruthRespFlavorVsGenJetPt"+eta_suffix+ "_" +cut + " bin variable      =  Flavor\n")
                fcfg.write("MCTruthRespFlavorVsGenJetPt"+eta_suffix+ "_" +cut + " bin edges         =  -1.5 -0.5 0.5 1.5 2.5 3.5\n")
                fcfg.write("MCTruthRespFlavorVsGenJetPt"+eta_suffix+ "_" +cut + " cut variable      =  AbsEta\n")
                fcfg.write("MCTruthRespFlavorVsGenJetPt"+eta_suffix+ "_" +cut + " cut edges         =  "+eta_min_list[index_eta] + " " + eta_max_list[index_eta] +"\n")
                fcfg.write("MCTruthRespFlavorVsGenJetPt"+eta_suffix+ "_" +cut + " cut2 variable     =  ThirdJetFractionPlain\n")
                fcfg.write("MCTruthRespFlavorVsGenJetPt"+eta_suffix+ "_" +cut + " cut2 edges        =  0.0 " + cut_no_list[index_cut]+ "\n")
                fcfg.write("MCTruthRespFlavorVsGenJetPt"+eta_suffix+ "_" +cut + " correction types  =  L2L3 ### Uncorrected; L2L3\n")
                fcfg.write("MCTruthRespFlavorVsGenJetPt"+eta_suffix+ "_" +cut + " profile types     =  Mean; GaussFitMean\n")
                fcfg.write("MCTruthRespFlavorVsGenJetPt"+eta_suffix+ "_" +cut + " legend label      =  L2L3:" + labelforcorrections + "\n")
                fcfg.write("MCTruthRespFlavorVsGenJetPt"+eta_suffix+ "_" +cut + " input samples     =  0:Z2star; 1:Z2; 2:Herwig\n")


        
        if(flavorDef=="phys\n"):
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
    

