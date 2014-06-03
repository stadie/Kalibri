#! /usr/bin/python

import os
import sys
import ConfigParser
from runDiJets_CommonConfigs import BinningValues, TriggerNamesThresholds, PUWeightingInfo, determineDataDir, determineDataDirMC, importDatatypesNewTrigger, configureJERsmearing

Usage="""
##############################################################################################
##############################################################################################
####  This is the steering script to create appropriate Kalibri-style config files for    ####
####  the dijet L2Residual analysis and run Kalibri (the junk executable) on these        ####
####  config-files.                                                                       ####
####                                                                                      ####
####  In order to allow easy batch execution of multiple Kalibri instances defined in     ####
####  this way, there are a number of command line parameters, that can be used           ####
####                                                                                      ####
####  Usage:                                                                              ####
####    ./scripts/runDiJets_Residuals.py ro run the script                                ####
####                                                                                      ####
####  Cmd-line options (8 arguments need to be supplied if used):                         ####
####    DIR_JETALGO       : is a suffix to the output folder name (can be used for        ####
####                        extra information)                                            ####
####    PF_CALO_JPT       : chooses the jet type (for PF, akFastPF-files are read         ####
####                        in, see below - does not make a difference when JEC           ####
####                        is overridden)                                                ####
####    MC                : Choose MC-generation                                          ####
####    MC_type           : Choose specific MC-type, here the path to the n-tupels        ####
####                        is overridden by the following datadirmc option               ####
####    BINNING           : choose binning in eta, currently only ""kostas"" and          ####
####                        ""k_HFfix"" are properly defined herekostas                   ####
####    datadirmc         : path to the MC/MC_type n-tupels                               ####
####    DATATYPE          : label for the dataset used to determine the path to the       ####
####                        data n-tupels                                                 ####
####    CORRECT_JETS_L1   : switch to decide whether L1-corrections should be applied     ####
####                        before TwoJetsPtBalanceEvents are created                     ####
####    SINGLEJET         : switch to decide whether single jet triggers/turn-ons         ####
####                        and pt-variables are used.                                    ####
##############################################################################################
##############################################################################################


Example for use of Cmd-line options (8 arguments need to be supplied if used):
  DIR_JETALGO       : 52MC_AK5                                                                                                                              
  PF_CALO_JPT       : JPT                                                                                                                                   
  MC                : Su12                                                                                                                                  
  MC_type           : Z2star52                                                                                                                              
  BINNING           : kostas                                                                                                                                
  datadirmc         : /scratch/hh/current/cms/user/kirschen/2012_Jets_v3/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6_Summer12-PU_S7_START52_V9-v1/merged   
  DATATYPE          : TEST                                                                                                                                  
  CORRECT_JETS_L1   : true                                                                                                                                  
  SINGLEJET         : 0                                                                                                                                     

The line for execution would read:
./scripts/runDiJets_Residuals.py 52MC_AK5 JPT Su12 Z2star52 kostas /scratch/hh/current/cms/user/kirschen/2012_Jets_v3/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6_Summer12-PU_S7_START52_V9-v1/merged TEST true 0                                                                                                                                  

##############################################################################################
##############################################################################################

"""
print Usage


def writeCfg(filename):
    config_1="""
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
SuppDiJet Max cut on L2L3-corrected asymmetry = 0.7 # was 0.4
Et cut on jet              = 0.0
Et cut on gamma            = 0.0
Et cut on Z                = 20.0
Et cut on tower            = 0.0
Et cut on cluster          = 0.0
Et cut on track            = 0.0
Et cut on n+1 Jet          = 9999.0
Min Delta Phi              = 0.0
SuppDiJet Min Delta Phi    = 2.7
Eta max cut on jet         = 5.2 # 4.7 for old JEC# was 5.2
Relative Rest Jet Cut      = 0.2      #NonLeadingJetsEt / PhotonEt
Min had fraction           = -1.05    #Default: 0.07
Max had fraction           = 2.05    #1.05    #Default: 0.95
Et genJet min              = 0.0
Et genJet max              = 7000.0
DeltaR cut on jet matching = 9999
Max cut on relative n+1 Jet Et = 1.0 # was 0.4
SuppDiJet Max cut on relative n+1 Jet Et = 0.4 # was 0.4
Max cut on relative Soft Jet Et = 200.0
#---------------------------------------------------------------------------------
#   Input / Output
#---------------------------------------------------------------------------------

#now later in config-file...  jet correction name   = Spring11_AK5Calo
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
JER - Assert J1J2 In Same Eta Bin = false
create plots                     = true
plots output directory           = plots
#plots format                      = pdf
plots only to root-file = true
export all XY projections = true
HF - Use HF as reference region for TwoJetsPtBalanceEvent = false
HF - Use HF as reference region for TwoJetsPtBalanceEvent = false
#plots only to root-file = true
#export all XY projections = false

# JetTruthEvent plots
create JetTruthEvent plots    =  false
create TwoJetsPtBalanceEvent plots = true
#create TwoJetsPtBalanceEvent PU weighting plots = true


"""

    fcfg = open(filename, "w")
    fcfg.write(config_1)

    binning_values=BinningValues(BINNING,False)
    abs_binning_values=BinningValues(BINNING,True)


    plot_list=['AsymmetryVsPt','AbsAsymmetryVsPt','OneBinAsymmetryVsPt','OneBinAbsAsymmetryVsPt','MPFVsPtAve','AbsMPFVsPtAve','OneBinMPFVsPtAve','OneBinAbsMPFVsPtAve','MPFT1VsPtAve','AbsMPFT1VsPtAve','OneBinMPFT1VsPtAve','OneBinAbsMPFT1VsPtAve','GenAsymmetryVsMeanPt','AbsGenAsymmetryVsMeanPt','GenLeadAsymmetryVsLeadPt','AbsGenLeadAsymmetryVsLeadPt','GenBarrAsymmetryVsBarrPt','AbsGenBarrAsymmetryVsBarrPt','AbsFlavorVsPtAve','AbsMPFVsPtAveFlavor']
    if(BINNING=="JEC_Mikko"):
        cut_list=['30','20','15','10']
        cut_no_list=['.30','.20','.15','.10']
    else:
        cut_list=['40','30','20','10']
        cut_no_list=['.40','.30','.20','.10']
    cut_flav_list=['-0.5','0.5','1.5','2.5','3.5','4.5']
#    cut_list=['40','35','30','25','20','15','10','05']
#    cut_no_list=['.40','.35','.30','.25','.20','.15','.10','.05']

    fcfg.write("TwoJetsPtBalanceEvent plots cut_list =   ")
    for index_cut, cut in enumerate(cut_list):
        fcfg.write(cut + " ")
    fcfg.write("\n")

    fcfg.write("TwoJetsPtBalanceEvent plots cut_no_list =   ")
    for index_cut, cut in enumerate(cut_no_list):
        fcfg.write(cut + " ")
        
    fcfg.write("\n")
        
    

    fcfg.write("TwoJetsPtBalanceEvent plots names =   ")
    for index_cut, cut in enumerate(cut_list):
        for index_plot, plot in enumerate(plot_list):
            fcfg.write(plot + cut + ";")
    fcfg.write("AsymmetryVsMCNPUTruth20_all_eta;AsymmetryVsNPV20_all_eta;AsymmetryVsPt20_all_eta;AsymmetryVsEta;AsymmetryVsEta_all_pt;MPFVsEta_all_pt;AsymmetryVsMCNPUVtx20_all_eta;AsymmetryVsNPV20_pt_bin_all_eta;AsymmetryVsMCNPUVtx20_pt_bin_all_eta;AsymmetryVsMCNPUTruth20_pt_bin_all_eta;AsymmetryVsTJF;AbsAsymmetryVsTJF;AsymmetryVsNPV;AbsAsymmetryVsNPV")
    fcfg.write("\n")
#    for index_samples, samples in enumerate(samples_all):
    for index_cut, cut in enumerate(cut_list):
        #    print "prepare " + samples
	
	##### Asymmetry #####
        fcfg.write(plot_list[0] + cut + " x variable        =  MeanPt; log\n")
        fcfg.write(plot_list[0] + cut + " x edges           =  20 " + RawTrigger_thresholds + "2000\n")
        fcfg.write(plot_list[0] + cut + " y variable        =  Asymmetry\n")
        fcfg.write(plot_list[0] + cut + " y edges           =  31 -0.70 0.70 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3 -0.70 0.70 0.7 1.3\n")
        fcfg.write(plot_list[0] + cut + " bin variable      =  Eta\n")
        fcfg.write(plot_list[0] + cut + " bin edges         =  " + binning_values + "\n")
        fcfg.write(plot_list[0] + cut + " cut variable      =  ThirdJetFractionPlain\n")
        fcfg.write(plot_list[0] + cut + " cut edges         =  0.0 " + cut_no_list[index_cut]+ "\n")
        fcfg.write(plot_list[0] + cut + " correction types  =  L2L3; L2L3Res\n")
        fcfg.write(plot_list[0] + cut + " 1 correction types=  L2L3\n")
        fcfg.write(plot_list[0] + cut + " profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans; IQMean; RatioOfIQMeans\n")
        fcfg.write(plot_list[0] + cut + " input samples     =  0:data;1:MC\n\n")

        fcfg.write(plot_list[1] + cut + " x variable        =  MeanPt; log\n")
        fcfg.write(plot_list[1] + cut + " x edges           =  20 " + RawTrigger_thresholds + "2000\n")
        fcfg.write(plot_list[1] + cut + " y variable        =  Asymmetry\n")
        fcfg.write(plot_list[1] + cut + " y edges           =  31 -0.70 0.70 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3 -0.70 0.70 0.7 1.3 0.0 0.5 0.0 0.5\n")
        fcfg.write(plot_list[1] + cut + " bin variable      =  AbsEta\n")
        fcfg.write(plot_list[1] + cut + " bin edges         =  " + abs_binning_values + "\n")
        fcfg.write(plot_list[1] + cut + " cut variable      =  ThirdJetFractionPlain\n")
        fcfg.write(plot_list[1] + cut + " cut edges         =  0.0 " + cut_no_list[index_cut]+ "\n")
        fcfg.write(plot_list[1] + cut + " correction types  =  L2L3; L2L3Res\n")
        fcfg.write(plot_list[1] + cut + " 1 correction types=  L2L3\n")
        fcfg.write(plot_list[1] + cut + " profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans; IQMean; RatioOfIQMeans; StandardDeviation; GaussFitWidth\n")
        fcfg.write(plot_list[1] + cut + " input samples     =  0:data;1:MC\n\n")

        fcfg.write(plot_list[2] + cut + " x variable        =  MeanPt; log\n")
        fcfg.write(plot_list[2] + cut + " x edges           =  1 20 2000\n")
        fcfg.write(plot_list[2] + cut + " y variable        =  Asymmetry\n")
        fcfg.write(plot_list[2] + cut + " y edges           =  31 -0.70 0.70 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3 -0.70 0.70 0.7 1.3\n")
        fcfg.write(plot_list[2] + cut + " bin variable      =  Eta\n")
        fcfg.write(plot_list[2] + cut + " bin edges         =  " + binning_values + "\n")
        fcfg.write(plot_list[2] + cut + " cut variable      =  ThirdJetFractionPlain\n")
        fcfg.write(plot_list[2] + cut + " cut edges         =  0.0 " + cut_no_list[index_cut]+ "\n")
        fcfg.write(plot_list[2] + cut + " correction types  =  L2L3; L2L3Res\n")
        fcfg.write(plot_list[2] + cut + " 1 correction types=  L2L3\n")
        fcfg.write(plot_list[2] + cut + " profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans; IQMean; RatioOfIQMeans\n")
        fcfg.write(plot_list[2] + cut + " input samples     =  0:data;1:MC\n\n")

        fcfg.write(plot_list[3] + cut + " x variable        =  MeanPt; log\n")
        fcfg.write(plot_list[3] + cut + " x edges           =  1 20 2000\n")
        fcfg.write(plot_list[3] + cut + " y variable        =  Asymmetry\n")
        fcfg.write(plot_list[3] + cut + " y edges           =  31 -0.70 0.70 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3 -0.70 0.70 0.7 1.3\n")
        fcfg.write(plot_list[3] + cut + " bin variable      =  AbsEta\n")
        fcfg.write(plot_list[3] + cut + " bin edges         =  " + abs_binning_values + "\n")
        fcfg.write(plot_list[3] + cut + " cut variable      =  ThirdJetFractionPlain\n")
        fcfg.write(plot_list[3] + cut + " cut edges         =  0.0 " + cut_no_list[index_cut]+ "\n")
        fcfg.write(plot_list[3] + cut + " correction types  =  L2L3; L2L3Res\n")
        fcfg.write(plot_list[3] + cut + " 1 correction types=  L2L3\n")
        fcfg.write(plot_list[3] + cut + " profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans; IQMean; RatioOfIQMeans\n")
        fcfg.write(plot_list[3] + cut + " input samples     =  0:data;1:MC\n\n")
	
	##### MPF raw #####
        fcfg.write(plot_list[4] + cut + " x variable        =  MeanPt; log\n")
        fcfg.write(plot_list[4] + cut + " x edges           =  20 " + RawTrigger_thresholds + "2000\n")
        fcfg.write(plot_list[4] + cut + " y variable        =  MPFResponse\n")
        fcfg.write(plot_list[4] + cut + " y edges           =  31 0.50 1.50 0.5 1.5 0.5 1.5 0.5 1.5\n")
        fcfg.write(plot_list[4] + cut + " bin variable      =  Eta\n")
        fcfg.write(plot_list[4] + cut + " bin edges         =  " + binning_values + "\n")
        fcfg.write(plot_list[4] + cut + " cut variable      =  ThirdJetFractionPlain\n")
        fcfg.write(plot_list[4] + cut + " cut edges         =  0.0 " + cut_no_list[index_cut]+ "\n")
        fcfg.write(plot_list[4] + cut + " correction types  =  L2L3; L2L3Res\n")
        fcfg.write(plot_list[4] + cut + " 1 correction types=  L2L3\n")
        fcfg.write(plot_list[4] + cut + " profile types     =  Mean; GaussFitMean; IQMean;\n")
        fcfg.write(plot_list[4] + cut + " input samples     =  0:data;1:MC\n\n")

        fcfg.write(plot_list[5] + cut + " x variable        =  MeanPt; log\n")
        fcfg.write(plot_list[5] + cut + " x edges           =  20 " + RawTrigger_thresholds + "2000\n")
        fcfg.write(plot_list[5] + cut + " y variable        =  MPFResponse\n")
        fcfg.write(plot_list[5] + cut + " y edges           =  31 0.50 1.50 0.5 1.5 0.5 1.5 0.5 1.5\n")
        fcfg.write(plot_list[5] + cut + " bin variable      =  AbsEta\n")
        fcfg.write(plot_list[5] + cut + " bin edges         =  " + abs_binning_values + "\n")
        fcfg.write(plot_list[5] + cut + " cut variable      =  ThirdJetFractionPlain\n")
        fcfg.write(plot_list[5] + cut + " cut edges         =  0.0 " + cut_no_list[index_cut]+ "\n")
        fcfg.write(plot_list[5] + cut + " correction types  =  L2L3; L2L3Res\n")
        fcfg.write(plot_list[5] + cut + " 1 correction types=  L2L3\n")
        fcfg.write(plot_list[5] + cut + " profile types     =  Mean; GaussFitMean; IQMean;\n")
        fcfg.write(plot_list[5] + cut + " input samples     =  0:data;1:MC\n\n")

        fcfg.write(plot_list[6] + cut + " x variable        =  MeanPt; log\n")
        fcfg.write(plot_list[6] + cut + " x edges           =  1 20 2000\n")
        fcfg.write(plot_list[6] + cut + " y variable        =  MPFResponse\n")
        fcfg.write(plot_list[6] + cut + " y edges           =  31 0.50 1.50 0.5 1.5 0.5 1.5 0.5 1.5\n")
        fcfg.write(plot_list[6] + cut + " bin variable      =  Eta\n")
        fcfg.write(plot_list[6] + cut + " bin edges         =  " + binning_values + "\n")
        fcfg.write(plot_list[6] + cut + " cut variable      =  ThirdJetFractionPlain\n")
        fcfg.write(plot_list[6] + cut + " cut edges         =  0.0 " + cut_no_list[index_cut]+ "\n")
        fcfg.write(plot_list[6] + cut + " correction types  =  L2L3; L2L3Res\n")
        fcfg.write(plot_list[6] + cut + " 1 correction types=  L2L3\n")
        fcfg.write(plot_list[6] + cut + " profile types     =  Mean; GaussFitMean; IQMean;\n")
        fcfg.write(plot_list[6] + cut + " input samples     =  0:data;1:MC\n\n")

        fcfg.write(plot_list[7] + cut + " x variable        =  MeanPt; log\n")
        fcfg.write(plot_list[7] + cut + " x edges           =  1 20 2000\n")
        fcfg.write(plot_list[7] + cut + " y variable        =  MPFResponse\n")
        fcfg.write(plot_list[7] + cut + " y edges           =  31 0.50 1.50 0.5 1.5 0.5 1.5 0.5 1.5\n")
        fcfg.write(plot_list[7] + cut + " bin variable      =  AbsEta\n")
        fcfg.write(plot_list[7] + cut + " bin edges         =  " + abs_binning_values + "\n")
        fcfg.write(plot_list[7] + cut + " cut variable      =  ThirdJetFractionPlain\n")
        fcfg.write(plot_list[7] + cut + " cut edges         =  0.0 " + cut_no_list[index_cut]+ "\n")
        fcfg.write(plot_list[7] + cut + " correction types  =  L2L3; L2L3Res\n")
        fcfg.write(plot_list[7] + cut + " 1 correction types=  L2L3\n")
        fcfg.write(plot_list[7] + cut + " profile types     =  Mean; GaussFitMean; IQMean;\n")
        fcfg.write(plot_list[7] + cut + " input samples     =  0:data;1:MC\n\n")
	
	##### MPF T1 #####
        fcfg.write(plot_list[8] + cut + " x variable        =  MeanPt; log\n")
        fcfg.write(plot_list[8] + cut + " x edges           =  20 " + RawTrigger_thresholds + "2000\n")
        fcfg.write(plot_list[8] + cut + " y variable        =  MPFMETT1Response\n")
        fcfg.write(plot_list[8] + cut + " y edges           =  31 0.50 1.50 0.5 1.5 0.5 1.5 0.5 1.5\n")
        fcfg.write(plot_list[8] + cut + " bin variable      =  Eta\n")
        fcfg.write(plot_list[8] + cut + " bin edges         =  " + binning_values + "\n")
        fcfg.write(plot_list[8] + cut + " cut variable      =  ThirdJetFractionPlain\n")
        fcfg.write(plot_list[8] + cut + " cut edges         =  0.0 " + cut_no_list[index_cut]+ "\n")
        fcfg.write(plot_list[8] + cut + " correction types  =  L2L3; L2L3Res\n")
        fcfg.write(plot_list[8] + cut + " 1 correction types=  L2L3\n")
        fcfg.write(plot_list[8] + cut + " profile types     =  Mean; GaussFitMean; IQMean;\n")
        fcfg.write(plot_list[8] + cut + " input samples     =  0:data;1:MC\n\n")

        fcfg.write(plot_list[9] + cut + " x variable        =  MeanPt; log\n")
        fcfg.write(plot_list[9] + cut + " x edges           =  20 " + RawTrigger_thresholds + "2000\n")
        fcfg.write(plot_list[9] + cut + " y variable        =  MPFMETT1Response\n")
        fcfg.write(plot_list[9] + cut + " y edges           =  31 0.50 1.50 0.5 1.5 0.5 1.5 0.5 1.5\n")
        fcfg.write(plot_list[9] + cut + " bin variable      =  AbsEta\n")
        fcfg.write(plot_list[9] + cut + " bin edges         =  " + abs_binning_values + "\n")
        fcfg.write(plot_list[9] + cut + " cut variable      =  ThirdJetFractionPlain\n")
        fcfg.write(plot_list[9] + cut + " cut edges         =  0.0 " + cut_no_list[index_cut]+ "\n")
        fcfg.write(plot_list[9] + cut + " correction types  =  L2L3; L2L3Res\n")
        fcfg.write(plot_list[9] + cut + " 1 correction types=  L2L3\n")
        fcfg.write(plot_list[9] + cut + " profile types     =  Mean; GaussFitMean; IQMean;\n")
        fcfg.write(plot_list[9] + cut + " input samples     =  0:data;1:MC\n\n")

        fcfg.write(plot_list[10] + cut + " x variable        =  MeanPt; log\n")
        fcfg.write(plot_list[10] + cut + " x edges           =  1 20 2000\n")
        fcfg.write(plot_list[10] + cut + " y variable        =  MPFMETT1Response\n")
        fcfg.write(plot_list[10] + cut + " y edges           =  31 0.50 1.50 0.5 1.5 0.5 1.5 0.5 1.5\n")
        fcfg.write(plot_list[10] + cut + " bin variable      =  Eta\n")
        fcfg.write(plot_list[10] + cut + " bin edges         =  " + binning_values + "\n")
        fcfg.write(plot_list[10] + cut + " cut variable      =  ThirdJetFractionPlain\n")
        fcfg.write(plot_list[10] + cut + " cut edges         =  0.0 " + cut_no_list[index_cut]+ "\n")
        fcfg.write(plot_list[10] + cut + " correction types  =  L2L3; L2L3Res\n")
        fcfg.write(plot_list[10] + cut + " 1 correction types=  L2L3\n")
        fcfg.write(plot_list[10] + cut + " profile types     =  Mean; GaussFitMean; IQMean;\n")
        fcfg.write(plot_list[10] + cut + " input samples     =  0:data;1:MC\n\n")

        fcfg.write(plot_list[11] + cut + " x variable        =  MeanPt; log\n")
        fcfg.write(plot_list[11] + cut + " x edges           =  1 20 2000\n")
        fcfg.write(plot_list[11] + cut + " y variable        =  MPFMETT1Response\n")
        fcfg.write(plot_list[11] + cut + " y edges           =  31 0.50 1.50 0.5 1.5 0.5 1.5 0.5 1.5\n")
        fcfg.write(plot_list[11] + cut + " bin variable      =  AbsEta\n")
        fcfg.write(plot_list[11] + cut + " bin edges         =  " + abs_binning_values + "\n")
        fcfg.write(plot_list[11] + cut + " cut variable      =  ThirdJetFractionPlain\n")
        fcfg.write(plot_list[11] + cut + " cut edges         =  0.0 " + cut_no_list[index_cut]+ "\n")
        fcfg.write(plot_list[11] + cut + " correction types  =  L2L3; L2L3Res\n")
        fcfg.write(plot_list[11] + cut + " 1 correction types=  L2L3\n")
        fcfg.write(plot_list[11] + cut + " profile types     =  Mean; GaussFitMean; IQMean;\n")
        fcfg.write(plot_list[11] + cut + " input samples     =  0:data;1:MC\n\n")

	
	##### Generator reference plots #####
	fcfg.write(plot_list[12] + cut + " x variable        =  MeanPt; log\n")
        fcfg.write(plot_list[12] + cut + " x edges           =  20 " + RawTrigger_thresholds + "2000\n")
        fcfg.write(plot_list[12] + cut + " y variable        =  GenAsymmetry\n")
        fcfg.write(plot_list[12] + cut + " y edges           =  31 -0.70 0.70 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3 -0.70 0.70 0.7 1.3\n")
        fcfg.write(plot_list[12] + cut + " bin variable      =  Eta\n")
        fcfg.write(plot_list[12] + cut + " bin edges         =  " + binning_values + "\n")
        fcfg.write(plot_list[12] + cut + " cut variable      =  ThirdJetFractionPlain\n")
        fcfg.write(plot_list[12] + cut + " cut edges         =  0.0 " + cut_no_list[index_cut]+ "\n")
	fcfg.write(plot_list[12] + cut + " distributions     =  L2L3; L2L3Res\n")
	fcfg.write(plot_list[12] + cut + " 1 distributions   =  L2L3\n")
        fcfg.write(plot_list[12] + cut + " correction types  =  L2L3; L2L3Res\n")
	fcfg.write(plot_list[12] + cut + " 1 correction types=  L2L3\n")
        fcfg.write(plot_list[12] + cut + " profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans; IQMean; RatioOfIQMeans\n")
        fcfg.write(plot_list[12] + cut + " input samples     =  0:data;1:MC\n\n")
	
	fcfg.write(plot_list[13] + cut + " x variable        =  MeanPt; log\n")
        fcfg.write(plot_list[13] + cut + " x edges           =  20 " + RawTrigger_thresholds + "2000\n")
        fcfg.write(plot_list[13] + cut + " y variable        =  GenAsymmetry\n")
        fcfg.write(plot_list[13] + cut + " y edges           =  31 -0.70 0.70 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3 -0.70 0.70 0.7 1.3\n")
        fcfg.write(plot_list[13] + cut + " bin variable      =  AbsEta\n")
        fcfg.write(plot_list[13] + cut + " bin edges         =  " + abs_binning_values + "\n")
        fcfg.write(plot_list[13] + cut + " cut variable      =  ThirdJetFractionPlain\n")
        fcfg.write(plot_list[13] + cut + " cut edges         =  0.0 " + cut_no_list[index_cut]+ "\n")
	fcfg.write(plot_list[13] + cut + " distributions     =  L2L3; L2L3Res\n")
	fcfg.write(plot_list[13] + cut + " 1 distributions   =  L2L3\n")
        fcfg.write(plot_list[13] + cut + " correction types  =  L2L3; L2L3Res\n")
	fcfg.write(plot_list[13] + cut + " 1 correction types=  L2L3\n")
        fcfg.write(plot_list[13] + cut + " profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans; IQMean; RatioOfIQMeans\n")
        fcfg.write(plot_list[13] + cut + " input samples     =  0:data;1:MC\n\n")

	fcfg.write(plot_list[14] + cut + " x variable        =  JetLeadPt; log\n")
        fcfg.write(plot_list[14] + cut + " x edges           =  20 " + RawTrigger_thresholds + "2000\n")
        fcfg.write(plot_list[14] + cut + " y variable        =  GenLeadAsymmetry\n")
        fcfg.write(plot_list[14] + cut + " y edges           =  31 -0.70 0.70 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3 -0.70 0.70 0.7 1.3\n")
        fcfg.write(plot_list[14] + cut + " bin variable      =  Eta\n")
        fcfg.write(plot_list[14] + cut + " bin edges         =  " + binning_values + "\n")
        fcfg.write(plot_list[14] + cut + " cut variable      =  ThirdJetFractionPlain\n")
        fcfg.write(plot_list[14] + cut + " cut edges         =  0.0 " + cut_no_list[index_cut]+ "\n")
	fcfg.write(plot_list[14] + cut + " distributions     =  L2L3; L2L3Res\n")
	fcfg.write(plot_list[14] + cut + " 1 distributions   =  L2L3\n")
        fcfg.write(plot_list[14] + cut + " correction types  =  L2L3; L2L3Res\n")
	fcfg.write(plot_list[14] + cut + " 1 correction types=  L2L3\n")
        fcfg.write(plot_list[14] + cut + " profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans; IQMean; RatioOfIQMeans\n")
        fcfg.write(plot_list[14] + cut + " input samples     =  0:data;1:MC\n\n")
	
	fcfg.write(plot_list[15] + cut + " x variable        =  JetLeadPt; log\n")
        fcfg.write(plot_list[15] + cut + " x edges           =  20 " + RawTrigger_thresholds + "2000\n")
        fcfg.write(plot_list[15] + cut + " y variable        =  GenLeadAsymmetry\n")
        fcfg.write(plot_list[15] + cut + " y edges           =  31 -0.70 0.70 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3 -0.70 0.70 0.7 1.3\n")
        fcfg.write(plot_list[15] + cut + " bin variable      =  AbsEta\n")
        fcfg.write(plot_list[15] + cut + " bin edges         =  " + abs_binning_values + "\n")
        fcfg.write(plot_list[15] + cut + " cut variable      =  ThirdJetFractionPlain\n")
        fcfg.write(plot_list[15] + cut + " cut edges         =  0.0 " + cut_no_list[index_cut]+ "\n")
	fcfg.write(plot_list[15] + cut + " distributions     =  L2L3; L2L3Res\n")
	fcfg.write(plot_list[15] + cut + " 1 distributions   =  L2L3\n")
        fcfg.write(plot_list[15] + cut + " correction types  =  L2L3; L2L3Res\n")
	fcfg.write(plot_list[15] + cut + " 1 correction types=  L2L3\n")
        fcfg.write(plot_list[15] + cut + " profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans; IQMean; RatioOfIQMeans\n")
        fcfg.write(plot_list[15] + cut + " input samples     =  0:data;1:MC\n\n")

	fcfg.write(plot_list[16] + cut + " x variable        =  Jet2Pt; log\n")
        fcfg.write(plot_list[16] + cut + " x edges           =  20 " + RawTrigger_thresholds + "2000\n")
        fcfg.write(plot_list[16] + cut + " y variable        =  GenBarrAsymmetry\n")
        fcfg.write(plot_list[16] + cut + " y edges           =  31 -0.70 0.70 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3 -0.70 0.70 0.7 1.3\n")
        fcfg.write(plot_list[16] + cut + " bin variable      =  Eta\n")
        fcfg.write(plot_list[16] + cut + " bin edges         =  " + binning_values + "\n")
        fcfg.write(plot_list[16] + cut + " cut variable      =  ThirdJetFractionPlain\n")
        fcfg.write(plot_list[16] + cut + " cut edges         =  0.0 " + cut_no_list[index_cut]+ "\n")
	fcfg.write(plot_list[16] + cut + " distributions     =  L2L3; L2L3Res\n")
	fcfg.write(plot_list[16] + cut + " 1 distributions   =  L2L3\n")
        fcfg.write(plot_list[16] + cut + " correction types  =  L2L3; L2L3Res\n")
	fcfg.write(plot_list[16] + cut + " 1 correction types=  L2L3\n")
        fcfg.write(plot_list[16] + cut + " profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans; IQMean; RatioOfIQMeans\n")
        fcfg.write(plot_list[16] + cut + " input samples     =  0:data;1:MC\n\n")
	
	fcfg.write(plot_list[17] + cut + " x variable        =  Jet2Pt; log\n")
        fcfg.write(plot_list[17] + cut + " x edges           =  20 " + RawTrigger_thresholds + "2000\n")
        fcfg.write(plot_list[17] + cut + " y variable        =  GenBarrAsymmetry\n")
        fcfg.write(plot_list[17] + cut + " y edges           =  31 -0.70 0.70 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3 -0.70 0.70 0.7 1.3\n")
        fcfg.write(plot_list[17] + cut + " bin variable      =  AbsEta\n")
        fcfg.write(plot_list[17] + cut + " bin edges         =  " + abs_binning_values + "\n")
        fcfg.write(plot_list[17] + cut + " cut variable      =  ThirdJetFractionPlain\n")
        fcfg.write(plot_list[17] + cut + " cut edges         =  0.0 " + cut_no_list[index_cut]+ "\n")
	fcfg.write(plot_list[17] + cut + " distributions     =  L2L3; L2L3Res\n")
	fcfg.write(plot_list[17] + cut + " 1 distributions   =  L2L3\n")
        fcfg.write(plot_list[17] + cut + " correction types  =  L2L3; L2L3Res\n")
	fcfg.write(plot_list[17] + cut + " 1 correction types=  L2L3\n")
        fcfg.write(plot_list[17] + cut + " profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans; IQMean; RatioOfIQMeans\n")
        fcfg.write(plot_list[17] + cut + " input samples     =  0:data;1:MC\n\n")
	
	##### Flavour plots #####
	
        fcfg.write(plot_list[18] + cut + " x variable        =  MeanPt; log\n")
        fcfg.write(plot_list[18] + cut + " x edges           =  20 " + RawTrigger_thresholds + "2000\n")
        fcfg.write(plot_list[18] + cut + " y variable        =  Flavor\n")
        fcfg.write(plot_list[18] + cut + " y edges           =  5 -1.5 4.5 -1.5 4.5 -1.5 4.5 -1.5 4.5\n")
        fcfg.write(plot_list[18] + cut + " bin variable      =  AbsEta\n")
        fcfg.write(plot_list[18] + cut + " bin edges         =  " + abs_binning_values + "\n")
        fcfg.write(plot_list[18] + cut + " cut variable      =  ThirdJetFractionPlain\n")
        fcfg.write(plot_list[18] + cut + " cut edges         =  0.0 " + cut_no_list[index_cut]+ "\n")
        fcfg.write(plot_list[18] + cut + " correction types  =  L2L3; L2L3Res\n")
        fcfg.write(plot_list[18] + cut + " 1 correction types=  L2L3\n")
        fcfg.write(plot_list[18] + cut + " profile types     =  Mean; GaussFitMean; IQMean\n")		
        fcfg.write(plot_list[18] + cut + " input samples     =  0:data;1:MC\n\n")	
	
        fcfg.write(plot_list[19] + cut + " x variable        =  MeanPt; log\n")
        fcfg.write(plot_list[19] + cut + " x edges           =  20 " + RawTrigger_thresholds + "2000\n")
        fcfg.write(plot_list[19] + cut + " y variable        =  MPFResponse\n")
        fcfg.write(plot_list[19] + cut + " y edges           =  31 0.50 1.50 0.5 1.5 0.5 1.5 0.5 1.5\n")
        fcfg.write(plot_list[19] + cut + " bin variable      =  AbsEta\n")
        fcfg.write(plot_list[19] + cut + " bin edges         =  " + abs_binning_values + "\n")
        fcfg.write(plot_list[19] + cut + " cut variable      =  Flavor\n")
        fcfg.write(plot_list[19] + cut + " cut edges         =  "+ cut_flav_list[index_cut]+cut_flav_list[index_cut+1]+ "\n")
        fcfg.write(plot_list[19] + cut + " correction types  =  L2L3; L2L3Res\n")
        fcfg.write(plot_list[19] + cut + " 1 correction types=  L2L3\n")
        fcfg.write(plot_list[19] + cut + " profile types     =  Mean; GaussFitMean; IQMean\n")		
        fcfg.write(plot_list[19] + cut + " input samples     =  0:data;1:MC\n\n")	
	
        	
	
#    fcfg.write("\n")
#    config_2="""
        config_2="""


AsymmetryVsPt20_all_eta x variable        =  MeanPt; log
AsymmetryVsPt20_all_eta x edges           =  15 20 2000
AsymmetryVsPt20_all_eta y variable        =  Asymmetry
AsymmetryVsPt20_all_eta y edges           =  51 -0.30 0.30 -0.5 0.5
AsymmetryVsPt20_all_eta bin variable      =  Eta
AsymmetryVsPt20_all_eta bin edges         =  -5.192 5.192
AsymmetryVsPt20_all_eta cut variable      =  ThirdJetFractionPlain
AsymmetryVsPt20_all_eta cut edges         =  0.0 0.2
AsymmetryVsPt20_all_eta correction types  =  L2L3; L2L3Res
AsymmetryVsPt20_all_eta 1 correction types  =  L2L3
AsymmetryVsPt20_all_eta profile types     =  Mean
AsymmetryVsPt20_all_eta input samples     =  0:data;1:MC

AsymmetryVsEta x variable        =  Eta
AsymmetryVsEta x edges           =  26 -5.2 5.2
AsymmetryVsEta y variable        =  Asymmetry
AsymmetryVsEta y edges           =  31 -0.70 0.70 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3 -0.70 0.70 0.7 1.3
AsymmetryVsEta bin variable      =  MeanPt
AsymmetryVsEta bin edges         =  20 30 50 80 120 200 360 500 900 7000
AsymmetryVsEta cut variable      =  ThirdJetFractionPlain
AsymmetryVsEta cut edges         =  0.0 0.20
AsymmetryVsEta correction types  =  L2L3; L2L3Res
AsymmetryVsEta 1 correction types  =  L2L3
AsymmetryVsEta profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans; IQMean; RatioOfIQMeans
#AsymmetryVsEta legend label      =  L2L3:CMS default
AsymmetryVsEta input samples     =  0:data;1:MC

AsymmetryVsEta_all_pt x variable        =  Eta
AsymmetryVsEta_all_pt x edges           =  26 -5.2 5.2
AsymmetryVsEta_all_pt y variable        =  Asymmetry
AsymmetryVsEta_all_pt y edges           =  31 -0.70 0.70 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3 -0.70 0.70 0.7 1.3
AsymmetryVsEta_all_pt bin variable      =  MeanPt
AsymmetryVsEta_all_pt bin edges         =  20 7000
AsymmetryVsEta_all_pt cut variable      =  ThirdJetFractionPlain
AsymmetryVsEta_all_pt cut edges         =  0.0 0.20
AsymmetryVsEta_all_pt correction types  =  L2L3; L2L3Res
AsymmetryVsEta_all_pt 1 correction types  =  L2L3
AsymmetryVsEta_all_pt profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans; IQMean; RatioOfIQMeans
AsymmetryVsEta_all_pt input samples     =  0:data;1:MC


MPFVsEta_all_pt x variable        =  Eta
MPFVsEta_all_pt x edges           =  26 -5.2 5.2
MPFVsEta_all_pt y variable        =  MPFResponse
MPFVsEta_all_pt y edges           =  31 0.50 1.50 0.5 1.5 0.5 1.5 0.5 1.5
MPFVsEta_all_pt bin variable      =  Jet2Pt
MPFVsEta_all_pt bin edges         =  20 7000
MPFVsEta_all_pt cut variable      =  ThirdJetFractionPlain
MPFVsEta_all_pt cut edges         =  0.0 0.20
MPFVsEta_all_pt correction types  =  L2L3; L2L3Res
MPFVsEta_all_pt 1 correction types  =  L2L3
MPFVsEta_all_pt profile types     =  Mean; GaussFitMean; IQMean
MPFVsEta_all_pt input samples     =  0:data;1:MC

AsymmetryVsNPV20_all_eta x variable        =  VtxN
AsymmetryVsNPV20_all_eta x edges           =  44 0.0 44.0
AsymmetryVsNPV20_all_eta y variable        =  Asymmetry
AsymmetryVsNPV20_all_eta y edges           =  31 -0.70 0.70 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3 -0.70 0.70 0.7 1.3
AsymmetryVsNPV20_all_eta bin variable      =  Eta
AsymmetryVsNPV20_all_eta bin edges         =  -5.192 5.192
AsymmetryVsNPV20_all_eta cut variable      =  ThirdJetFractionPlain
AsymmetryVsNPV20_all_eta cut edges         =  0.0 0.2
AsymmetryVsNPV20_all_eta correction types  =  L2L3; L2L3Res
AsymmetryVsNPV20_all_eta 1 correction types  =  L2L3
AsymmetryVsNPV20_all_eta profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans; IQMean; RatioOfIQMeans
AsymmetryVsNPV20_all_eta input samples     =  0:data;1:MC


AsymmetryVsMCNPUVtx20_all_eta x variable        =  MCNPUVtx
AsymmetryVsMCNPUVtx20_all_eta x edges           =  44 0.0 44.0
AsymmetryVsMCNPUVtx20_all_eta y variable        =  Asymmetry
AsymmetryVsMCNPUVtx20_all_eta y edges           =  31 -0.70 0.70 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3 -0.70 0.70 0.7 1.3
AsymmetryVsMCNPUVtx20_all_eta bin variable      =  Eta
AsymmetryVsMCNPUVtx20_all_eta bin edges         =  -5.192 5.192
AsymmetryVsMCNPUVtx20_all_eta cut variable      =  ThirdJetFractionPlain
AsymmetryVsMCNPUVtx20_all_eta cut edges         =  0.0 0.2
AsymmetryVsMCNPUVtx20_all_eta correction types  =  L2L3; L2L3Res
AsymmetryVsMCNPUVtx20_all_eta 1 correction types  =  L2L3
AsymmetryVsMCNPUVtx20_all_eta profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans; IQMean; RatioOfIQMeans
AsymmetryVsMCNPUVtx20_all_eta input samples     =  0:data;1:MC

AsymmetryVsMCNPUTruth20_all_eta x variable        =  MCNPUTruth
AsymmetryVsMCNPUTruth20_all_eta x edges           =  44 0.0 44.0
AsymmetryVsMCNPUTruth20_all_eta y variable        =  Asymmetry
AsymmetryVsMCNPUTruth20_all_eta y edges           =  31 -0.70 0.70 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3 -0.70 0.70 0.7 1.3
AsymmetryVsMCNPUTruth20_all_eta bin variable      =  Eta
AsymmetryVsMCNPUTruth20_all_eta bin edges         =  -5.192 5.192
AsymmetryVsMCNPUTruth20_all_eta cut variable      =  ThirdJetFractionPlain
AsymmetryVsMCNPUTruth20_all_eta cut edges         =  0.0 0.2
AsymmetryVsMCNPUTruth20_all_eta correction types  =  L2L3; L2L3Res
AsymmetryVsMCNPUTruth20_all_eta 1 correction types  =  L2L3
AsymmetryVsMCNPUTruth20_all_eta profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans; IQMean; RatioOfIQMeans
AsymmetryVsMCNPUTruth20_all_eta input samples     =  0:data;1:MC


AsymmetryVsNPV20_pt_bin_all_eta x variable        =  VtxN
AsymmetryVsNPV20_pt_bin_all_eta x edges           =  44 0.0 44.0
AsymmetryVsNPV20_pt_bin_all_eta y variable        =  Asymmetry
AsymmetryVsNPV20_pt_bin_all_eta y edges           =  31 -0.70 0.70 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3 -0.70 0.70 0.7 1.3
AsymmetryVsNPV20_pt_bin_all_eta bin variable      =  MeanPt
#AsymmetryVsNPV20_pt_bin_all_eta bin edges         =  20 30 50 80 120 200 360 500 900 7000
AsymmetryVsNPV20_pt_bin_all_eta bin edges         =  45 90 150 210 270 340 420 900 7000
AsymmetryVsNPV20_pt_bin_all_eta cut variable      =  ThirdJetFractionPlain
AsymmetryVsNPV20_pt_bin_all_eta cut edges         =  0.0 0.2
AsymmetryVsNPV20_pt_bin_all_eta correction types  =  L2L3; L2L3Res
AsymmetryVsNPV20_pt_bin_all_eta 1 correction types  =  L2L3
AsymmetryVsNPV20_pt_bin_all_eta profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans; IQMean; RatioOfIQMeans
AsymmetryVsNPV20_pt_bin_all_eta input samples     =  0:data;1:MC


AsymmetryVsMCNPUVtx20_pt_bin_all_eta x variable        =  MCNPUVtx
AsymmetryVsMCNPUVtx20_pt_bin_all_eta x edges           =  44 0.0 44.0
AsymmetryVsMCNPUVtx20_pt_bin_all_eta y variable        =  Asymmetry
AsymmetryVsMCNPUVtx20_pt_bin_all_eta y edges           =  31 -0.70 0.70 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3 -0.70 0.70 0.7 1.3
AsymmetryVsMCNPUVtx20_pt_bin_all_eta bin variable      =  MeanPt
#AsymmetryVsMCNPUVtx20_pt_bin_all_eta bin edges         =  20 30 50 80 120 200 360 500 900 7000
AsymmetryVsMCNPUVtx20_pt_bin_all_eta bin edges         =  45 90 150 210 270 340 420 900 7000
AsymmetryVsMCNPUVtx20_pt_bin_all_eta cut variable      =  ThirdJetFractionPlain
AsymmetryVsMCNPUVtx20_pt_bin_all_eta cut edges         =  0.0 0.2
AsymmetryVsMCNPUVtx20_pt_bin_all_eta correction types  =  L2L3; L2L3Res
AsymmetryVsMCNPUVtx20_pt_bin_all_eta 1 correction types  =  L2L3
AsymmetryVsMCNPUVtx20_pt_bin_all_eta profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans; IQMean; RatioOfIQMeans
AsymmetryVsMCNPUVtx20_pt_bin_all_eta input samples     =  0:data;1:MC

AsymmetryVsMCNPUTruth20_pt_bin_all_eta x variable        =  MCNPUTruth
AsymmetryVsMCNPUTruth20_pt_bin_all_eta x edges           =  44 0.0 44.0
AsymmetryVsMCNPUTruth20_pt_bin_all_eta y variable        =  Asymmetry
AsymmetryVsMCNPUTruth20_pt_bin_all_eta y edges           =  31 -0.70 0.70 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3 -0.70 0.70 0.7 1.3
AsymmetryVsMCNPUTruth20_pt_bin_all_eta bin variable      =  MeanPt
#AsymmetryVsMCNPUTruth20_pt_bin_all_eta bin edges         =  20 30 50 80 120 200 360 500 900 7000
AsymmetryVsMCNPUTruth20_pt_bin_all_eta bin edges         =  45 90 150 210 270 340 420 900 7000
AsymmetryVsMCNPUTruth20_pt_bin_all_eta cut variable      =  ThirdJetFractionPlain
AsymmetryVsMCNPUTruth20_pt_bin_all_eta cut edges         =  0.0 0.2
AsymmetryVsMCNPUTruth20_pt_bin_all_eta correction types  =  L2L3; L2L3Res
AsymmetryVsMCNPUTruth20_pt_bin_all_eta 1 correction types  =  L2L3
AsymmetryVsMCNPUTruth20_pt_bin_all_eta profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans; IQMean; RatioOfIQMeans
AsymmetryVsMCNPUTruth20_pt_bin_all_eta input samples     =  0:data;1:MC


"""

    fcfg.write(config_2)


    fcfg.write("AsymmetryVsTJF x variable        =  ThirdJetFractionPlain\n")
    fcfg.write("AsymmetryVsTJF x edges           =  8 0.00 0.40\n")
    fcfg.write("AsymmetryVsTJF y variable        =  Asymmetry\n")
    fcfg.write("AsymmetryVsTJF y edges           =  31 -0.70 0.70 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3 -0.70 0.70 0.7 1.3\n")
    fcfg.write("AsymmetryVsTJF bin variable      =  Eta\n")
    fcfg.write("AsymmetryVsTJF bin edges         =  " + binning_values + "\n")
    fcfg.write("AsymmetryVsTJF correction types  =  L2L3; L2L3Res\n")
    fcfg.write("AsymmetryVsTJF 1 correction types  =  L2L3\n")
    fcfg.write("AsymmetryVsTJF profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans; IQMean; RatioOfIQMeans\n")
    fcfg.write("#AsymmetryVTJF legend label      =  L2L3:CMS default\n")
    fcfg.write("AsymmetryVsTJF input samples     =  0:data;1:MC\n")
    fcfg.write("\n\n")


    fcfg.write("AbsAsymmetryVsTJF x variable        =  ThirdJetFractionPlain\n")
    fcfg.write("AbsAsymmetryVsTJF x edges           =  8 0.00 0.40\n")
    fcfg.write("AbsAsymmetryVsTJF y variable        =  Asymmetry\n")
    fcfg.write("AbsAsymmetryVsTJF y edges           =  31 -0.70 0.70 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3 -0.70 0.70 0.7 1.3\n")
    fcfg.write("AbsAsymmetryVsTJF bin variable      =  AbsEta\n")
    fcfg.write("AbsAsymmetryVsTJF bin edges         =  " + abs_binning_values + "\n")
    fcfg.write("AbsAsymmetryVsTJF correction types  =  L2L3; L2L3Res\n")
    fcfg.write("AbsAsymmetryVsTJF 1 correction types  =  L2L3\n")
    fcfg.write("AbsAsymmetryVsTJF profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans; IQMean; RatioOfIQMeans\n")
    fcfg.write("#AsymmetryVTJF legend label      =  L2L3:CMS default\n")
    fcfg.write("AbsAsymmetryVsTJF input samples     =  0:data;1:MC\n")
    fcfg.write("\n\n")



    fcfg.write("AsymmetryVsNPV x variable        =  VtxN\n")
    fcfg.write("AsymmetryVsNPV x edges           =  6 0.0 12.0\n")
    fcfg.write("AsymmetryVsNPV y variable        =  Asymmetry\n")
    fcfg.write("AsymmetryVsNPV y edges           =  31 -0.70 0.70 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3 -0.70 0.70 0.7 1.3\n")
    fcfg.write("AsymmetryVsNPV bin variable      =  Eta\n")
    fcfg.write("AsymmetryVsNPV bin edges         =  " + binning_values + "\n")
    fcfg.write("AsymmetryVsNPV cut variable      =  ThirdJetFractionPlain\n")
    fcfg.write("AsymmetryVsNPV cut edges         =  0.0 0.20\n")
    fcfg.write("AsymmetryVsNPV correction types  =  L2L3; L2L3Res\n")
    fcfg.write("AsymmetryVsNPV 1 correction types  =  L2L3\n")
    fcfg.write("AsymmetryVsNPV profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans; IQMean; RatioOfIQMeans\n")
    fcfg.write("#AsymmetryVNPV legend label      =  L2L3:CMS default\n")
    fcfg.write("AsymmetryVsNPV input samples     =  0:data;1:MC\n")
    fcfg.write("\n\n")

    fcfg.write("AbsAsymmetryVsNPV x variable        =  VtxN\n")
    fcfg.write("AbsAsymmetryVsNPV x edges           =  6 0.0 12.0\n")
    fcfg.write("AbsAsymmetryVsNPV y variable        =  Asymmetry\n")
    fcfg.write("AbsAsymmetryVsNPV y edges           =  31 -0.70 0.70 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3 -0.70 0.70 0.7 1.3\n")
    fcfg.write("AbsAsymmetryVsNPV bin variable      =  AbsEta\n")
    fcfg.write("AbsAsymmetryVsNPV bin edges         =  " + abs_binning_values + "\n")
    fcfg.write("AbsAsymmetryVsNPV cut variable      =  ThirdJetFractionPlain\n")
    fcfg.write("AbsAsymmetryVsNPV cut edges         =  0.0 0.20\n")
    fcfg.write("AbsAsymmetryVsNPV correction types  =  L2L3; L2L3Res\n")
    fcfg.write("AbsAsymmetryVsNPV 1 correction types  =  L2L3\n")
    fcfg.write("AbsAsymmetryVsNPV profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans; IQMean; RatioOfIQMeans\n")
    fcfg.write("#AbsAsymmetryVNPV legend label      =  L2L3:CMS default\n")
    fcfg.write("AbsAsymmetryVsNPV input samples     =  0:data;1:MC\n")
    fcfg.write("\n\n")





    if(CORRECTION!="ntuple"):
        fcfg.write("jet correction source = JetMETCor\n")
        fcfg.write("jet correction name   = "+CORRECTIONS+"\n")
        fcfg.write("MC jet correction name   = "+CORRECTIONS+"_MC \n")
    fcfg.write("Di-Jet input file = dijetlist\n")
    fcfg.write("Di-Jet Control1 input file = mcdijetlist\n")
    fcfg.write("Output file       = "+output+"\n");
    fcfg.write("Number of Threads = "+str(nthreads)+"\n")
#    fcfg.write("Number of IO Threads = 5\n")
    fcfg.write("Number of IO Threads = "+str(niothreads)+"\n")
    fcfg.write("Min cut on run number = "+ str(MinRunNumber) +"\n")
    fcfg.write("Max cut on run number = "+ str(MaxRunNumber) +"\n")



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

    if(SINGLEJET==1):
        fcfg.write("Use single jet triggers = true\n")
#    trigger_names = TriggerNamesThresholds(DATAYEAR,USE_NEW_TRIGGERS_AND_FASTPF,SINGLEJET,jettype,"names")
    fcfg.write(trigger_names)
#    trigger_thresholds = TriggerNamesThresholds(DATAYEAR,USE_NEW_TRIGGERS_AND_FASTPF,SINGLEJET,jettype,"thresholds")
    fcfg.write(trigger_thresholds)
        
    if(input != ""):
        fcfg.write("input calibration = Kalibri; "+input+"\n");

    fcfg.write("correct jets L1 = "+ CORRECT_JETS_L1 + "\n");
    fcfg.write("correct JEC and scale JER = "+ SMEAR_JER + "\n");
    fcfg.write("correct JEC and scale JER name = "+ configureJERsmearing(jettype) + "\n");
    fcfg.write("fire all triggers = " + MC_fire_all_triggers +"\n");
    fcfg.write("DiJetEventCuts = true\n");
    fcfg.write("EventWeightProcessor = true\n");
    fcfg.write("EventBinning = false\n");
    fcfg.write("DiJetEventWeighting = true\n");
    fcfg.write("PUTruthReweighting = true\n");
    fcfg.write("PU weighting = false   \n");

    PU_weighting_info = PUWeightingInfo(DATATYPE,MC_type)
    fcfg.write(PU_weighting_info);
    
    fcfg.write("MAX n reco Vtx             = " + str(nMaxRecoVtx) +"\n");
    fcfg.write("MIN n reco Vtx             = " + str(nMinRecoVtx) +"\n");
    fcfg.write("set weights to one         = " + str(weights_eq_one) +"\n\n\n");

    fcfg.write("output dir name         = " + dirname +"\n");

    fcfg.close()
    return




#kostas binning:
#0: 0.0 1: 0.261 2: 0.522 3: 0.783 4: 0.957 5: 1.131 6: 1.305 7: 1.479
#8: 1.93 9: 2.322 10: 2.411 11: 2.5 12: 2.853 13: 2.964 14: 3.139 15: 3.489 5.191







##################################
## one suggestion would be to replace all those configs set here and modified by cmd-line options to "real" config-files
## an example is created by writeconfigs.py and called SampleConfigRunDiJets.cfg. It could be expanded and adapted to any needs.
##################################

config = ConfigParser.SafeConfigParser()
config.read('/afs/naf.desy.de/user/k/kirschen/public/SampleConfigRunDiJets.cfg')

if len(sys.argv) > 1:
    print "there are", len(sys.argv)-1, "arguments:"
    for arg in sys.argv[1:]:
        print arg
else:
    print "there are no arguments!"

MAINSECTIONTOREAD="DEFAULT"
EXTRASECTIONTOREAD="PFJets"
if len(sys.argv) > 3:
    MAINSECTIONTOREAD=sys.argv[1]
    EXTRASECTIONTOREAD=sys.argv[2]

#deactivated in order not to interfere with below config

#MC                     = config.get(SECTIONTOREAD, 'MC')
#DO_MC_ONLY_STUDY       = config.get(MAINSECTIONTOREAD, 'DO_MC_ONLY_STUDY')
DO_MC_ONLY_STUDY       = "false"
#PF_CALO_JPT       = config.get(EXTRAMAINSECTIONTOREAD, 'PF_CALO_JPT')
##################################
## end of dummy part as suggestion
##################################



##################################
##################################
##################################
# main program starts here!!!
# change these variables to steer the process
# start in CalibCore root folder
# with ./scripts/runDiJets_Residuals.py
# or using starter script for batch-processing like
# run_batch_dijetresiduals.py
##################################
##################################
##################################
##################################
## is a suffix to the output folder name (can be used for extra information)
##################################
DIR_JETALGO="RD-V2-TESTING"
##################################
## chooses the jet type (for PF, akFastPF-files are read in, see below - does not make a difference when JEC is overridden)
##################################
PF_CALO_JPT="PFCHS"
##################################
## chooses the jet algorithm - used to pick the corresponding n-tupel .root-files (ak5 is default)
##################################
jetalgo="ak5"
##################################
## Override JEC from text files as defined in JetMETCorFactorsFactory.cc; set to "ntuple" to use n-tuple corrections
##################################
CORRECTION="2013Winter14V1_AK5"
##################################
## Switch to decide whether L1 corrections should be applied or not (default definitely "true" in 2012 ;) )
##################################
CORRECT_JETS_L1="true"
##################################
## Switch to decide whether JER should be smeared in control events (MC)
##################################
SMEAR_JER="true"
##################################
## DATAYEAR variable used to determine trigger thresholds, datasets, ...
##################################
DATAYEAR = "2013"
##################################
## Detailled datasample, similar influence as above
##################################
DATATYPE = "2013ABCD_ReReco"
##################################
## choose binning in eta, currently only "kostas" and "k_HFfix" are properly defined here
##################################
BINNING="Fine"
##################################
## Use single jet triggers if =1 (influences trigger thresholds and trigger pt variables in runtime, look for useSingleJetTriggers_ in code)
##################################
SINGLEJET=0
##################################
## Choose MC-generation 
##################################
MC = "Su12"
##################################
## Choose specific MC-type, determines where to look for n-tupels to read in
##################################
MC_type="Z253_RunDependent"#"Z253_V11_T1T2"#_PUS6S7"
##################################
## Choose minimum run number to read in, important for 2011 dataset, where MinRunNumber=163337 in order to get debugged corrected pt dijetave-triggers
##################################
MinRunNumber=1#193117
##################################
## Choose maximum run number
##################################
MaxRunNumber=300000
##################################
## Use with care - was used to compare MC with MC that had inconsistent weights, sets weights of "data" equal to 1 if set to true
##################################
weights_eq_one="false"
##################################
## Use with care - was used to compare MC with MC that had inconsistent trigger information, currently only implemented for dijetavetriggers (corrected and uncorrected triggers), would need update in DiJetReader.cc
##################################
MC_fire_all_triggers="false"
##################################
## allows to cut on number of reconstructed vertices (e.g. in order to compare a low and high PU sample)
##################################
nMaxRecoVtx =500000
if(nMaxRecoVtx < 100):
    DIR_JETALGO = DIR_JETALGO + "_PU_" + str(nMaxRecoVtx) 
##################################
## same as above, but cut on min number
##################################
nMinRecoVtx =0
if(nMinRecoVtx > 0):
    DIR_JETALGO = DIR_JETALGO + "_PUmin_" + str(nMinRecoVtx) 


##################################
#more things to edit/influence:
#
# - whether to save plots in eps as well
# - whether to save x/y-projections of each 2D-histo
# - plot_list and cut_list 
# - whether eventprocessors are active or not
# - whether plots are created at the intermediate (eventprocessor steps)
# - nevents to e.g. run over only a few events
# - MAIN_dirname as it is hardcoded "root-dir" of plotting and exporting
# - DATATYPES_NEW_TRIGGER - right now needs to add datatypes here for reading in ak5FastPF-files...
# ...much more
##################################


if(DO_MC_ONLY_STUDY=="true"):
    MinRunNumber=-1

    #    MC_fire_all_triggers="true"
    DATATYPE="Z2wPU_DMC"

    

USE_NEW_TRIGGERS_AND_FASTPF=0

DATATYPES_NEW_TRIGGER=importDatatypesNewTrigger()

for new_trigger in DATATYPES_NEW_TRIGGER:
    if(DATATYPE==new_trigger):
        USE_NEW_TRIGGERS_AND_FASTPF=1


#keep following cmd-line option part for compatibility

if len(sys.argv) > 1:
    print "there are", len(sys.argv)-1, "arguments:"
    for arg in sys.argv[1:]:
        print arg
else:
    print "there are no arguments!"


if (len(sys.argv) > 8):
    print "redefine parameters from cmdline-options... "
    DIR_JETALGO=sys.argv[1]
    PF_CALO_JPT=sys.argv[2]
    MC=sys.argv[3]
    MC_type=sys.argv[4]
    BINNING=sys.argv[5]
    DATATYPE=sys.argv[7]
    CORRECT_JETS_L1=sys.argv[8]
    SINGLEJET=int(sys.argv[9])
    print "DIR_JETALGO="+ DIR_JETALGO + " PF_CALO_JPT="+PF_CALO_JPT+ " MC="+MC+" MC-type=" + MC_type+" BINNING="+BINNING+" DATATYPE="+DATATYPE+" CORRECT_JETS_L1="+CORRECT_JETS_L1+" SINGLEJET="+str(SINGLEJET)
    print "... done."
#   datadirmc=sys.argv[6] #is done below


jettype = jetalgo+PF_CALO_JPT
jettype_import=jettype
CORRECTIONSUFFIX=PF_CALO_JPT
if(USE_NEW_TRIGGERS_AND_FASTPF or DATATYPE=="42X_uncorr"):
    if(PF_CALO_JPT=="PF" or PF_CALO_JPT=="PFCHS"):
      if(jetalgo=="ak5"):
        CORRECTIONSUFFIX="Fast"+PF_CALO_JPT
        if(PF_CALO_JPT=="PF"):
            jettype_import=jetalgo+"Fast"+PF_CALO_JPT
#    jettype_import=jetalgo+"Fast"+PF_CALO_JPT


CORRECTIONS=CORRECTION+CORRECTIONSUFFIX

datadir   = determineDataDir(DATAYEAR,DATATYPE)
datadirmc = determineDataDirMC(MC,MC_type)


trigger_thresholds = TriggerNamesThresholds(DATAYEAR,USE_NEW_TRIGGERS_AND_FASTPF,SINGLEJET,jettype,"thresholds")
trigger_names = TriggerNamesThresholds(DATAYEAR,USE_NEW_TRIGGERS_AND_FASTPF,SINGLEJET,jettype,"names")
RawTrigger_thresholds = TriggerNamesThresholds(DATAYEAR,USE_NEW_TRIGGERS_AND_FASTPF,SINGLEJET,jettype,"RawThresholds")


if (len(sys.argv) > 6):
    print "even override datadirmc from cmdline-options... "
    datadirmc=sys.argv[6]
    print "datadirmc="+ datadirmc 
    print "... done."




nthreads = 21
niothreads = 7
nevents = 10000
#nthreads = 1
#niothreads = 1
MAIN_dirname = "/nfs/dust/cms/user/rathjd/KalibriNew/Kalibri/"+DATAYEAR+DATATYPE+"_CORR" + CORRECTION +"_MC_"+MC+MC_type+"_kostas_"+ DIR_JETALGO
dirname = MAIN_dirname + "/dijetsFall10_TuneZ2_AK5"+PF_CALO_JPT+"_weighted_residuals_"+BINNING
useconstraint = False
batch = False
doBinnedFit = False
doUnbinnedFit = True




if not os.path.exists(MAIN_dirname):
    os.system("mkdir "+MAIN_dirname)

if os.path.exists(dirname):
    os.system("rm -rf "+dirname+"/*")
else:
    os.system("mkdir "+dirname)



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
os.system("ls "+datadir+"/*"+jettype_import+"*.root > "+dirname+"/dijetlist");
os.system("ls "+datadirmc+"/*"+jettype_import+"*.root > "+dirname+"/mcdijetlist");

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

