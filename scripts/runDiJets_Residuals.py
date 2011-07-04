#! /usr/bin/python

import os

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
Et cut on jet              = 0.0
Et cut on gamma            = 0.0
Et cut on Z                = 20.0
Et cut on tower            = 0.0
Et cut on cluster          = 0.0
Et cut on track            = 0.0
Et cut on n+1 Jet          = 9999.0
Min Delta Phi              = 2.7
Eta max cut on jet         = 5.2  
Relative Rest Jet Cut      = 0.2      #NonLeadingJetsEt / PhotonEt
Min had fraction           = -1.05    #Default: 0.07
Max had fraction           = 1.05    #Default: 0.95
Et genJet min              = 0.0
Et genJet max              = 7000.0
DeltaR cut on jet matching = 9999
Max cut on relative n+1 Jet Et = 200.0 # was 0.4
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
create plots                     = true
plots output directory           = plots
#plots format                      = pdf
plots only to root-file = true

# JetTruthEvent plots
create JetTruthEvent plots    =  false
create TwoJetsPtBalanceEvent plots = true


"""

    fcfg = open(filename, "w")
    fcfg.write(config_1)

    if(BINNING=="kostas"):
        binning_values = "-5.191 -3.489 -3.139 -2.964 -2.853 -2.5 -2.411 -2.322 -1.93 -1.479 -1.305 -1.131 -0.957 -0.783 -0.522 -0.261 0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191"
        abs_binning_values = "0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191"
    if(BINNING=="k_HFfix"):
        binning_values = "-5.191 -2.964 -2.853 -2.5 -2.411 -2.322 -1.93 -1.479 -1.305 -1.131 -0.957 -0.783 -0.522 -0.261 0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 5.191"
        abs_binning_values = "0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 5.191"


    plot_list=['AsymmetryVsPt','AbsAsymmetryVsPt','OneBinAsymmetryVsPt','OneBinAbsAsymmetryVsPt']
    cut_list=['40','30','20','10']
    cut_no_list=['.40','.30','.20','.10']
#    cut_list=['40','35','30','25','20','15','10','05']
#    cut_no_list=['.40','.35','.30','.25','.20','.15','.10','.05']
        
    fcfg.write("TwoJetsPtBalanceEvent plots names =   ")
    for index_cut, cut in enumerate(cut_list):
        fcfg.write(plot_list[0] + cut + ";")
        fcfg.write(plot_list[1] + cut + ";")
        fcfg.write(plot_list[2] + cut + ";")
        fcfg.write(plot_list[3] + cut + ";")
    fcfg.write("AsymmetryVsPt20_all_eta;AsymmetryVsEta;AsymmetryVsTJF;AbsAsymmetryVsTJF")
    fcfg.write("\n")
##    for index_samples, samples in enumerate(samples_all):
    for index_cut, cut in enumerate(cut_list):
        #    print "prepare " + samples
        fcfg.write(plot_list[0] + cut + " x variable        =   MeanPt; log\n")
        fcfg.write(plot_list[0] + cut + " x edges           =  15 20 2000\n")
        fcfg.write(plot_list[0] + cut + " y variable        =  Asymmetry\n")
        fcfg.write(plot_list[0] + cut + " y edges           =  51 -0.30 0.30 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3\n")
        fcfg.write(plot_list[0] + cut + " bin variable      =  Eta\n")
        fcfg.write(plot_list[0] + cut + " bin edges         =  " + binning_values + "\n")
        fcfg.write(plot_list[0] + cut + " cut variable      =  ThirdJetFractionPlain\n")
        fcfg.write(plot_list[0] + cut + " cut edges         = 0.0 " + cut_no_list[index_cut]+ "\n")
        fcfg.write(plot_list[0] + cut + " correction types  =  L2L3; L2L3Res\n")
        fcfg.write(plot_list[0] + cut + " 1 correction types  =  L2L3\n")
        fcfg.write(plot_list[0] + cut + " profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans\n")
        fcfg.write(plot_list[0] + cut + " input samples     =  0:data;1:MC\n\n")

        fcfg.write(plot_list[1] + cut + " x variable        =   MeanPt; log\n")
        fcfg.write(plot_list[1] + cut + " x edges           =  15 20 2000\n")
        fcfg.write(plot_list[1] + cut + " y variable        =  Asymmetry\n")
        fcfg.write(plot_list[1] + cut + " y edges           =  51 -0.30 0.30 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3\n")
        fcfg.write(plot_list[1] + cut + " bin variable      =  AbsEta\n")
        fcfg.write(plot_list[1] + cut + " bin edges         =  " + abs_binning_values + "\n")
        fcfg.write(plot_list[1] + cut + " cut variable      =  ThirdJetFractionPlain\n")
        fcfg.write(plot_list[1] + cut + " cut edges         = 0.0 " + cut_no_list[index_cut]+ "\n")
        fcfg.write(plot_list[1] + cut + " correction types  =  L2L3; L2L3Res\n")
        fcfg.write(plot_list[1] + cut + " 1 correction types  =  L2L3\n")
        fcfg.write(plot_list[1] + cut + " profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans\n")
        fcfg.write(plot_list[1] + cut + " input samples     =  0:data;1:MC\n\n")

        fcfg.write(plot_list[2] + cut + " x variable        =   MeanPt; log\n")
        fcfg.write(plot_list[2] + cut + " x edges           =  1 20 2000\n")
        fcfg.write(plot_list[2] + cut + " y variable        =  Asymmetry\n")
        fcfg.write(plot_list[2] + cut + " y edges           =  51 -0.30 0.30 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3\n")
        fcfg.write(plot_list[2] + cut + " bin variable      =  Eta\n")
        fcfg.write(plot_list[2] + cut + " bin edges         =  " + binning_values + "\n")
        fcfg.write(plot_list[2] + cut + " cut variable      =  ThirdJetFractionPlain\n")
        fcfg.write(plot_list[2] + cut + " cut edges         = 0.0 " + cut_no_list[index_cut]+ "\n")
        fcfg.write(plot_list[2] + cut + " correction types  =  L2L3; L2L3Res\n")
        fcfg.write(plot_list[2] + cut + " 1 correction types  =  L2L3\n")
        fcfg.write(plot_list[2] + cut + " profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans\n")
        fcfg.write(plot_list[2] + cut + " input samples     =  0:data;1:MC\n\n")

        fcfg.write(plot_list[3] + cut + " x variable        =   MeanPt; log\n")
        fcfg.write(plot_list[3] + cut + " x edges           =  1 20 2000\n")
        fcfg.write(plot_list[3] + cut + " y variable        =  Asymmetry\n")
        fcfg.write(plot_list[3] + cut + " y edges           =  51 -0.30 0.30 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3\n")
        fcfg.write(plot_list[3] + cut + " bin variable      =  AbsEta\n")
        fcfg.write(plot_list[3] + cut + " bin edges         =  " + abs_binning_values + "\n")
        fcfg.write(plot_list[3] + cut + " cut variable      =  ThirdJetFractionPlain\n")
        fcfg.write(plot_list[3] + cut + " cut edges         = 0.0 " + cut_no_list[index_cut]+ "\n")
        fcfg.write(plot_list[3] + cut + " correction types  =  L2L3; L2L3Res\n")
        fcfg.write(plot_list[3] + cut + " 1 correction types  =  L2L3\n")
        fcfg.write(plot_list[3] + cut + " profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans\n")
        fcfg.write(plot_list[3] + cut + " input samples     =  0:data;1:MC\n\n")


        config_2="""

AsymmetryVsPt20_all_eta x variable        =  MeanPt; log
AsymmetryVsPt20_all_eta x edges           =  15 20 2000
AsymmetryVsPt20_all_eta y variable        =  Asymmetry
AsymmetryVsPt20_all_eta y edges           =  51 -0.30 0.30 -0.5 0.5
AsymmetryVsPt20_all_eta cut variable      =  ThirdJetFractionPlain
AsymmetryVsPt20_all_eta cut edges         =  0.0 0.2
AsymmetryVsPt20_all_eta correction types  =  L2L3; L2L3Res
AsymmetryVsPt20_all_eta 1 correction types  =  L2L3
AsymmetryVsPt20_all_eta profile types     =  Mean
AsymmetryVsPt20_all_eta input samples     =  0:data;1:MC

AsymmetryVsEta x variable        =  Eta
AsymmetryVsEta x edges           =  26 -5.2 5.2
AsymmetryVsEta y variable        =  Asymmetry
AsymmetryVsEta y edges           =  51 -0.30 0.30 -0.5 0.8 -0.5 0.8 0.7 1.3 0.7 1.3
AsymmetryVsEta bin variable      =  MeanPt
AsymmetryVsEta bin edges         =  20 30 50 80 120 200 360 500 900 7000
AsymmetryVsEta cut variable      =  ThirdJetFractionPlain
AsymmetryVsEta cut edges         =  0.0 0.20
AsymmetryVsEta correction types  =  L2L3; L2L3Res
AsymmetryVsEta 1 correction types  =  L2L3
AsymmetryVsEta profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AsymmetryVsEta legend label      =  L2L3:CMS default
AsymmetryVsEta input samples     =  0:data;1:MC


"""

    fcfg.write(config_2)


    fcfg.write("AsymmetryVsTJF x variable        =  ThirdJetFractionPlain\n")
    fcfg.write("AsymmetryVsTJF x edges           =  8 0.00 0.40\n")
    fcfg.write("AsymmetryVsTJF y variable        =  Asymmetry\n")
    fcfg.write("AsymmetryVsTJF y edges           =  51 -0.30 0.30 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3\n")
    fcfg.write("AsymmetryVsTJF bin variable      =  Eta\n")
    fcfg.write("AsymmetryVsTJF bin edges         =  " + binning_values + "\n")
    fcfg.write("AsymmetryVsTJF correction types  =  L2L3; L2L3Res\n")
    fcfg.write("AsymmetryVsTJF 1 correction types  =  L2L3\n")
    fcfg.write("AsymmetryVsTJF profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans\n")
    fcfg.write("#AsymmetryVTJF legend label      =  L2L3:CMS default\n")
    fcfg.write("AsymmetryVsTJF input samples     =  0:data;1:MC\n")
    fcfg.write("\n\n")


    fcfg.write("AbsAsymmetryVsTJF x variable        =  ThirdJetFractionPlain\n")
    fcfg.write("AbsAsymmetryVsTJF x edges           =  8 0.00 0.40\n")
    fcfg.write("AbsAsymmetryVsTJF y variable        =  Asymmetry\n")
    fcfg.write("AbsAsymmetryVsTJF y edges           =  51 -0.30 0.30 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3\n")
    fcfg.write("AbsAsymmetryVsTJF bin variable      =  AbsEta\n")
    fcfg.write("AbsAsymmetryVsTJF bin edges         =  " + abs_binning_values + "\n")
    fcfg.write("AbsAsymmetryVsTJF correction types  =  L2L3; L2L3Res\n")
    fcfg.write("AbsAsymmetryVsTJF 1 correction types  =  L2L3\n")
    fcfg.write("AbsAsymmetryVsTJF profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans\n")
    fcfg.write("#AsymmetryVTJF legend label      =  L2L3:CMS default\n")
    fcfg.write("AbsAsymmetryVsTJF input samples     =  0:data;1:MC\n")
    fcfg.write("\n\n")






    if(CORRECTION!="ntuple"):
        fcfg.write("jet correction source = JetMETCor\n")
        fcfg.write("jet correction name   = "+CORRECTIONS+"\n")
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

    if(DATAYEAR == "2011"):
        if(DATAYEAR == "2011" and DATATYPE=="PrRe62pb" or DATATYPE=="42X_corr" or DATATYPE=="42X_PrRe" or DATATYPE=="42X_combPrRe_ReRe"):
            fcfg.write("Di-Jet trigger names = HLT_DiJetAve30;HLT_DiJetAve60;HLT_DiJetAve80;HLT_DiJetAve110;HLT_DiJetAve150;HLT_DiJetAve190;HLT_DiJetAve240;HLT_DiJetAve300;HLT_DiJetAve370\n")
#conservative thresholds....
#            if(jettype == "ak5Calo"):
#                fcfg.write("Di-Jet trigger thresholds = 45 85 105 130 175 220 270 335 405\n") #Matthias preliminary
#            if(jettype == "ak5PF"):
#                fcfg.write("Di-Jet trigger thresholds = 45 85 105 130 175 220 270 335 405\n") #Matthias preliminary

#new default thresholds....
            if(jettype == "ak5Calo"):
                fcfg.write("Di-Jet trigger thresholds = 35 69 89 120 163 204 256 318 390\n") #Matthias preliminary
            if(jettype == "ak5PF"):
                fcfg.write("Di-Jet trigger thresholds = 40 75 100 135 175 220 273 335 405 \n") #Matthias preliminary
            if(jettype == "ak5JPT"):
                fcfg.write("Di-Jet trigger thresholds = 40 75 100 135 175 220 273 335 405 \n") #Matthias preliminary

##aggressive thresholds old
#            if(jettype == "ak5Calo"):
#                fcfg.write("Di-Jet trigger thresholds = 31 62 81 111 156 199 251 313 385\n") #Matthias preliminary
#            if(jettype == "ak5PF"):
#                fcfg.write("Di-Jet trigger thresholds = 40 80 104 130 172 216 269 333 405 \n") #Matthias preliminary
        else:
            fcfg.write("Di-Jet trigger names = HLT_DiJetAve15U;HLT_DiJetAve30U;HLT_DiJetAve50U;HLT_DiJetAve70U;HLT_DiJetAve100U;HLT_DiJetAve140U;HLT_DiJetAve180U;HLT_DiJetAve300U\n")
            if(jettype == "ak5Calo"):
                fcfg.write("Di-Jet trigger thresholds = 38 59 86 111 147 196 249 389\n") #Matthias
            if(jettype == "ak5PF"):
                fcfg.write("Di-Jet trigger thresholds = 43 70 100 127 168 214 279 423 \n") #Matthias
            if(jettype == "ak5JPT"):
                fcfg.write("Di-Jet trigger thresholds = 43 66 96 124 165 220 285 430 \n") #Extrapolation from PF  for last two bins Matthias
            if(jettype == "ak7PF"):
                fcfg.write("Di-Jet trigger thresholds = 47 79 111 140 187 240\n") # CMS AN-2010/371 not updated for 2011
            if(jettype == "ak7Calo"):
                fcfg.write("Di-Jet trigger thresholds = 42 47 97 127 168 223\n") # CMS AN-2010/371 not updated for 2011
    else:
        fcfg.write("Di-Jet trigger names = HLT_DiJetAve15U;HLT_DiJetAve30U;HLT_DiJetAve50U;HLT_DiJetAve70U;HLT_DiJetAve100U;HLT_DiJetAve140U\n")
        if(jettype == "ak5Calo"):
            fcfg.write("Di-Jet trigger thresholds = 38 59 86 111 147 196 \n")
        if(jettype == "ak5PF"):
            fcfg.write("Di-Jet trigger thresholds = 43 70 100 127 168 214 \n")
        if(jettype == "ak5JPT"):
           fcfg.write("Di-Jet trigger thresholds = 43 66 96 124 165 220\n") # CMS AN-2010/371
        if(jettype == "ak7JPT"):
            fcfg.write("Di-Jet trigger thresholds = 46 74 107 138 187 245\n") # CMS AN-2010/371
        if(jettype == "ak7PF"):
            fcfg.write("Di-Jet trigger thresholds = 47 79 111 140 187 240\n") # CMS AN-2010/371 
        if(jettype == "ak7Calo"):
            fcfg.write("Di-Jet trigger thresholds = 42 47 97 127 168 223\n") # CMS AN-2010/371 

        
    if(input != ""):
        fcfg.write("input calibration = Kalibri; "+input+"\n");

    fcfg.write("correct jets L1 = true\n");
    fcfg.write("fire all triggers = " + MC_fire_all_triggers +"\n");
    fcfg.close()
    return

#kostas binning:
#0: 0.0 1: 0.261 2: 0.522 3: 0.783 4: 0.957 5: 1.131 6: 1.305 7: 1.479
#8: 1.93 9: 2.322 10: 2.411 11: 2.5 12: 2.853 13: 2.964 14: 3.139 15: 3.489 5.191

#main program starts here!!!
#change these variables to steer the fit
DIR_JETALGO="AK5_entriescut_asymm_hardcut"
jetalgo="ak5"
PF_CALO_JPT="JPT"
CORRECTION="Su11_He_AK5"
#CORRECTION="Fall10_HenningJER_AK5"
#CORRECTION="ntuple"
#CORRECTION="ntuple"change other stuff...
#42x no corrs CORRECTIONS=CORRECTION+PF_CALO_JPT
CORRECTIONS=CORRECTION+PF_CALO_JPT
jettype = jetalgo+PF_CALO_JPT
jettype_import=jettype
DATAYEAR="2011"
DATATYPE="42X_combPrRe_ReRe"
MC="Su11"
MC_type="Z2wPU"
MC_fire_all_triggers="false"
#MinRunNumber=-1
MinRunNumber=163337
MaxRunNumber=1000000000
BINNING="k_HFfix"
#BINNING="kostas"


if(DATAYEAR == "2010"):
    if(DATATYPE=="ReReco"):
        datadir = "/scratch/hh/current/cms/user/stadie/DiJetNov4ReReco_v1C/merged"
#    if(DATATYPE=="Skim"):
#        datadir = "/scratch/hh/current/cms/user/stadie/DiJetNov4Skim_v1B"
    if(DATATYPE=="Skim"):
        datadir = "/scratch/hh/current/cms/user/stadie/2010/DiJetNov4Skim_v1C/merged"
if(DATAYEAR == "2011"):
    if(DATATYPE=="PrReco"):
        datadir = "/scratch/hh/current/cms/user/stadie/JetRun2011APromptRecoD/merged"
    if(DATATYPE=="PrRe42pb"):
        datadir = "/scratch/hh/current/cms/user/stadie/2011/JetRun2011APromptReco/Cert_160404-163369/merged"
    if(DATATYPE=="PrRe62pb"):
        datadir = "/afs/naf.desy.de/user/m/mschrode/lustre/data/Jet_Run2011A-PromptReco-v2_163337-163757"
    if(DATATYPE=="42X_corr"):
        datadir = "/scratch/hh/current/cms/user/stadie/2011/Jet2011AMay10ReReco_Cert_160404-163869/merged"
    if(DATATYPE=="42X_uncorr"):
        datadir = "/scratch/hh/current/cms/user/stadie/2011/Jet2011AMay10ReReco_Cert_160404-163869/merged"
    if(DATATYPE=="42X_PrRe"):
        datadir = "/scratch/hh/current/cms/user/stadie/2011/Jet2011APromptRecoV4_Cert_160404-165970/merged"
    if(DATATYPE=="42X_combPrRe_ReRe"):
#        datadir = "/afs/naf.desy.de/user/k/kirschen/scratch/2011_06_L2L3_Residuals_42X/combination_42XRereco_p_PrReco"
        datadir = "/afs/naf.desy.de/user/k/kirschen/scratch/2011_06_L2L3_Residuals_42X/combine_May10ReReco_and_166861"

##/scratch/hh/current/cms/user/stadie/JetRun2011APromptRecoD/merged



if(MC == "Su11"):
    if(MC_type=="D6TwPU"):
        datadirmc = "/scratch/hh/current/cms/user/stadie/2011/QCD_Pt-15to3000_TuneD6T_Flat_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/v4"
    if(MC_type=="Z2wPU"):
        datadirmc = "/scratch/hh/current/cms/user/stadie/2011/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2/v4"

if(MC == "S11"):
    if(MC_type=="wPU"):
        datadirmc = "/scratch/hh/current/cms/user/stadie/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1Amerged"
    if(MC_type=="wPU_SmConst"):
        datadirmc = "/scratch/hh/current/cms/user/kirschen/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1Amerged_smeared3"
    if(MC_type=="wPU_n"):
        datadirmc = "/scratch/hh/current/cms/user/stadie/2011/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1/B/merged"
    if(MC_type=="wPU_v2"):
        datadirmc = "/scratch/hh/current/cms/user/stadie/2011/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Spring11-PU_S2_START311_V2-v2_AODSIM/A/merged"

if(MC == "F10"):
    if(MC_type=="Z2wPU"):
        datadirmc="/scratch/hh/current/cms/user/stadie/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1Amerged" 
    if(MC_type=="Z2"):
        datadirmc = "/scratch/hh/current/cms/user/stadie/2010/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Fall10-START38_V12-v1E/merged"
    if(MC_type=="Z2_SmConst"):
        datadirmc = "/scratch/hh/current/cms/user/stadie/2010/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Fall10-START38_V12-v1E/smeared3"
#    if(MC_type=="Z2wPU"):





nthreads = 2
niothreads = 2
nevents =  -1
MAIN_dirname = "/afs/naf.desy.de/user/k/kirschen/scratch/2011_06_L2L3_Residuals_42X/"+DATAYEAR+DATATYPE+"_CORR" + CORRECTION +"_MC_"+MC+MC_type+"_kostas_"+ DIR_JETALGO
dirname = MAIN_dirname + "/dijetsFall10_TuneZ2_AK5"+PF_CALO_JPT+"_weighted_residuals_"+BINNING
useconstraint = False
batch = False
doBinnedFit = False
doUnbinnedFit = True

if(DATATYPE=="42X_corr" or DATATYPE=="42X_uncorr"  or DATATYPE=="42X_PrRe" or DATATYPE=="42X_combPrRe_ReRe"):
    if(PF_CALO_JPT=="PF"):
        jettype_import=jetalgo+"Fast"+PF_CALO_JPT


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


