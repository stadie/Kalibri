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
Eta max cut on jet         = 5.2  
Relative Rest Jet Cut      = 0.2      #NonLeadingJetsEt / PhotonEt
Min had fraction           = -0.05    #Default: 0.07
Max had fraction           = 1.05    #Default: 0.95
Et genJet min              = 0.0
Et genJet max              = 7000.0
DeltaR cut on jet matching = 9999
Max cut on relative n+1 Jet Et = 200.0 # was 0.4
Max cut on relative Soft Jet Et = 200.0

#---------------------------------------------------------------------------------
#   Input / Output
#---------------------------------------------------------------------------------

#jet correction source = JetMETCor
#jet correction name   = Spring10_AK5TRK
#jet correction source = JetMETCor
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

TwoJetsPtBalanceEvent plots names =   AsymmetryVsPt20_all_eta;AsymmetryVsEta;AsymmetryVsPt40;AsymmetryVsPt35;AsymmetryVsPt30;AsymmetryVsPt25;AsymmetryVsPt20;AsymmetryVsPt15;AsymmetryVsPt10;AsymmetryVsPt05;AsymmetryVsTJF;AbsAsymmetryVsPt40;AbsAsymmetryVsPt35;AbsAsymmetryVsPt30;AbsAsymmetryVsPt25;AbsAsymmetryVsPt20;AbsAsymmetryVsPt15;AbsAsymmetryVsPt10;AbsAsymmetryVsPt05;AbsAsymmetryVsTJF;OneBinAsymmetryVsPt40;OneBinAsymmetryVsPt35;OneBinAsymmetryVsPt30;OneBinAsymmetryVsPt25;OneBinAsymmetryVsPt20;OneBinAsymmetryVsPt15;OneBinAsymmetryVsPt10;OneBinAsymmetryVsPt05;OneBinAbsAsymmetryVsPt40;OneBinAbsAsymmetryVsPt35;OneBinAbsAsymmetryVsPt30;OneBinAbsAsymmetryVsPt25;OneBinAbsAsymmetryVsPt20;OneBinAbsAsymmetryVsPt15;OneBinAbsAsymmetryVsPt10;OneBinAbsAsymmetryVsPt05

AsymmetryVsPt20_all_eta x variable        =  MeanPt; log
AsymmetryVsPt20_all_eta x edges           =  30 20 2000
AsymmetryVsPt20_all_eta y variable        =  Asymmetry
AsymmetryVsPt20_all_eta y edges           =  51 -0.85 0.85 -0.5 0.5
AsymmetryVsPt20_all_eta cut variable      =  ThirdJetFractionPlain
AsymmetryVsPt20_all_eta cut edges         =  0.0 0.2
AsymmetryVsPt20_all_eta correction types  =  L2L3; L2L3Res
AsymmetryVsPt20_all_eta 1 correction types  =  L2L3
AsymmetryVsPt20_all_eta profile types     =  Mean
AsymmetryVsPt20_all_eta input samples     =  0:data;1:MC

AsymmetryVsEta x variable        =  Eta
AsymmetryVsEta x edges           =  26 -5.2 5.2
AsymmetryVsEta y variable        =  Asymmetry
AsymmetryVsEta y edges           =  51 -0.85 0.85 -0.5 0.8 -0.5 0.8 0.7 1.3 0.7 1.3
AsymmetryVsEta bin variable      =  MeanPt
AsymmetryVsEta bin edges         =  20 30 50 80 120 200 360 500 900 7000
AsymmetryVsEta cut variable      =  ThirdJetFractionPlain
AsymmetryVsEta cut edges         =  0.0 0.20
AsymmetryVsEta correction types  =  L2L3; L2L3Res
AsymmetryVsEta 1 correction types  =  L2L3
AsymmetryVsEta profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AsymmetryVsEta legend label      =  L2L3:CMS default
AsymmetryVsEta input samples     =  0:data;1:MC

AsymmetryVsPt40 x variable        =  MeanPt; log
AsymmetryVsPt40 x edges           =  30 20 2000
AsymmetryVsPt40 y variable        =  Asymmetry
AsymmetryVsPt40 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
AsymmetryVsPt40 bin variable      =  Eta
#AsymmetryVsPt40 bin edges        =  0.0 1.3 2.6 3.0 5.2
AsymmetryVsPt40 bin edges         =  -5.191 -3.489 -3.139 -2.964 -2.853 -2.5 -2.411 -2.322 -1.93 -1.479 -1.305 -1.131 -0.957 -0.783 -0.522 -0.261 0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191 
#AsymmetryVsPt40 bin edges         =  -6.0 -4.0 -3.5 -3.2 -3.0 -2.8 -2.6 -2.4 -2.2 -2.0 -1.8 -1.5 -1.4 -1.3 -1.2 -1.1 -0.9 -0.6 -0.3 0.0 0.3 0.6 0.9 1.1 1.2 1.3 1.4 1.5 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.2 3.5 4.0 6.0 
#AsymmetryVsPt40 bin edges         = -5.191 -3.489 -3.139 -2.964 -2.853 -2.5 -2.411 -2.322 -1.93 -1.479 -1.305 -1.131 -0.957 -0.783 -0.522 -0.261 0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191 
AsymmetryVsPt40 cut variable      =  ThirdJetFractionPlain
AsymmetryVsPt40 cut edges         =  0.0 0.40
AsymmetryVsPt40 correction types  =  L2L3; L2L3Res
AsymmetryVsPt40 1 correction types  =  L2L3
AsymmetryVsPt40 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AsymmetryVsPt40 legend label      =  L2L3:CMS default
AsymmetryVsPt40 input samples     =  0:data;1:MC

AsymmetryVsPt35 x variable        =  MeanPt; log
AsymmetryVsPt35 x edges           =  30 20 2000
AsymmetryVsPt35 y variable        =  Asymmetry
AsymmetryVsPt35 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
AsymmetryVsPt35 bin variable      =  Eta
#AsymmetryVsPt35 bin edges        =  0.0 1.3 2.6 3.0 5.2
AsymmetryVsPt35 bin edges         =  -5.191 -3.489 -3.139 -2.964 -2.853 -2.5 -2.411 -2.322 -1.93 -1.479 -1.305 -1.131 -0.957 -0.783 -0.522 -0.261 0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191 
AsymmetryVsPt35 cut variable      =  ThirdJetFractionPlain
AsymmetryVsPt35 cut edges         =  0.0 0.35
AsymmetryVsPt35 correction types  =  L2L3; L2L3Res
AsymmetryVsPt35 1 correction types  =  L2L3
AsymmetryVsPt35 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AsymmetryVsPt35 legend label      =  L2L3:CMS default
AsymmetryVsPt35 input samples     =  0:data;1:MC

AsymmetryVsPt30 x variable        =  MeanPt; log
AsymmetryVsPt30 x edges           =  30 20 2000
AsymmetryVsPt30 y variable        =  Asymmetry
AsymmetryVsPt30 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
AsymmetryVsPt30 bin variable      =  Eta
#AsymmetryVsPt30 bin edges        =  0.0 1.3 2.6 3.0 5.2
AsymmetryVsPt30 bin edges         =  -5.191 -3.489 -3.139 -2.964 -2.853 -2.5 -2.411 -2.322 -1.93 -1.479 -1.305 -1.131 -0.957 -0.783 -0.522 -0.261 0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191 
AsymmetryVsPt30 cut variable      =  ThirdJetFractionPlain
AsymmetryVsPt30 cut edges         =  0.0 0.30
AsymmetryVsPt30 correction types  =  L2L3; L2L3Res
AsymmetryVsPt30 1 correction types  =  L2L3
AsymmetryVsPt30 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AsymmetryVsPt30 legend label      =  L2L3:CMS default
AsymmetryVsPt30 input samples     =  0:data;1:MC

AsymmetryVsPt25 x variable        =  MeanPt; log
AsymmetryVsPt25 x edges           =  30 20 2000
AsymmetryVsPt25 y variable        =  Asymmetry
AsymmetryVsPt25 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
AsymmetryVsPt25 bin variable      =  Eta
#AsymmetryVsPt25 bin edges        =  0.0 1.3 2.6 3.0 5.2
AsymmetryVsPt25 bin edges         =  -5.191 -3.489 -3.139 -2.964 -2.853 -2.5 -2.411 -2.322 -1.93 -1.479 -1.305 -1.131 -0.957 -0.783 -0.522 -0.261 0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191 
AsymmetryVsPt25 cut variable      =  ThirdJetFractionPlain
AsymmetryVsPt25 cut edges         =  0.0 0.25
AsymmetryVsPt25 correction types  =  L2L3; L2L3Res
AsymmetryVsPt25 1 correction types  =  L2L3
AsymmetryVsPt25 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AsymmetryVsPt25 legend label      =  L2L3:CMS default
AsymmetryVsPt25 input samples     =  0:data;1:MC

AsymmetryVsPt20 x variable        =  MeanPt; log
AsymmetryVsPt20 x edges           =  30 20 2000
AsymmetryVsPt20 y variable        =  Asymmetry
AsymmetryVsPt20 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
AsymmetryVsPt20 bin variable      =  Eta
#AsymmetryVsPt20 bin edges        =  0.0 1.3 2.6 3.0 5.2
AsymmetryVsPt20 bin edges         =  -5.191 -3.489 -3.139 -2.964 -2.853 -2.5 -2.411 -2.322 -1.93 -1.479 -1.305 -1.131 -0.957 -0.783 -0.522 -0.261 0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191 
AsymmetryVsPt20 cut variable      =  ThirdJetFractionPlain
AsymmetryVsPt20 cut edges         =  0.0 0.2
AsymmetryVsPt20 correction types  =  L2L3; L2L3Res
AsymmetryVsPt20 1 correction types  =  L2L3
AsymmetryVsPt20 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AsymmetryVsPt20 legend label      =  L2L3:CMS default
AsymmetryVsPt20 input samples     =  0:data;1:MC

AsymmetryVsPt15 x variable        =  MeanPt; log
AsymmetryVsPt15 x edges           =  30 20 2000
AsymmetryVsPt15 y variable        =  Asymmetry
AsymmetryVsPt15 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
AsymmetryVsPt15 bin variable      =  Eta
#AsymmetryVsPt15 bin edges        =  0.0 1.3 2.6 3.0 5.2
AsymmetryVsPt15 bin edges         =  -5.191 -3.489 -3.139 -2.964 -2.853 -2.5 -2.411 -2.322 -1.93 -1.479 -1.305 -1.131 -0.957 -0.783 -0.522 -0.261 0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191 
AsymmetryVsPt15 cut variable      =  ThirdJetFractionPlain
AsymmetryVsPt15 cut edges         =  0.0 0.15
AsymmetryVsPt15 correction types  =  L2L3; L2L3Res
AsymmetryVsPt15 1 correction types  =  L2L3
AsymmetryVsPt15 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AsymmetryVsPt15 legend label      =  L2L3:CMS default
AsymmetryVsPt15 input samples     =  0:data;1:MC

AsymmetryVsPt10 x variable        =  MeanPt; log
AsymmetryVsPt10 x edges           =  30 20 2000
AsymmetryVsPt10 y variable        =  Asymmetry
AsymmetryVsPt10 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
AsymmetryVsPt10 bin variable      =  Eta
#AsymmetryVsPt10 bin edges        =  0.0 1.3 2.6 3.0 5.2
AsymmetryVsPt10 bin edges         =  -5.191 -3.489 -3.139 -2.964 -2.853 -2.5 -2.411 -2.322 -1.93 -1.479 -1.305 -1.131 -0.957 -0.783 -0.522 -0.261 0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191 
AsymmetryVsPt10 cut variable      =  ThirdJetFractionPlain
AsymmetryVsPt10 cut edges         =  0.0 0.1
AsymmetryVsPt10 correction types  =  L2L3; L2L3Res
AsymmetryVsPt10 1 correction types  =  L2L3
AsymmetryVsPt10 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AsymmetryVsPt10 legend label      =  L2L3:CMS default
AsymmetryVsPt10 input samples     =  0:data;1:MC

AsymmetryVsPt05 x variable        =  MeanPt; log
AsymmetryVsPt05 x edges           =  30 20 2000
AsymmetryVsPt05 y variable        =  Asymmetry
AsymmetryVsPt05 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
AsymmetryVsPt05 bin variable      =  Eta
#AsymmetryVsPt05 bin edges        =  0.0 1.3 2.6 3.0 5.2
AsymmetryVsPt05 bin edges         =  -5.191 -3.489 -3.139 -2.964 -2.853 -2.5 -2.411 -2.322 -1.93 -1.479 -1.305 -1.131 -0.957 -0.783 -0.522 -0.261 0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191 
AsymmetryVsPt05 cut variable      =  ThirdJetFractionPlain
AsymmetryVsPt05 cut edges         =  0.0 0.05
AsymmetryVsPt05 correction types  =  L2L3; L2L3Res
AsymmetryVsPt05 1 correction types  =  L2L3
AsymmetryVsPt05 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AsymmetryVsPt05 legend label      =  L2L3:CMS default
AsymmetryVsPt05 input samples     =  0:data;1:MC

AsymmetryVsTJF x variable        =  ThirdJetFractionPlain
AsymmetryVsTJF x edges           =  8 0.00 0.40
AsymmetryVsTJF y variable        =  Asymmetry
AsymmetryVsTJF y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
AsymmetryVsTJF bin variable      =  Eta
AsymmetryVsTJF bin edges         =  -5.191 -3.489 -3.139 -2.964 -2.853 -2.5 -2.411 -2.322 -1.93 -1.479 -1.305 -1.131 -0.957 -0.783 -0.522 -0.261 0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191 
AsymmetryVsTJF correction types  =  L2L3; L2L3Res
AsymmetryVsTJF 1 correction types  =  L2L3
AsymmetryVsTJF profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AsymmetryVTJF legend label      =  L2L3:CMS default
AsymmetryVsTJF input samples     =  0:data;1:MC


AbsAsymmetryVsPt40 x variable        =  MeanPt; log
AbsAsymmetryVsPt40 x edges           =  30 20 2000
AbsAsymmetryVsPt40 y variable        =  Asymmetry
AbsAsymmetryVsPt40 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
AbsAsymmetryVsPt40 bin variable      =  AbsEta
#AbsAsymmetryVsPt40 bin edges        =  0.0 1.3 2.6 3.0 5.2
AbsAsymmetryVsPt40 bin edges         =  0.0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191
#AbsAsymmetryVsPt40 bin edges         =  0.0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191
AbsAsymmetryVsPt40 cut variable      =  ThirdJetFractionPlain
AbsAsymmetryVsPt40 cut edges         =  0.0 0.40
AbsAsymmetryVsPt40 correction types  =  L2L3; L2L3Res
AbsAsymmetryVsPt40 1 correction types  =  L2L3
AbsAsymmetryVsPt40 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AbsAsymmetryVsPt40 legend label      =  L2L3:CMS default
AbsAsymmetryVsPt40 input samples     =  0:data;1:MC

AbsAsymmetryVsPt35 x variable        =  MeanPt; log
AbsAsymmetryVsPt35 x edges           =  30 20 2000
AbsAsymmetryVsPt35 y variable        =  Asymmetry
AbsAsymmetryVsPt35 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
AbsAsymmetryVsPt35 bin variable      =  AbsEta
#AbsAsymmetryVsPt35 bin edges        =  0.0 1.3 2.6 3.0 5.2
AbsAsymmetryVsPt35 bin edges         =  0.0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191
AbsAsymmetryVsPt35 cut variable      =  ThirdJetFractionPlain
AbsAsymmetryVsPt35 cut edges         =  0.0 0.35
AbsAsymmetryVsPt35 correction types  =  L2L3; L2L3Res
AbsAsymmetryVsPt35 1 correction types  =  L2L3
AbsAsymmetryVsPt35 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AbsAsymmetryVsPt35 legend label      =  L2L3:CMS default
AbsAsymmetryVsPt35 input samples     =  0:data;1:MC

AbsAsymmetryVsPt30 x variable        =  MeanPt; log
AbsAsymmetryVsPt30 x edges           =  30 20 2000
AbsAsymmetryVsPt30 y variable        =  Asymmetry
AbsAsymmetryVsPt30 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
AbsAsymmetryVsPt30 bin variable      =  AbsEta
#AbsAsymmetryVsPt30 bin edges        =  0.0 1.3 2.6 3.0 5.2
AbsAsymmetryVsPt30 bin edges         =  0.0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191
AbsAsymmetryVsPt30 cut variable      =  ThirdJetFractionPlain
AbsAsymmetryVsPt30 cut edges         =  0.0 0.30
AbsAsymmetryVsPt30 correction types  =  L2L3; L2L3Res
AbsAsymmetryVsPt30 1 correction types  =  L2L3
AbsAsymmetryVsPt30 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AbsAsymmetryVsPt30 legend label      =  L2L3:CMS default
AbsAsymmetryVsPt30 input samples     =  0:data;1:MC

AbsAsymmetryVsPt25 x variable        =  MeanPt; log
AbsAsymmetryVsPt25 x edges           =  30 20 2000
AbsAsymmetryVsPt25 y variable        =  Asymmetry
AbsAsymmetryVsPt25 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
AbsAsymmetryVsPt25 bin variable      =  AbsEta
#AbsAsymmetryVsPt25 bin edges        =  0.0 1.3 2.6 3.0 5.2
AbsAsymmetryVsPt25 bin edges         =  0.0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191
AbsAsymmetryVsPt25 cut variable      =  ThirdJetFractionPlain
AbsAsymmetryVsPt25 cut edges         =  0.0 0.25
AbsAsymmetryVsPt25 correction types  =  L2L3; L2L3Res
AbsAsymmetryVsPt25 1 correction types  =  L2L3
AbsAsymmetryVsPt25 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AbsAsymmetryVsPt25 legend label      =  L2L3:CMS default
AbsAsymmetryVsPt25 input samples     =  0:data;1:MC

AbsAsymmetryVsPt20 x variable        =  MeanPt; log
AbsAsymmetryVsPt20 x edges           =  30 20 2000
AbsAsymmetryVsPt20 y variable        =  Asymmetry
AbsAsymmetryVsPt20 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
AbsAsymmetryVsPt20 bin variable      =  AbsEta
#AbsAsymmetryVsPt20 bin edges        =  0.0 1.3 2.6 3.0 5.2
AbsAsymmetryVsPt20 bin edges         =  0.0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191
AbsAsymmetryVsPt20 cut variable      =  ThirdJetFractionPlain
AbsAsymmetryVsPt20 cut edges         =  0.0 0.2
AbsAsymmetryVsPt20 correction types  =  L2L3; L2L3Res
AbsAsymmetryVsPt20 1 correction types  =  L2L3
AbsAsymmetryVsPt20 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AbsAsymmetryVsPt20 legend label      =  L2L3:CMS default
AbsAsymmetryVsPt20 input samples     =  0:data;1:MC

AbsAsymmetryVsPt15 x variable        =  MeanPt; log
AbsAsymmetryVsPt15 x edges           =  30 20 2000
AbsAsymmetryVsPt15 y variable        =  Asymmetry
AbsAsymmetryVsPt15 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
AbsAsymmetryVsPt15 bin variable      =  AbsEta
#AbsAsymmetryVsPt15 bin edges        =  0.0 1.3 2.6 3.0 5.2
AbsAsymmetryVsPt15 bin edges         =  0.0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191
AbsAsymmetryVsPt15 cut variable      =  ThirdJetFractionPlain
AbsAsymmetryVsPt15 cut edges         =  0.0 0.15
AbsAsymmetryVsPt15 correction types  =  L2L3; L2L3Res
AbsAsymmetryVsPt15 1 correction types  =  L2L3
AbsAsymmetryVsPt15 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AbsAsymmetryVsPt15 legend label      =  L2L3:CMS default
AbsAsymmetryVsPt15 input samples     =  0:data;1:MC

AbsAsymmetryVsPt10 x variable        =  MeanPt; log
AbsAsymmetryVsPt10 x edges           =  30 20 2000
AbsAsymmetryVsPt10 y variable        =  Asymmetry
AbsAsymmetryVsPt10 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
AbsAsymmetryVsPt10 bin variable      =  AbsEta
#AbsAsymmetryVsPt10 bin edges        =  0.0 1.3 2.6 3.0 5.2
AbsAsymmetryVsPt10 bin edges         =  0.0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191
AbsAsymmetryVsPt10 cut variable      =  ThirdJetFractionPlain
AbsAsymmetryVsPt10 cut edges         =  0.0 0.1
AbsAsymmetryVsPt10 correction types  =  L2L3; L2L3Res
AbsAsymmetryVsPt10 1 correction types  =  L2L3
AbsAsymmetryVsPt10 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AbsAsymmetryVsPt10 legend label      =  L2L3:CMS default
AbsAsymmetryVsPt10 input samples     =  0:data;1:MC

AbsAsymmetryVsPt05 x variable        =  MeanPt; log
AbsAsymmetryVsPt05 x edges           =  30 20 2000
AbsAsymmetryVsPt05 y variable        =  Asymmetry
AbsAsymmetryVsPt05 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
AbsAsymmetryVsPt05 bin variable      =  AbsEta
#AbsAsymmetryVsPt05 bin edges        =  0.0 1.3 2.6 3.0 5.2
AbsAsymmetryVsPt05 bin edges         =  0.0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191
AbsAsymmetryVsPt05 cut variable      =  ThirdJetFractionPlain
AbsAsymmetryVsPt05 cut edges         =  0.0 0.05
AbsAsymmetryVsPt05 correction types  =  L2L3; L2L3Res
AbsAsymmetryVsPt05 1 correction types  =  L2L3
AbsAsymmetryVsPt05 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AbsAsymmetryVsPt05 legend label      =  L2L3:CMS default
AbsAsymmetryVsPt05 input samples     =  0:data;1:MC

AbsAsymmetryVsTJF x variable        =  ThirdJetFractionPlain
AbsAsymmetryVsTJF x edges           =  8 0.00 0.40
AbsAsymmetryVsTJF y variable        =  Asymmetry
AbsAsymmetryVsTJF y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
AbsAsymmetryVsTJF bin variable      =  AbsEta
AbsAsymmetryVsTJF bin edges         =  0.0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191
AbsAsymmetryVsTJF correction types  =  L2L3; L2L3Res
AbsAsymmetryVsTJF 1 correction types  =  L2L3
AbsAsymmetryVsTJF profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#AsymmetryVTJF legend label      =  L2L3:CMS default
AbsAsymmetryVsTJF input samples     =  0:data;1:MC







OneBinAsymmetryVsPt40 x variable        =  MeanPt; log
OneBinAsymmetryVsPt40 x edges           =  1 20 2000
OneBinAsymmetryVsPt40 y variable        =  Asymmetry
OneBinAsymmetryVsPt40 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
OneBinAsymmetryVsPt40 bin variable      =  Eta
#OneBinAsymmetryVsPt40 bin edges        =  0.0 1.3 2.6 3.0 5.2
OneBinAsymmetryVsPt40 bin edges         =  -5.191 -3.489 -3.139 -2.964 -2.853 -2.5 -2.411 -2.322 -1.93 -1.479 -1.305 -1.131 -0.957 -0.783 -0.522 -0.261 0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191  
#OneBinAsymmetryVsPt40 bin edges         = -6.0 -4.0 -3.5 -2.7 -2.4 -2.1 -1.8 -1.5 -1.3 -1.1 -0.9 -0.6 -0.3 0.0 0.3 0.6 0.9 1.1 1.3 1.5 1.8 2.1 2.4 2.7 3.0 3.5 4.0 6.0 
OneBinAsymmetryVsPt40 cut variable      =  ThirdJetFractionPlain
OneBinAsymmetryVsPt40 cut edges         =  0.0 0.40
OneBinAsymmetryVsPt40 correction types  =  L2L3; L2L3Res
OneBinAsymmetryVsPt40 1 correction types  =  L2L3
OneBinAsymmetryVsPt40 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#OneBinAsymmetryVsPt40 legend label      =  L2L3:CMS default
OneBinAsymmetryVsPt40 input samples     =  0:data;1:MC

OneBinAsymmetryVsPt35 x variable        =  MeanPt; log
OneBinAsymmetryVsPt35 x edges           =  1 20 2000
OneBinAsymmetryVsPt35 y variable        =  Asymmetry
OneBinAsymmetryVsPt35 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
OneBinAsymmetryVsPt35 bin variable      =  Eta
#OneBinAsymmetryVsPt35 bin edges        =  0.0 1.3 2.6 3.0 5.2
OneBinAsymmetryVsPt35 bin edges         =  -5.191 -3.489 -3.139 -2.964 -2.853 -2.5 -2.411 -2.322 -1.93 -1.479 -1.305 -1.131 -0.957 -0.783 -0.522 -0.261 0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191  
OneBinAsymmetryVsPt35 cut variable      =  ThirdJetFractionPlain
OneBinAsymmetryVsPt35 cut edges         =  0.0 0.35
OneBinAsymmetryVsPt35 correction types  =  L2L3; L2L3Res
OneBinAsymmetryVsPt35 1 correction types  =  L2L3
OneBinAsymmetryVsPt35 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#OneBinAsymmetryVsPt35 legend label      =  L2L3:CMS default
OneBinAsymmetryVsPt35 input samples     =  0:data;1:MC

OneBinAsymmetryVsPt30 x variable        =  MeanPt; log
OneBinAsymmetryVsPt30 x edges           =  1 20 2000
OneBinAsymmetryVsPt30 y variable        =  Asymmetry
OneBinAsymmetryVsPt30 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
OneBinAsymmetryVsPt30 bin variable      =  Eta
#OneBinAsymmetryVsPt30 bin edges        =  0.0 1.3 2.6 3.0 5.2
OneBinAsymmetryVsPt30 bin edges         =  -5.191 -3.489 -3.139 -2.964 -2.853 -2.5 -2.411 -2.322 -1.93 -1.479 -1.305 -1.131 -0.957 -0.783 -0.522 -0.261 0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191  
OneBinAsymmetryVsPt30 cut variable      =  ThirdJetFractionPlain
OneBinAsymmetryVsPt30 cut edges         =  0.0 0.30
OneBinAsymmetryVsPt30 correction types  =  L2L3; L2L3Res
OneBinAsymmetryVsPt30 1 correction types  =  L2L3
OneBinAsymmetryVsPt30 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#OneBinAsymmetryVsPt30 legend label      =  L2L3:CMS default
OneBinAsymmetryVsPt30 input samples     =  0:data;1:MC

OneBinAsymmetryVsPt25 x variable        =  MeanPt; log
OneBinAsymmetryVsPt25 x edges           =  1 20 2000
OneBinAsymmetryVsPt25 y variable        =  Asymmetry
OneBinAsymmetryVsPt25 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
OneBinAsymmetryVsPt25 bin variable      =  Eta
#OneBinAsymmetryVsPt25 bin edges        =  0.0 1.3 2.6 3.0 5.2
OneBinAsymmetryVsPt25 bin edges         =  -5.191 -3.489 -3.139 -2.964 -2.853 -2.5 -2.411 -2.322 -1.93 -1.479 -1.305 -1.131 -0.957 -0.783 -0.522 -0.261 0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191  
OneBinAsymmetryVsPt25 cut variable      =  ThirdJetFractionPlain
OneBinAsymmetryVsPt25 cut edges         =  0.0 0.25
OneBinAsymmetryVsPt25 correction types  =  L2L3; L2L3Res
OneBinAsymmetryVsPt25 1 correction types  =  L2L3
OneBinAsymmetryVsPt25 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#OneBinAsymmetryVsPt25 legend label      =  L2L3:CMS default
OneBinAsymmetryVsPt25 input samples     =  0:data;1:MC

OneBinAsymmetryVsPt20 x variable        =  MeanPt; log
OneBinAsymmetryVsPt20 x edges           =  1 20 2000
OneBinAsymmetryVsPt20 y variable        =  Asymmetry
OneBinAsymmetryVsPt20 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
OneBinAsymmetryVsPt20 bin variable      =  Eta
#OneBinAsymmetryVsPt20 bin edges        =  0.0 1.3 2.6 3.0 5.2
OneBinAsymmetryVsPt20 bin edges         =  -5.191 -3.489 -3.139 -2.964 -2.853 -2.5 -2.411 -2.322 -1.93 -1.479 -1.305 -1.131 -0.957 -0.783 -0.522 -0.261 0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191  
OneBinAsymmetryVsPt20 cut variable      =  ThirdJetFractionPlain
OneBinAsymmetryVsPt20 cut edges         =  0.0 0.2
OneBinAsymmetryVsPt20 correction types  =  L2L3; L2L3Res
OneBinAsymmetryVsPt20 1 correction types  =  L2L3
OneBinAsymmetryVsPt20 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#OneBinAsymmetryVsPt20 legend label      =  L2L3:CMS default
OneBinAsymmetryVsPt20 input samples     =  0:data;1:MC

OneBinAsymmetryVsPt15 x variable        =  MeanPt; log
OneBinAsymmetryVsPt15 x edges           =  1 20 2000
OneBinAsymmetryVsPt15 y variable        =  Asymmetry
OneBinAsymmetryVsPt15 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
OneBinAsymmetryVsPt15 bin variable      =  Eta
#OneBinAsymmetryVsPt15 bin edges        =  0.0 1.3 2.6 3.0 5.2
OneBinAsymmetryVsPt15 bin edges         =  -5.191 -3.489 -3.139 -2.964 -2.853 -2.5 -2.411 -2.322 -1.93 -1.479 -1.305 -1.131 -0.957 -0.783 -0.522 -0.261 0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191  
OneBinAsymmetryVsPt15 cut variable      =  ThirdJetFractionPlain
OneBinAsymmetryVsPt15 cut edges         =  0.0 0.15
OneBinAsymmetryVsPt15 correction types  =  L2L3; L2L3Res
OneBinAsymmetryVsPt15 1 correction types  =  L2L3
OneBinAsymmetryVsPt15 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#OneBinAsymmetryVsPt15 legend label      =  L2L3:CMS default
OneBinAsymmetryVsPt15 input samples     =  0:data;1:MC

OneBinAsymmetryVsPt10 x variable        =  MeanPt; log
OneBinAsymmetryVsPt10 x edges           =  1 20 2000
OneBinAsymmetryVsPt10 y variable        =  Asymmetry
OneBinAsymmetryVsPt10 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
OneBinAsymmetryVsPt10 bin variable      =  Eta
#OneBinAsymmetryVsPt10 bin edges        =  0.0 1.3 2.6 3.0 5.2
OneBinAsymmetryVsPt10 bin edges         =  -5.191 -3.489 -3.139 -2.964 -2.853 -2.5 -2.411 -2.322 -1.93 -1.479 -1.305 -1.131 -0.957 -0.783 -0.522 -0.261 0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191  
OneBinAsymmetryVsPt10 cut variable      =  ThirdJetFractionPlain
OneBinAsymmetryVsPt10 cut edges         =  0.0 0.1
OneBinAsymmetryVsPt10 correction types  =  L2L3; L2L3Res
OneBinAsymmetryVsPt10 1 correction types  =  L2L3
OneBinAsymmetryVsPt10 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#OneBinAsymmetryVsPt10 legend label      =  L2L3:CMS default
OneBinAsymmetryVsPt10 input samples     =  0:data;1:MC

OneBinAsymmetryVsPt05 x variable        =  MeanPt; log
OneBinAsymmetryVsPt05 x edges           =  1 20 2000
OneBinAsymmetryVsPt05 y variable        =  Asymmetry
OneBinAsymmetryVsPt05 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
OneBinAsymmetryVsPt05 bin variable      =  Eta
#OneBinAsymmetryVsPt05 bin edges        =  0.0 1.3 2.6 3.0 5.2
OneBinAsymmetryVsPt05 bin edges         =  -5.191 -3.489 -3.139 -2.964 -2.853 -2.5 -2.411 -2.322 -1.93 -1.479 -1.305 -1.131 -0.957 -0.783 -0.522 -0.261 0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191  
OneBinAsymmetryVsPt05 cut variable      =  ThirdJetFractionPlain
OneBinAsymmetryVsPt05 cut edges         =  0.0 0.05
OneBinAsymmetryVsPt05 correction types  =  L2L3; L2L3Res
OneBinAsymmetryVsPt05 1 correction types  =  L2L3
OneBinAsymmetryVsPt05 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#OneBinAsymmetryVsPt05 legend label      =  L2L3:CMS default
OneBinAsymmetryVsPt05 input samples     =  0:data;1:MC


OneBinAbsAsymmetryVsPt40 x variable        =  MeanPt; log
OneBinAbsAsymmetryVsPt40 x edges           =  1 20 2000
OneBinAbsAsymmetryVsPt40 y variable        =  Asymmetry
OneBinAbsAsymmetryVsPt40 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
OneBinAbsAsymmetryVsPt40 bin variable      =  AbsEta
#OneBinAbsAsymmetryVsPt40 bin edges        =  0.0 1.3 2.6 3.0 5.2
OneBinAbsAsymmetryVsPt40 bin edges         =  0.0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191
#OneBinAbsAsymmetryVsPt40 bin edges         =  0.0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191
OneBinAbsAsymmetryVsPt40 cut variable      =  ThirdJetFractionPlain
OneBinAbsAsymmetryVsPt40 cut edges         =  0.0 0.40
OneBinAbsAsymmetryVsPt40 correction types  =  L2L3; L2L3Res
OneBinAbsAsymmetryVsPt40 1 correction types  =  L2L3
OneBinAbsAsymmetryVsPt40 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#OneBinAbsAsymmetryVsPt40 legend label      =  L2L3:CMS default
OneBinAbsAsymmetryVsPt40 input samples     =  0:data;1:MC

OneBinAbsAsymmetryVsPt35 x variable        =  MeanPt; log
OneBinAbsAsymmetryVsPt35 x edges           =  1 20 2000
OneBinAbsAsymmetryVsPt35 y variable        =  Asymmetry
OneBinAbsAsymmetryVsPt35 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
OneBinAbsAsymmetryVsPt35 bin variable      =  AbsEta
#OneBinAbsAsymmetryVsPt35 bin edges        =  0.0 1.3 2.6 3.0 5.2
OneBinAbsAsymmetryVsPt35 bin edges         =  0.0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191
OneBinAbsAsymmetryVsPt35 cut variable      =  ThirdJetFractionPlain
OneBinAbsAsymmetryVsPt35 cut edges         =  0.0 0.35
OneBinAbsAsymmetryVsPt35 correction types  =  L2L3; L2L3Res
OneBinAbsAsymmetryVsPt35 1 correction types  =  L2L3
OneBinAbsAsymmetryVsPt35 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#OneBinAbsAsymmetryVsPt35 legend label      =  L2L3:CMS default
OneBinAbsAsymmetryVsPt35 input samples     =  0:data;1:MC

OneBinAbsAsymmetryVsPt30 x variable        =  MeanPt; log
OneBinAbsAsymmetryVsPt30 x edges           =  1 20 2000
OneBinAbsAsymmetryVsPt30 y variable        =  Asymmetry
OneBinAbsAsymmetryVsPt30 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
OneBinAbsAsymmetryVsPt30 bin variable      =  AbsEta
#OneBinAbsAsymmetryVsPt30 bin edges        =  0.0 1.3 2.6 3.0 5.2
OneBinAbsAsymmetryVsPt30 bin edges         =  0.0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191
OneBinAbsAsymmetryVsPt30 cut variable      =  ThirdJetFractionPlain
OneBinAbsAsymmetryVsPt30 cut edges         =  0.0 0.30
OneBinAbsAsymmetryVsPt30 correction types  =  L2L3; L2L3Res
OneBinAbsAsymmetryVsPt30 1 correction types  =  L2L3
OneBinAbsAsymmetryVsPt30 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#OneBinAbsAsymmetryVsPt30 legend label      =  L2L3:CMS default
OneBinAbsAsymmetryVsPt30 input samples     =  0:data;1:MC

OneBinAbsAsymmetryVsPt25 x variable        =  MeanPt; log
OneBinAbsAsymmetryVsPt25 x edges           =  1 20 2000
OneBinAbsAsymmetryVsPt25 y variable        =  Asymmetry
OneBinAbsAsymmetryVsPt25 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
OneBinAbsAsymmetryVsPt25 bin variable      =  AbsEta
#OneBinAbsAsymmetryVsPt25 bin edges        =  0.0 1.3 2.6 3.0 5.2
OneBinAbsAsymmetryVsPt25 bin edges         =  0.0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191
OneBinAbsAsymmetryVsPt25 cut variable      =  ThirdJetFractionPlain
OneBinAbsAsymmetryVsPt25 cut edges         =  0.0 0.25
OneBinAbsAsymmetryVsPt25 correction types  =  L2L3; L2L3Res
OneBinAbsAsymmetryVsPt25 1 correction types  =  L2L3
OneBinAbsAsymmetryVsPt25 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#OneBinAbsAsymmetryVsPt25 legend label      =  L2L3:CMS default
OneBinAbsAsymmetryVsPt25 input samples     =  0:data;1:MC

OneBinAbsAsymmetryVsPt20 x variable        =  MeanPt; log
OneBinAbsAsymmetryVsPt20 x edges           =  1 20 2000
OneBinAbsAsymmetryVsPt20 y variable        =  Asymmetry
OneBinAbsAsymmetryVsPt20 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
OneBinAbsAsymmetryVsPt20 bin variable      =  AbsEta
#OneBinAbsAsymmetryVsPt20 bin edges        =  0.0 1.3 2.6 3.0 5.2
OneBinAbsAsymmetryVsPt20 bin edges         =  0.0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191
OneBinAbsAsymmetryVsPt20 cut variable      =  ThirdJetFractionPlain
OneBinAbsAsymmetryVsPt20 cut edges         =  0.0 0.2
OneBinAbsAsymmetryVsPt20 correction types  =  L2L3; L2L3Res
OneBinAbsAsymmetryVsPt20 1 correction types  =  L2L3
OneBinAbsAsymmetryVsPt20 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#OneBinAbsAsymmetryVsPt20 legend label      =  L2L3:CMS default
OneBinAbsAsymmetryVsPt20 input samples     =  0:data;1:MC

OneBinAbsAsymmetryVsPt15 x variable        =  MeanPt; log
OneBinAbsAsymmetryVsPt15 x edges           =  1 20 2000
OneBinAbsAsymmetryVsPt15 y variable        =  Asymmetry
OneBinAbsAsymmetryVsPt15 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
OneBinAbsAsymmetryVsPt15 bin variable      =  AbsEta
#OneBinAbsAsymmetryVsPt15 bin edges        =  0.0 1.3 2.6 3.0 5.2
OneBinAbsAsymmetryVsPt15 bin edges         =  0.0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191
OneBinAbsAsymmetryVsPt15 cut variable      =  ThirdJetFractionPlain
OneBinAbsAsymmetryVsPt15 cut edges         =  0.0 0.15
OneBinAbsAsymmetryVsPt15 correction types  =  L2L3; L2L3Res
OneBinAbsAsymmetryVsPt15 1 correction types  =  L2L3
OneBinAbsAsymmetryVsPt15 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#OneBinAbsAsymmetryVsPt15 legend label      =  L2L3:CMS default
OneBinAbsAsymmetryVsPt15 input samples     =  0:data;1:MC

OneBinAbsAsymmetryVsPt10 x variable        =  MeanPt; log
OneBinAbsAsymmetryVsPt10 x edges           =  1 20 2000
OneBinAbsAsymmetryVsPt10 y variable        =  Asymmetry
OneBinAbsAsymmetryVsPt10 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
OneBinAbsAsymmetryVsPt10 bin variable      =  AbsEta
#OneBinAbsAsymmetryVsPt10 bin edges        =  0.0 1.3 2.6 3.0 5.2
OneBinAbsAsymmetryVsPt10 bin edges         =  0.0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191
OneBinAbsAsymmetryVsPt10 cut variable      =  ThirdJetFractionPlain
OneBinAbsAsymmetryVsPt10 cut edges         =  0.0 0.1
OneBinAbsAsymmetryVsPt10 correction types  =  L2L3; L2L3Res
OneBinAbsAsymmetryVsPt10 1 correction types  =  L2L3
OneBinAbsAsymmetryVsPt10 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#OneBinAbsAsymmetryVsPt10 legend label      =  L2L3:CMS default
OneBinAbsAsymmetryVsPt10 input samples     =  0:data;1:MC

OneBinAbsAsymmetryVsPt05 x variable        =  MeanPt; log
OneBinAbsAsymmetryVsPt05 x edges           =  1 20 2000
OneBinAbsAsymmetryVsPt05 y variable        =  Asymmetry
OneBinAbsAsymmetryVsPt05 y edges           =  51 -0.85 0.85 -0.5 0.5 -0.5 0.5 0.7 1.3 0.7 1.3
OneBinAbsAsymmetryVsPt05 bin variable      =  AbsEta
#OneBinAbsAsymmetryVsPt05 bin edges        =  0.0 1.3 2.6 3.0 5.2
OneBinAbsAsymmetryVsPt05 bin edges         =  0.0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191
OneBinAbsAsymmetryVsPt05 cut variable      =  ThirdJetFractionPlain
OneBinAbsAsymmetryVsPt05 cut edges         =  0.0 0.05
OneBinAbsAsymmetryVsPt05 correction types  =  L2L3; L2L3Res
OneBinAbsAsymmetryVsPt05 1 correction types  =  L2L3
OneBinAbsAsymmetryVsPt05 profile types     =  Mean; GaussFitMean; RatioOfMeans; RatioOfGaussFitMeans
#OneBinAbsAsymmetryVsPt05 legend label      =  L2L3:CMS default
OneBinAbsAsymmetryVsPt05 input samples     =  0:data;1:MC



"""
    fcfg = open(filename, "w")
    fcfg.write(config)
    fcfg.write("jet correction source = JetMETCor\n")
    fcfg.write("jet correction name   = "+CORRECTIONS+"\n")
    fcfg.write("Di-Jet input file = dijetlist\n")
    fcfg.write("Di-Jet Control1 input file = mcdijetlist\n")
    fcfg.write("Output file       = "+output+"\n");
    fcfg.write("Number of Threads = "+str(nthreads)+"\n")
#    fcfg.write("Number of IO Threads = 5\n")
    fcfg.write("Number of IO Threads = "+str(niothreads)+"\n")
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


#main program starts here!!!
#change these variables to steer the fit
DIR_JETALGO="/afs/naf.desy.de/user/k/kirschen/scratch/2011_02_Fall10_residuals_kostas_validation/NEW_after_JetMET_10_corr_10_MC_10skim_data_AK5"
CORRECTIONS="Fall10_AK5Calo"
jetalgo="ak5"
PF_CALO_JPT="Calo"
jettype = jetalgo+PF_CALO_JPT
MC_fire_all_triggers="false"
#datadir = "/scratch/hh/current/cms/user/stadie/QCDFlat_Pt15to3000Spring10-START3X_V26_S09-v1C"
#datadir = "/scratch/hh/current/cms/user/stadie/JetMET_Run2010A-PromptReco-v4_DCSONLY"
#datadir = "/scratch/hh/current/cms/user/stadie/Spring10QCDDiJetPt*"
#sollte C sein.. Hartmut sagen...
#datadir = "/scratch/hh/current/cms/user/stadie/DiJetNov4ReReco_v1C/merged"
#calib skim 2010
datadir = "/scratch/hh/current/cms/user/stadie/DiJetNov4Skim_v1A/merged"
#datadir = "/scratch/hh/current/cms/user/stadie/DiJetNov4Skim_v1B"

#2011-Daten
#good
#datadir = "/scratch/hh/current/cms/user/stadie/JetRun2011APromptRecoB/merged"
#bad
#datadir = "/scratch/hh/current/cms/user/stadie/JetRun2011APromptRecoA/merged"

#PU-sample Z2PU
#datadirmc="/scratch/hh/current/cms/user/stadie/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1Amerged" 
#Z2 smeared
#datadirmc = "/scratch/hh/current/cms/user/stadie/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Fall10-START38_V12-v1Dsmeared"
#Z2 smeared const term (Mikko)
#datadirmc = "/scratch/hh/current/cms/user/stadie/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Fall10-START38_V12-v1Dsmeared3"
#Z2 smeared
datadirmc = "/scratch/hh/current/cms/user/stadie/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Fall10-START38_V12-v1Dmerged"
#datadirmc = "/scratch/hh/current/cms/user/stadie/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Fall10-START38_V12-v1D"
#D6T
#datadirmc = "/scratch/hh/current/cms/user/stadie/QCD_Pt_15to3000_TuneD6T_Flat_7TeV_pythia6_Fall10-START38_V12-v1Amerged"
#Herwig++
#datadirmc = "/scratch/hh/current/cms/user/stadie/QCD_Pt-15To3000_Tune23_Flat_7TeV-herwigpp_Fall10-START38_V12-v1Amerged"

# Z2 Spring11
#datadirmc = "/scratch/hh/current/cms/user/stadie/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6_Spring11-PU_S1_START311_V1G1-v1Amerged"
#older sample by Matthias datadirmc = "/scratch/hh/current/cms/user/mschrode/mc/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6-Spring11-PU_S1_START311_V1G1-v1"


nthreads = 2
niothreads = 3
nevents =  -1
dirname = DIR_JETALGO+"/dijetsFall10_TuneZ2_AK5"+PF_CALO_JPT+"_weighted_residuals_kostas"
useconstraint = False
batch = False
doBinnedFit = False
doUnbinnedFit = True


#write configs and run the fit
print "producing validation plots for "+jettype+" using "+datadir+" and "+datadirmc+" as MC sample\n"


if not os.path.exists(DIR_JETALGO):
    os.system("mkdir "+DIR_JETALGO)
    
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


