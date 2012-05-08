#! /usr/bin/python

import os

#################################################################
##usage:                                                       ##
##adapt values accordingly and execute in                      ##
##CalibCore root dir:                                          ##
##./L2andJERScripts/run_batch_dijetresiduals.py                ##
##to start batch execution of different samples                ##
#################################################################

MC="Su12" 
BINNING_list= ['kostas','k_HFfix']
PF_CALO_JPT_list=['Calo','PF','JPT']
CORRECT_JETS_L1="true"
SINGLEJET=0

datadir_mc_list=['/scratch/hh/current/cms/user/kirschen/2012_Jets_v2/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6_Summer12-PU_S7_START50_V15-v1/merged','/scratch/hh/current/cms/user/kirschen/2012_Jets_v3/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6_Summer12-PU_S7_START52_V9-v1/merged']
DIR_JETALGO_list=['50MC_AK5','52MC_AK5']
MC_type_list=['Z2star','Z2star52']
DATATYPE_list=['TEST','TEST']




for index_binning, BINNING in enumerate(BINNING_list):
    for ALGO in PF_CALO_JPT_list:
        for index_samples, samples in enumerate(datadir_mc_list):
            PF_CALO_JPT=ALGO
            print "trying to run " + samples
            log_name = "logs/AUTOMATIC_RUN_BATCH_LOG_AK5_MC_"+MC+MC_type_list[index_samples]+"_"+BINNING+"_"+PF_CALO_JPT+"_"+DATATYPE_list[index_samples]+"_L1" +CORRECT_JETS_L1+"_SINGLEJET" +str(SINGLEJET)+"_"
            os.system("rm " +log_name +DIR_JETALGO_list[index_samples]+ ".txt" )
            print "./scripts/runDiJets_Residuals.py "+DIR_JETALGO_list[index_samples]+" "+PF_CALO_JPT+" " +MC +" " +MC_type_list[index_samples]+" " +BINNING+" " + samples +" " +DATATYPE_list[index_samples]+" " +CORRECT_JETS_L1+" " +str(SINGLEJET)+ " &> " + log_name +DIR_JETALGO_list[index_samples]+ ".txt"
            os.system("./scripts/runDiJets_Residuals.py "+DIR_JETALGO_list[index_samples]+" "+PF_CALO_JPT+" " +MC +" " +MC_type_list[index_samples]+" " +BINNING+" " + samples +" " +DATATYPE_list[index_samples]+" " +CORRECT_JETS_L1+" " +str(SINGLEJET)+ " &> " + log_name +DIR_JETALGO_list[index_samples]+ ".txt" )
            print "finished running " + samples

