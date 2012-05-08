#! /usr/bin/python

############################################################################
############################################################################
############################################################################
#### This script automates the extrapolation and subsequent nice        ####
#### plotting starting from Kalibriplots.root-files. It collects        ####
#### resulting .txt-correction files and saves all "final" plots        ####
#### in a folder with the current date.                                 ####
#### Usage:                                                             ####
####      -edit config part below (especially samples)                  ####
####      -execute './run_do_flex_systematicsMPF.py                     ####
####                                                                    ####
#### Simplified workflow of what is happening in this script            ####
####      -do the extrapolation and determination of residuals          ####
####       with several calls of do_flex_extrapol                       ####
####      -read in the resulting .root-files with a (large) number      ####
####       of plotting macros to cut down the number of plots           ####
####      -collect .txt-files with residuals with simpler naming        ####
####                                                                    ####
####                                                                    ####
####                                                                    ####
####                                                                    ####
####                                                                    ####
#### Warning:                                                           ####
####      -for the final plots a number of similar, but semi-independent####
####       scripts is called, some of them "only" for historical reasons####
####       They should be replaced in order to be more                  ####
####       maintainable. Same holds true for do_flex_extrapol-part      ####
####      - An alternative way for extraction of residuals is to use    ####
####       Extrapolation.cc which does a pt-dependent krad  (equivalent ####
####       to resolution) and is much more concise and easy to          ####
####       understand. As it is much more modular and builds partly upon####
####       Kalibri-code it is much safer to adapt to newer needs.       ####
####                                                                    ####
####                                                                    ####
#### Macros called in the script:                                       ####
####                                                                    ####
#### plotMPF_vs_relresponse.C compare_algo_Residuals.C compare_kFSR.C   ####
#### compare_kFSR_em.C compare_PT_dep.C compare_PT_dep_pt_equiv.C       ####
#### compare_Residuals.C compare_Residuals_em.C                         ####
#### compare_Residuals_kFSR_eq_one.C compare_RES_VAL.C                  ####
#### compare_RES_VAL_Data_MC.C compare_Slope_Residuals.C                ####
#### compare_Systematics.C compare_VAL_Residuals.C                      ####
####                                                                    ####
#### compare_Systematics.C has many dependencies and should also be     ####
#### replaced... list:                                                  ####
#### plotHistsAndRatio.C plotRatios.C plotRatios_wFit.C                 ####
#### plotRatios_wAbsUncFit.C plotHistswithFittedSystematics.C           ####
#### plotHistswSomeUncAndAbsUncertainty.C                               ####
############################################################################
############################################################################
############################################################################



import os
import time
import time
from datetime import date

sleep_time=20.
sleep_tens=1

def wait():
    for i in range(sleep_tens):
        print "waiting " +str(i*10) + " of " + str(sleep_tens*10)+" seconds"
        time.sleep(10)
    print "waiting done"

#########################################
#####needed config part begins here######
#########################################

reference_samples_all = ['2012TEST_CORRFinal2011_AK5_MC_Su12Z2Star_kostas_TrueReweighting_MPF_AK5']
nice_labels_reference_samples_all = ['2012']

samples_all = ['2012TEST_CORRFinal2011_AK5_MC_Su12Z2Star_kostas_TrueReweighting_MPF_AK5','2012TEST_CORRFinal2011_AK5_MC_Su12Z2Star_kostas_TrueReweighting_AK5']
nice_labels_samples_all = ['2012','2012_RR']


delete_folders=0
rerun_do_flex=0

choose_binning="kostas"
algo_types_all = ['PF']#,'Calo','JPT']

MEAN_or_GAUSSFITMEAN="Mean"
##MPF_or_rel_response is automatically determined later in the script depending upon whether the name in the samples_list contains MPF or not. If it doesn't, 
MPF_or_rel_response_="rel_response"
pub_style="false"

#########################################
#######config part ends here#############
#########################################



#########################################
#####prepare strings for script-call#####
#########################################


algo_list = "\"" #nice_labels_reference_samples_all[0]

for index_algos, algos in enumerate(algo_types_all):
    if(index_algos==0):
        algo_list = "\"" + algos
    else:
        algo_list = algo_list+"\",\"" + algos

algo_list = algo_list + "\"" #nice_labels_reference_samples_all[0]
  
for index_samples, samples in enumerate(samples_all):
    if(choose_binning=="k_HFfix"):
        nice_labels_samples_all[index_samples] = nice_labels_samples_all[index_samples]+"_HF"

for index_reference_samples, reference_samples in enumerate(reference_samples_all):
    if(choose_binning=="k_HFfix"):
        nice_labels_reference_samples_all[index_reference_samples] = nice_labels_reference_samples_all[index_reference_samples]+"_HF"


comparison_systematics = "\"" #nice_labels_reference_samples_all[0]
comparison_systematics_labels = "\"" #nice_labels_reference_samples_all[0]

for index_samples, samples in enumerate(samples_all):
    if(index_samples==0):
        comparison_systematics = "\"" + samples
        comparison_systematics_labels = "\"" + nice_labels_samples_all[index_samples]
    else:
        comparison_systematics = comparison_systematics + "\",\"" + samples
        comparison_systematics_labels = comparison_systematics_labels + "\",\"" + nice_labels_samples_all[index_samples]

comparison_systematics = comparison_systematics + "\""
comparison_systematics_labels = comparison_systematics_labels + "\""


print comparison_systematics
print comparison_systematics_labels

   
if(delete_folders):
    for samples in samples_all:
        print "prepare " + samples
        if os.path.exists(samples):
            os.system("rm -r " + samples)


#########################################
#####call of do_flex_extrapol       #####
#########################################

for algo_types in algo_types_all:
    print algo_types
    if(rerun_do_flex):
        os.system("rm PTDEPENDENCE_"+choose_binning+"_"+algo_types+".root")
    for index_samples, samples in enumerate(samples_all):
        MPF_or_rel_response_="rel_response"
        if(samples.find('MPF')!=-1):
            print MPF_or_rel_response_ + " ... found MPF: "+ str(samples.find('MPF'))
            MPF_or_rel_response_="MPF"
            print MPF_or_rel_response_
        MEAN_or_GAUSSFITMEAN="Mean"
        if(samples.find('GAUSS')!=-1):
            print  MEAN_or_GAUSSFITMEAN + " ... found GAUSS: "+ str(samples.find('GAUSS'))
            MEAN_or_GAUSSFITMEAN="GaussFitMean"
            print MEAN_or_GAUSSFITMEAN
        if(samples.find('IQMean')!=-1):
            print  MEAN_or_GAUSSFITMEAN + " ... found IQMean: "+ str(samples.find('IQMean'))
            MEAN_or_GAUSSFITMEAN="IQMean"
            print MEAN_or_GAUSSFITMEAN
        print "process " + samples + " for " + algo_types
        os.system("root -l -b -q 'plotMPF_vs_relresponse.C(\"../"+samples+"/dijetsFall10_TuneZ2_AK5"+ algo_types +"_weighted_residuals_"+choose_binning+"/plots/KalibriPlots.root\",\""+ nice_labels_samples_all[index_samples]+"\",\""+algo_types+"\")'")
        if(rerun_do_flex):
            os.system("rm logs/log_do_flex*_" + samples+"_"+choose_binning+"_"+algo_types+".txt")
            print "process " + samples + " for " + algo_types + " step1 "
            #step1: epsexport, root_export=true, use_imported_kFSRAbs=false, NO_easy_mean, NO_fit, export_all_plots=true, NO_kFSR_eq_one
            os.system("root -l -b -q 'run_do_flex_extrapol.C+(\""+algo_types +"\",\"TuneZ2\",\"TuneZ2\",\".eps\",\"true\",\"false\",\""+choose_binning+"\",\"use_no_easy_mean\",\"use_no_fitted_kFSR\",\""+samples+"\",\""+MEAN_or_GAUSSFITMEAN+"\",true,\"kFSR_NOT_eq_one\",\""+MPF_or_rel_response_+"\")' &> logs/log_do_flex1_" + samples+"_"+choose_binning+"_"+algo_types+".txt&")
            os.system("fs flush")
            wait()
#            time.sleep(sleep_time)
            print "process " + samples + " for " + algo_types + " step2 "
            #step2: !!!NO_epsexport, root_export=true, use_imported_kFSRAbs=false, NO_easy_mean, NO_fit, export_all_plots=true, !!!kFSR_eq_one
            os.system("root -l -b -q 'run_do_flex_extrapol.C+(\""+algo_types +"\",\"TuneZ2\",\"TuneZ2\",\"no_img\",\"true\",\"false\", \""+choose_binning+"\",\"use_no_easy_mean\",\"use_no_fitted_kFSR\",\""+samples+"\",\""+MEAN_or_GAUSSFITMEAN+"\",true,\"kFSR_eq_one\",\""+MPF_or_rel_response_+"\")' &> logs/log_do_flex2_" + samples+"_"+choose_binning+"_"+algo_types+".txt&")
            os.system("fs flush")
            wait()
#            time.sleep(sleep_time)
            print "process " + samples + " for " + algo_types + " step3 "
            #step3: NO_epsexport, root_export=true, use_imported_kFSRAbs=false, !!!easy_mean, NO_fit, export_all_plots=true, kFSR_eq_one
            os.system("root -l -b -q 'run_do_flex_extrapol.C+(\""+algo_types +"\",\"TuneZ2\",\"TuneZ2\",\"no_img\",\"true\",\"false\", \""+choose_binning+"\",\"use_easy_mean\",\"use_no_fitted_kFSR\",\""+samples+"\",\""+MEAN_or_GAUSSFITMEAN+"\",true,\"kFSR_eq_one\",\""+MPF_or_rel_response_+"\")' &> logs/log_do_flex3_" + samples+"_"+choose_binning+"_"+algo_types+".txt&")
            os.system("fs flush")
            wait()
#            time.sleep(sleep_time)
            print "process " + samples + " for " + algo_types + " step4.1 "
            #step4_1: NO_epsexport, root_export=true, use_imported_kFSRAbs=false, easy_mean, NO_fit, export_all_plots=true, !!!NO_kFSR_eq_one
            os.system("root -l -b -q 'run_do_flex_extrapol.C+(\""+algo_types +"\",\"TuneZ2\",\"TuneZ2\",\".no_img\",\"true\",\"false\", \""+choose_binning+"\",\"use_easy_mean\",\"use_no_fitted_kFSR\",\""+samples+"\",\""+MEAN_or_GAUSSFITMEAN+"\",true,\"kFSR_NOT_eq_one\",\""+MPF_or_rel_response_+"\")' &> logs/log_do_flex4_1_" + samples+"_"+choose_binning+"_"+algo_types+".txt&")
            os.system("fs flush")
            wait()
#            time.sleep(sleep_time)
            print "process " + samples + " for " + algo_types + " step4.2 "
            #step4_2: !!!epsexport, !!!root_export=false, !!!use_imported_kFSRAbs=true, easy_mean, NO_fit, export_all_plots=true, NO_kFSR_eq_one
            os.system("root -l -b -q 'run_do_flex_extrapol.C+(\""+algo_types +"\",\"TuneZ2\",\"TuneZ2\",\".eps\",\"false\",\"true\", \""+choose_binning+"\",\"use_easy_mean\",\"use_no_fitted_kFSR\",\""+samples+"\",\""+MEAN_or_GAUSSFITMEAN+"\",true,\"kFSR_NOT_eq_one\",\""+MPF_or_rel_response_+"\")' &> logs/log_do_flex4_2_" + samples+"_"+choose_binning+"_"+algo_types+".txt&")
            os.system("fs flush")
            wait()
#            time.sleep(sleep_time)
            print "process " + samples + " for " + algo_types + " step5 "
            print "root -l -b -q 'run_do_flex_extrapol.C+(\""+algo_types +"\",\"TuneZ2\",\"TuneZ2\",\".eps\",\"false\",\"true\", \""+choose_binning+"\",\"use_easy_mean\",\"use_fitted_kFSR\",\""+samples+"\",\""+MEAN_or_GAUSSFITMEAN+"\",false,\"kFSR_NOT_eq_one\",\""+MPF_or_rel_response_+"\")' &> logs/log_do_flex5_" + samples+"_"+choose_binning+"_"+algo_types+".txt"
            #step5: epsexport, root_export=false, use_imported_kFSRAbs=true, easy_mean, !!!fit, !!!export_all_plots=false, NO_kFSR_eq_one
            os.system("root -l -b -q 'run_do_flex_extrapol.C+(\""+algo_types +"\",\"TuneZ2\",\"TuneZ2\",\".eps\",\"false\",\"true\", \""+choose_binning+"\",\"use_easy_mean\",\"use_fitted_kFSR\",\""+samples+"\",\""+MEAN_or_GAUSSFITMEAN+"\",false,\"kFSR_NOT_eq_one\",\""+MPF_or_rel_response_+"\")' &> logs/log_do_flex5_" + samples+"_"+choose_binning+"_"+algo_types+".txt&")
            os.system("fs flush")
            wait()
#            time.sleep(sleep_time)
            print "process " + samples + " for " + algo_types + " step6 "
            #step6: !!!NO_epsexport, root_export=false, use_imported_kFSRAbs=true, !!!NO_easy_mean, fit, export_all_plots=false, NO_kFSR_eq_one
            os.system("root -l -b -q 'run_do_flex_extrapol.C+(\""+algo_types +"\",\"TuneZ2\",\"TuneZ2\",\"no_img\",\"false\",\"true\", \""+choose_binning+"\",\"use_no_easy_mean\",\"use_fitted_kFSR\",\""+samples+"\",\""+MEAN_or_GAUSSFITMEAN+"\",false,\"kFSR_NOT_eq_one\",\""+MPF_or_rel_response_+"\")' &> logs/log_do_flex6_" + samples+"_"+choose_binning+"_"+algo_types+".txt&")
            os.system("fs flush")
            wait()
            print "process " + samples + " for " + algo_types + " step7 "
            #step7: NO_epsexport, root_export=false, use_imported_kFSRAbs=true, NO_easy_mean, !!!NO_fit, export_all_plots=false, NO_kFSR_eq_one
            os.system("root -l -b -q 'run_do_flex_extrapol.C+(\""+algo_types +"\",\"TuneZ2\",\"TuneZ2\",\"no_img\",\"false\",\"true\", \""+choose_binning+"\",\"use_no_easy_mean\",\"use_NO_fitted_kFSR\",\""+samples+"\",\""+MEAN_or_GAUSSFITMEAN+"\",false,\"kFSR_NOT_eq_one\",\""+MPF_or_rel_response_+"\")' &> logs/log_do_flex6_" + samples+"_"+choose_binning+"_"+algo_types+".txt&")
            os.system("fs flush")
            wait()
#            time.sleep(sleep_time)





print algo_types +  " done..."



print "................................"
print "................................"
print "......make nice plots now......."
print "................................"
print "................................"


print "................................"
print "................................"
print ".....wait 10 seconds for AFS...."
print "................................"
print "................................"

time.sleep(sleep_time)

print "................................"
print "................................"
print ".........waiting is done........"
print "................................"
print "................................"


#########################################
#####call of final plotting macros ######
#########################################

os.system("root -l -b -q MakeDateDir.h+")

for algo_types in algo_types_all:
#for algo_types in ['PF']:
    print algo_types
#    for samples in ['2011PrReco_CORRFall10_AK5_MC_F10Z2_SmConst_kostas_AK5', '2010Skim_CORRFall10_AK5_MC_F10Z2_SmConst_kostas_AK5']:
    for index_samples, samples in enumerate(samples_all):
        print "process " + samples + " for " + algo_types
        print index_samples
        for index_reference_samples, reference_samples in enumerate(reference_samples_all):
            #            print "root -b -q 'compare_kFSR.C(\""+reference_samples + "\",\""+samples+"\",\""+algo_types+"\",\"\")'"
            #            print "root -b -q 'compare_Residuals.C(\""+reference_samples + "\",\""+samples+"\",\""+algo_types+"\",\"\",\""+nice_labels_reference_samples_all[index_reference_samples]+"\",\""+nice_labels_samples_all[index_samples]+"\",\""+choose_binning+"\")'"
            os.system("root -b -q 'compare_kFSR.C(\""+reference_samples + "\",\""+samples+"\",\""+algo_types+"\",\"\",\""+nice_labels_reference_samples_all[index_reference_samples]+"\",\""+nice_labels_samples_all[index_samples]+"\",\""+choose_binning+"\","+pub_style+")'")
            os.system("root -b -q 'compare_kFSR_em.C(\""+reference_samples + "\",\""+samples+"\",\""+algo_types+"\",\"\",\""+nice_labels_reference_samples_all[index_reference_samples]+"\",\""+nice_labels_samples_all[index_samples]+"\",\""+choose_binning+"\")'")
            os.system("root -b -q 'compare_Slope_Residuals.C(\""+reference_samples + "\",\""+samples+"\",\""+algo_types+"\",\"\",\""+nice_labels_reference_samples_all[index_reference_samples]+"\",\""+nice_labels_samples_all[index_samples]+"\",\""+choose_binning+"\","+pub_style+")'")
            os.system("root -b -q 'compare_Residuals.C(\""+reference_samples + "\",\""+samples+"\",\""+algo_types+"\",\"\",\""+nice_labels_reference_samples_all[index_reference_samples]+"\",\""+nice_labels_samples_all[index_samples]+"\",\""+choose_binning+"\","+pub_style+")'")
            os.system("root -b -q 'compare_Residuals_kFSR_eq_one.C(\""+reference_samples + "\",\""+samples+"\",\""+algo_types+"\",\"\",\""+nice_labels_reference_samples_all[index_reference_samples]+"\",\""+nice_labels_samples_all[index_samples]+"\",\""+choose_binning+"\")'")
            os.system("root -b -q 'compare_Residuals_em.C(\""+reference_samples + "\",\""+samples+"\",\""+algo_types+"\",\"\",\""+nice_labels_reference_samples_all[index_reference_samples]+"\",\""+nice_labels_samples_all[index_samples]+"\",\""+choose_binning+"\")'")
            os.system("root -b -q 'compare_Residuals_ptdep_no_easymean_nofitkFSR.C(\""+reference_samples + "\",\""+samples+"\",\""+algo_types+"\",\"\",\""+nice_labels_reference_samples_all[index_reference_samples]+"\",\""+nice_labels_samples_all[index_samples]+"\",\""+choose_binning+"\")'")
            os.system("root -b -q 'compare_VAL_Residuals.C(\""+reference_samples + "\",\""+samples+"\",\""+algo_types+"\",\"\",\""+nice_labels_reference_samples_all[index_reference_samples]+"\",\""+nice_labels_samples_all[index_samples]+"\",\""+choose_binning+"\")'")
            os.system("root -b -q 'compare_RES_VAL.C(\""+samples+"\",\""+algo_types+"\",\"\",\""+nice_labels_samples_all[index_samples]+"\",\""+choose_binning+"\")'")
            os.system("root -b -q 'compare_RES_VAL_Data_MC.C(\""+samples+"\",\""+algo_types+"\",\"\",\""+nice_labels_samples_all[index_samples]+"\",\""+choose_binning+"\")'")
            os.system("root -b -q 'compare_Systematics.C+(\""+reference_samples + "\",\""+samples+ "\",\""+nice_labels_reference_samples_all[index_reference_samples]+"\",\""+nice_labels_samples_all[index_samples]+"\",\""+algo_types+"\",\"\",\""+choose_binning+"\","+pub_style+")'")
            print "root -b -q 'compare_Systematics.C+(\""+reference_samples + "\",\""+samples+ "\",\""+nice_labels_reference_samples_all[index_reference_samples]+"\",\""+nice_labels_samples_all[index_samples]+"\",\""+algo_types+"\",\"\",\""+choose_binning+"\","+pub_style+")'"
            os.system("root -b -q 'compare_PT_dep.C(\""+samples+"\",\""+algo_types+"\",\"\",\""+nice_labels_samples_all[index_samples]+"\",\""+choose_binning+"\","+pub_style+")'")
            os.system("root -b -q 'compare_PT_dep_pt_equiv.C(\""+samples+"\",\""+algo_types+"\",\"\",\""+nice_labels_samples_all[index_samples]+"\",\""+choose_binning+"\","+pub_style+")'")
            #            ("+comparison_systematics+","+comparison_systematics_labels+",\""+algo_types+"\",\"\",\""+choose_binning+"\","+pub_style+")'")


print "................................"
print "................................"
print ".......nice plots done.........."
print "................................"
print "................................"

for index_samples, samples in enumerate(samples_all):
#for algo_types in ['PF']:
    print algo_types
#    for samples in ['2011PrReco_CORRFall10_AK5_MC_F10Z2_SmConst_kostas_AK5', '2010Skim_CORRFall10_AK5_MC_F10Z2_SmConst_kostas_AK5']:
    for algo_types in algo_types_all:
        print "process " + samples + " for " + algo_types



d = date.today()
plot_dir_today=d.strftime("20%y_%m_%d")+"_plots"
#strftime("%d/%m/%y")
#2012_02_16_plots
print "................................"
print "................................"
print "..collect Residual .txt-files..."
print "................................"
print "................................"

for algo_types in algo_types_all:
    print algo_types
    os.system("root  -b -q 'compare_Systematics.C+("+comparison_systematics+","+comparison_systematics_labels+",\""+algo_types+"\",\"\",\""+choose_binning+"\","+pub_style+")'")

    print comparison_systematics
    print comparison_systematics_labels

    for index_samples, samples in enumerate(samples_all):
        print "process " + samples + " for " + algo_types
        print "produced residual for " + nice_labels_samples_all[index_samples] + "_L2L3Residual_AK5"+algo_types+".txt"
        os.system("cp " + samples +"/"+choose_binning+"_use_coarse_kFSRAbs_use_easy_mean_TuneZ2_L2L3Residual_AK5"+algo_types+".txt "+plot_dir_today+"/"+ nice_labels_samples_all[index_samples] + "_L2L3Residual_AK5"+algo_types+".txt")
        os.system("cp " + samples +"/"+choose_binning+"_use_coarse_kFSRAbs_TuneZ2_L2L3Residual_AK5"+algo_types+".txt "+plot_dir_today+"/"+ nice_labels_samples_all[index_samples] + "_pt_L2L3Residual_AK5"+algo_types+".txt")
        os.system("cp " + samples +"/"+choose_binning+"_use_coarse_kFSRAbs_use_easy_mean_TuneZ2_L2L3Residual_AK5"+algo_types+"_PTDEP.txt "+plot_dir_today+"/"+ nice_labels_samples_all[index_samples] + "_PTDEP_L2L3Residual_AK5"+algo_types+".txt")
        os.system("cp " + samples +"/"+choose_binning+"_use_coarse_kFSRAbs_TuneZ2_L2L3Residual_AK5"+algo_types+"_PTDEP.txt "+plot_dir_today+"/"+ nice_labels_samples_all[index_samples] + "_PTDEP_pt_L2L3Residual_AK5"+algo_types+".txt")
        os.system("cp " + samples +"/"+choose_binning+"_use_easy_mean_TuneZ2_L2L3Residual_AK5"+algo_types+".txt "+plot_dir_today+"/"+ nice_labels_samples_all[index_samples] + "_ALL_EASY_L2L3Residual_AK5"+algo_types+".txt")
        os.system("cp " + samples +"/"+choose_binning+"_use_coarse_kFSRAbs_use_easy_mean_TuneZ2_L2L3Residual_AK5"+algo_types+".txt "+plot_dir_today+"/"+ nice_labels_samples_all[index_samples] + "_ALL_EASY_imp_kFSRAbs_L2L3Residual_AK5"+algo_types+".txt")
        os.system("cp " + samples +"/"+choose_binning+"_use_coarse_kFSRAbs_use_easy_mean_TuneZ2_Abseta_L2L3Residual_AK5"+algo_types+".txt "+plot_dir_today+"/"+ nice_labels_samples_all[index_samples] + "_ALL_EASY_imp_kFSRAbs_Abseta_L2L3Residual_AK5"+algo_types+".txt")
#        os.system("cp " + samples +"/"+choose_binning+"_use_coarse_kFSRAbs_use_easy_mean_TuneZ2_L2L3Residual_AK5"+algo_types+".txt "+plot_dir_today+"/"+ nice_labels_samples_all[index_samples] + "_ALL_EASY_imp_kFSRAbs_L2L3Residual_AK5"+algo_types+".txt")
       

os.system("cp run_do_flex_select_systematicsMPF.py " + plot_dir_today + "/.")
