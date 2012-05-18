#! /usr/bin/python


def BinningValues(BINNING,AbsEta):
    if(BINNING=="kostas"):
        binning_values = "-5.191 -3.489 -3.139 -2.964 -2.853 -2.5 -2.411 -2.322 -1.93 -1.479 -1.305 -1.131 -0.957 -0.783 -0.522 -0.261 0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191"
        abs_binning_values = "0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 3.139 3.489 5.191"
    elif(BINNING=="k_HFfix"):
        binning_values = "-5.191 -2.964 -2.853 -2.5 -2.411 -2.322 -1.93 -1.479 -1.305 -1.131 -0.957 -0.783 -0.522 -0.261 0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 5.191"
        abs_binning_values = "0 0.261 0.522 0.783 0.957 1.131 1.305 1.479 1.93 2.322 2.411 2.5 2.853 2.964 5.191"
    elif(BINNING=="JER"):
        binning_values = "-5.2 -3.0 -2.0 -1.5 -1.0 -0.5 0 0.5 1.0 1.5 2.0 3.0 5.2"
        abs_binning_values = "0 0.5 1.0 1.5 2.0 3.0 5.2"
    else:
        print "Defining eta bins failed"

    if(AbsEta==True):
        return abs_binning_values
    else:
        return binning_values


def TriggerNamesThresholds(DATAYEAR,USE_NEW_TRIGGERS_AND_FASTPF,SINGLEJET,jettype,NamesThresholds):
    if(DATAYEAR == "2012"):
        if(DATAYEAR == "2012" and USE_NEW_TRIGGERS_AND_FASTPF==1 and SINGLEJET==1):
            print "Using single jet triggers now... " + str(SINGLEJET) + " " + str(USE_NEW_TRIGGERS_AND_FASTPF) + " " + str(DATAYEAR)
#            fcfg.write("Use single jet triggers = true\n")
            trigger_names = "Di-Jet trigger names = HLT_PFJet40;HLT_PFJet80;HLT_PFJet140;HLT_PFJet200;HLT_PFJet260;HLT_PFJet320;HLT_PFJet400\n"
            #Denis May 2011
            if(jettype == "ak5PF"):
                trigger_thresholds = "Di-Jet trigger thresholds = 72 112 186 255 318 388 472 \n"
            #ad-hoc values...
            if(jettype == "ak5Calo"):
                trigger_thresholds = "Di-Jet trigger thresholds = 50 95 155 220 280 350 430 \n"
            if(jettype == "ak5JPT"):
                trigger_thresholds = "Di-Jet trigger thresholds = 50 95 155 220 280 350 430 \n"
        elif(DATAYEAR == "2012" and USE_NEW_TRIGGERS_AND_FASTPF==1):
            trigger_names = "Di-Jet trigger names = HLT_DiPFJetAve40;HLT_DiPFJetAve80;HLT_DiPFJetAve140;HLT_DiPFJetAve200;HLT_DiPFJetAve260;HLT_DiPFJetAve320;HLT_DiPFJetAve400\n"
            #Denis May 2011
            if(jettype == "ak5PF"):
                trigger_thresholds = "Di-Jet trigger thresholds = 64 107 178 247 313 381 469 \n"
            #WARNING: ad-hoc values...
            if(jettype == "ak5Calo"):
                trigger_thresholds = "Di-Jet trigger thresholds = 50 95 155 220 280 350 430 \n"
            if(jettype == "ak5JPT"):
                trigger_thresholds = "Di-Jet trigger thresholds = 50 95 155 220 280 350 430 \n"




    elif(DATAYEAR == "2011"):
        if(DATAYEAR == "2011" and USE_NEW_TRIGGERS_AND_FASTPF==1 and SINGLEJET==1):
#            fcfg.write("Use single jet triggers = true\n")
            trigger_names = "Di-Jet trigger names = HLT_Jet30;HLT_Jet60;HLT_Jet80;HLT_Jet110;HLT_Jet150;HLT_Jet190;HLT_Jet240;HLT_Jet300;HLT_Jet370\n"
            #determined from full statistics 03/02/2012
            if(jettype == "ak5Calo"):
                trigger_thresholds = "Di-Jet trigger thresholds = 35 65 89 124 165 205 260 325 400\n"
            if(jettype == "ak5PF"):
                trigger_thresholds = "Di-Jet trigger thresholds = 45 80 110 140 185 230 285 350 430 \n"
            if(jettype == "ak5JPT"):
                trigger_thresholds = "Di-Jet trigger thresholds = 40 80 100 135 178 224 280 343 420 \n"
        elif(DATAYEAR == "2011" and USE_NEW_TRIGGERS_AND_FASTPF ==1):
            trigger_names = "Di-Jet trigger names = HLT_DiJetAve30;HLT_DiJetAve60;HLT_DiJetAve80;HLT_DiJetAve110;HLT_DiJetAve150;HLT_DiJetAve190;HLT_DiJetAve240;HLT_DiJetAve300;HLT_DiJetAve370\n"
            if(jettype == "ak7Calo"):
                trigger_thresholds = "Di-Jet trigger thresholds = 35 69 89 120 163 204 256 318 390\n"#Matthias preliminary
            if(jettype == "ak7PF"):
                trigger_thresholds = "Di-Jet trigger thresholds = 40 75 100 135 175 220 273 335 405 \n"#Matthias preliminary
            if(jettype == "ak7JPT"):
                trigger_thresholds = "Di-Jet trigger thresholds = 40 75 100 135 175 220 273 335 405 \n"#Matthias preliminary

            #new default thresholds... on full stats 02/03/2012
            if(jettype == "ak5Calo"):
                trigger_thresholds = "Di-Jet trigger thresholds = 37 65 89 120 165 205 260 320 395\n"
            if(jettype == "ak5PF"):
                trigger_thresholds = "Di-Jet trigger thresholds = 42 78 103 135 180 225 280 345 420 \n"
            if(jettype == "ak5JPT"):
                trigger_thresholds = "Di-Jet trigger thresholds = 42 75 95 130 175 220 275 340 415 \n"
        else:
            trigger_names = "Di-Jet trigger names = HLT_DiJetAve15U;HLT_DiJetAve30U;HLT_DiJetAve50U;HLT_DiJetAve70U;HLT_DiJetAve100U;HLT_DiJetAve140U;HLT_DiJetAve180U;HLT_DiJetAve300U\n"
            if(jettype == "ak5Calo"):
                trigger_thresholds = "Di-Jet trigger thresholds = 38 59 86 111 147 196 249 389\n"#Matthias
            if(jettype == "ak5PF"):
                trigger_thresholds = "Di-Jet trigger thresholds = 43 70 100 127 168 214 279 423 \n"#Matthias
            if(jettype == "ak5JPT"):
                trigger_thresholds = "Di-Jet trigger thresholds = 43 66 96 124 165 220 285 430 \n"#Extrapolation from PF  for last two bins Matthias
            if(jettype == "ak7PF"):
                trigger_thresholds = "Di-Jet trigger thresholds = 47 79 111 140 187 240\n"# CMS AN-2010/371 not updated for 2011
            if(jettype == "ak7Calo"):
                trigger_thresholds = "Di-Jet trigger thresholds = 42 47 97 127 168 223\n"# CMS AN-2010/371 not updated for 2011
    else:
        trigger_names = "Di-Jet trigger names = HLT_DiJetAve15U;HLT_DiJetAve30U;HLT_DiJetAve50U;HLT_DiJetAve70U;HLT_DiJetAve100U;HLT_DiJetAve140U\n"
        if(jettype == "ak5Calo"):
            trigger_thresholds = "Di-Jet trigger thresholds = 38 59 86 111 147 196 \n"
        if(jettype == "ak5PF"):
            trigger_thresholds = "Di-Jet trigger thresholds = 43 70 100 127 168 214 \n"
        if(jettype == "ak5JPT"):
           trigger_thresholds = "Di-Jet trigger thresholds = 43 66 96 124 165 220\n"# CMS AN-2010/371
        if(jettype == "ak7JPT"):
            trigger_thresholds = "Di-Jet trigger thresholds = 46 74 107 138 187 245\n"# CMS AN-2010/371
        if(jettype == "ak7PF"):
            trigger_thresholds = "Di-Jet trigger thresholds = 47 79 111 140 187 240\n"# CMS AN-2010/371 
        if(jettype == "ak7Calo"):
            trigger_thresholds = "Di-Jet trigger thresholds = 42 47 97 127 168 223\n"# CMS AN-2010/371 

    if(NamesThresholds=="names"):
        return trigger_names
    elif(NamesThresholds=="thresholds"):
        return trigger_thresholds
    else:
        print "Defining trigger thresholds failed"


def PUWeightingInfo(DATATYPE,MC_type):
    
    if(DATATYPE=="May10"):
        PU_weighting_info = "PU weighting era = Fall11\n PU weighting histogram = /afs/naf.desy.de/user/k/kirschen/PUDistributions/Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v3.pileup_v2.root \n"
    elif(DATATYPE=="PrReV4"):
        PU_weighting_info = "PU weighting era = Fall11\n PU weighting histogram = /afs/naf.desy.de/user/k/kirschen/PUDistributions/Cert_165088-167913_7TeV_PromptReco_JSON.pileup_v2.root \n"
    elif(DATATYPE=="Aug05"):
        PU_weighting_info = "PU weighting era = Fall11\n PU weighting histogram = /afs/naf.desy.de/user/k/kirschen/PUDistributions/Cert_170249-172619_7TeV_ReReco5Aug_Collisions11_JSON_v2.pileup_v2.root \n"
    elif(DATATYPE=="PrReV6"):
        PU_weighting_info = "PU weighting era = Fall11\n PU weighting histogram = /afs/naf.desy.de/user/k/kirschen/PUDistributions/Cert_172620-173692_PromptReco_JSON.pileup_v2.root \n"
    elif(DATATYPE=="11AReRe"):
        PU_weighting_info = "PU weighting era = Fall11\n PU weighting histogram = /afs/naf.desy.de/user/k/kirschen/PUDistributions/work/PU_2011A.root \n"
    elif(DATATYPE=="11BPrV1"):
        PU_weighting_info = "PU weighting era = Fall11\n PU weighting histogram = /afs/naf.desy.de/user/k/kirschen/PUDistributions/2011B_V1_merged_last_three_files.root \n"
    elif(DATATYPE=="11BReRe"):
        PU_weighting_info = "PU weighting era = Fall11\n PU weighting histogram = /afs/naf.desy.de/user/k/kirschen/PUDistributions/2011B_V1_merged_last_three_files.root \n"
    elif(DATATYPE=="Full2011"):
        PU_weighting_info = "PU weighting era = Fall11\n PU weighting histogram = /afs/naf.desy.de/user/k/kirschen/PUDistributions/alltogether.root \n"
    elif(DATATYPE=="42XFull2011"):
        PU_weighting_info = "PU weighting era = Fall11\n PU weighting histogram = /afs/naf.desy.de/user/k/kirschen/PUDistributions/alltogether.root \n"
    elif(DATATYPE=="Z2wPUsmeared_DMC"):
        PU_weighting_info = "PU weighting era = Fall11\n PU weighting histogram = /afs/naf.desy.de/user/k/kirschen/PUDistributions/EXPORT_shifted_PU-hiust.root \n"
    elif(DATATYPE=="TEST"):
        PU_weighting_info = "PU weighting era = Summer12\n PU weighting histogram = /scratch/hh/current/cms/user/kirschen/PUDistributions/Inclusive/MyDataPileupHistogramObservedAllHLT.root \n PU TruthWeighting = Cert_2012_190456-193336 \n PU TruthWeighting MC distribution = TrueSummer12 \n"
    elif(DATATYPE=="2012_193336" and MC_type=="Z2Star_PUS6S7"):
        PU_weighting_info = "PU TruthWeighting = Cert_2012_190456-193336 \n PU TruthWeighting MC distribution = TruePU_S6PlusS7 \n"
    else:
        print "Defining PU reweighting paths failed"
        PU_weighting_info = "PU weighting era = Flat10\n PU weighting histogram = /afs/naf.desy.de/user/k/kirschen/scratch/2011_06_L2L3_Residuals_42X/PUDist_Cert_160404-163869_7TeV_May10ReReco.root \n"

    return PU_weighting_info


def determineDataDir(DATAYEAR,DATATYPE):
    if(DATAYEAR == "2010"):
        if(DATATYPE=="ReReco"):
            datadir = "/scratch/hh/current/cms/user/stadie/DiJetNov4ReReco_v1C/merged"
    #    if(DATATYPE=="Skim"):
    #        datadir = "/scratch/hh/current/cms/user/stadie/DiJetNov4Skim_v1B"
        if(DATATYPE=="Skim"):
            datadir = "/scratch/hh/current/cms/user/stadie/2010/DiJetNov4Skim_v1C/merged"
    if(DATAYEAR == "2011"):
        if(DATATYPE=="Hpp"):
            datadir = "/scratch/hh/current/cms/user/stadie/2011/v6/QCD_Pt-15to3000_Tune23_Flat_7TeV_herwigpp_Summer11-PU_S3_START42_V11-v2/merged"
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
    #        datadir = "/afs/naf.desy.de/user/k/kirschen/scratch/2011_06_L2L3_Residuals_42X/combine_May10ReReco_and_166861"
    #        datadir = "/afs/naf.desy.de/user/k/kirschen/scratch/2011_06_L2L3_Residuals_42X/combine_May10ReReco_and_167913"
            datadir = "/afs/naf.desy.de/user/k/kirschen/scratch/2011_06_L2L3_Residuals_42X/combine_May10ReReco_and_170307"
        if(DATATYPE=="2fb_ReRe_PrRe"):
            datadir = "/afs/naf.desy.de/user/k/kirschen/scratch/2011_06_L2L3_Residuals_42X/combine_May10ReReco_05Aug_v4_v6_160404-173692"
        if(DATATYPE=="May10_pl_v4"):
            datadir = "/afs/naf.desy.de/user/k/kirschen/scratch/2011_06_L2L3_Residuals_42X/combine_May10ReReco_v4_160404-173692"
        if(DATATYPE=="Aug05_pl_v6"):
            datadir = "/afs/naf.desy.de/user/k/kirschen/scratch/2011_06_L2L3_Residuals_42X/combine_05Aug_v6_160404-173692"
        if(DATATYPE=="May10"):
            datadir = "/scratch/hh/current/cms/user/stadie/2011/Jet2011AMay10ReReco_Cert_160404-163869v2/merged"
        if(DATATYPE=="PrReV4"):
            datadir = "/scratch/hh/current/cms/user/stadie/2011/Jet2011APromptRecoV4_Cert_160404-173692/merged"
        if(DATATYPE=="Aug05"):
            datadir = "/scratch/hh/current/cms/user/stadie/2011/Jet2011A05Aug2011V1_Cert_160404-173692/merged"
        if(DATATYPE=="PrReV6"):
            datadir = "/scratch/hh/current/cms/user/stadie/2011/v8/Jet2011APromptRecoV6_Cert_160404-177515/merged"
        if(DATATYPE=="11BPrV1"):
            datadir = "/scratch/hh/current/cms/user/stadie/2011/v8/Jet2011BPromptRecoV1_Cert_160404-180252/merged"
        if(DATATYPE=="11AReRe"):
            datadir = "/scratch/hh/current/cms/user/kirschen/2011_Jets_v9/Jet2011A08NovReRecoV1_Cert_Nov08_160404-180252/merged"
        if(DATATYPE=="11BReRe"):
            datadir = "/scratch/hh/current/cms/user/kirschen/2011_Jets_v9/Jet2011BReRecoV1_Cert_Nov08_160404-180252/merged"
        if(DATATYPE=="Full2011"):
    #        datadir = "/afs/naf.desy.de/user/k/kirschen/scratch/2011_06_L2L3_Residuals_42X/combine_May10ReReco_05Aug_v4_v6_2011BV1_160404-180252"
    #        datadir = "/afs/naf.desy.de/user/k/kirschen/scratch/2012_01_L2L3Residuals/combine_May10ReReco_05Aug_v4_v6_2011B_160404_180252"
            datadir = "/afs/naf.desy.de/user/k/kirschen/scratch/2012_01_L2L3Residuals/combine_ReReco2011A_2011B_160404_180252"
        if(DATATYPE=="42XFull2011"):
            datadir = "/afs/naf.desy.de/user/k/kirschen/scratch/2012_01_L2L3Residuals/combine_May10ReReco_05Aug_v4_v6_2011B_160404_180252"
        if(DATATYPE=="Z2wPUsmeared_DMC"):
            #outdated#datadir = "/scratch/hh/current/cms/user/kirschen/v8_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Fall11-PU_S6_START42_V14B-v1_with_METcorr_nominal"
            datadir = "/scratch/hh/current/cms/user/kirschen/YOSSOF_JETMET_FIGURES_v8_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Fall11-PU_S6_START42_V14B-v1_with_METcorr_nominal"
        if(DATATYPE=="Z2wPU_DMC"):
            datadir = "/scratch/hh/current/cms/user/stadie/2011/v8/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Fall11-PU_S6_START42_V14B-v1/merged"
        if(DATATYPE=="Z2wPUSu11_DMC"):
            datadir = "/scratch/hh/current/cms/user/stadie/2011/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2/merged"
    
    
    if(DATAYEAR == "2012"):
        if(DATATYPE=="TEST"):
    #        datadir = "/scratch/hh/current/cms/user/kirschen/2012_Jets_v2/Jet2012APromptRecoV1_Cert_2012_190456-191276/merged"
    #        datadir = "/scratch/hh/current/cms/user/kirschen/2012_Jets_v3/Jet2012APromptRecoV1_Cert_2012_190456-191859/merged"
            datadir = "/scratch/hh/current/cms/user/kirschen/2012_Jets_v3/Jet2012APromptRecoV1_Cert_2012_190456-193336/merged"
        if(DATATYPE=="2012_193336"):
            datadir = "/scratch/hh/current/cms/user/kirschen/2012_Jets_v3/Jet2012APromptRecoV1_Cert_2012_190456-193336/merged"

    return datadir


def determineDataDirMC(MC,MC_type):
    
    if(MC == "Su12"):
        if(MC_type=="Z2Star_PUS6S7"):
            datadirmc = "/scratch/hh/current/cms/user/kirschen/2012_Jets_v3/Merge_PUS6_PUS7_QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12-_START52_V9-v1"
        if(MC_type=="Z2Star"):
            datadirmc = "/scratch/hh/current/cms/user/kirschen/2012_Jets_v2/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6_Summer12-PU_S7_START50_V15-v1/merged"
        if(MC_type=="Z2Star52"):
            datadirmc = "/scratch/hh/current/cms/user/kirschen/2012_Jets_v3/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6_Summer12-PU_S7_START52_V9-v1/merged"
        if(MC_type=="Z252"):
            datadirmc = "/scratch/hh/current/cms/user/kirschen/2012_Jets_v3/QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12-PU_S7_START52_V9-v1/merged"
    if(MC == "F11"):
        if(MC_type=="Z2wPU"):
            datadirmc = "/scratch/hh/current/cms/user/stadie/2011/v8/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Fall11-PU_S6_START42_V14B-v1/merged"
        if(MC_type=="Z2wPUsmeared"):
            datadirmc = "/scratch/hh/current/cms/user/kirschen/v8_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Fall11-PU_S6_START42_V14B-v1_with_METcorr_nominal"
        if(MC_type=="Z2wPUsm_Y"):
            datadirmc = "/scratch/hh/current/cms/user/kirschen/YOSSOF_v8_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Fall11-PU_S6_START42_V14B-v1_with_METcorr_nominal"
        if(MC_type=="Z2wPUsm_Y_ip"):
            datadirmc = "/scratch/hh/current/cms/user/kirschen/YOSSOF_FIXHF_and_thresholds_INTERPOLATION_v8_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Fall11-PU_S6_START42_V14B-v1_with_METcorr_nominal"
        if(MC_type=="Z2wPUsm_Y_f"):
            datadirmc = "/scratch/hh/current/cms/user/kirschen/YOSSOF_JETMET_FIGURES_v8_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Fall11-PU_S6_START42_V14B-v1_with_METcorr_nominal"
    
    if(MC == "Su11"):
        if(MC_type=="D6TwPU"):
            datadirmc = "/scratch/hh/current/cms/user/stadie/2011/QCD_Pt-15to3000_TuneD6T_Flat_7TeV-pythia6_Summer11-PU_S3_START42_V11-v2/v5"
        if(MC_type=="Z2wPU"):
            datadirmc = "/scratch/hh/current/cms/user/stadie/2011/QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2/merged"
        if(MC_type=="Hpp"):
            datadirmc = "/scratch/hh/current/cms/user/stadie/2011/QCD_Pt-15to3000_Tune23_Flat_7TeV_herwigpp_Summer11-PU_S3_START42_V11-v2/merged"
        if(MC_type=="Z2wPUsmeared"):
            #changed to with METcor 7.11.11: datadirmc = "/scratch/hh/current/cms/user/kirschen/v6_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2_smeared_Matthias_no_METcorr_nominal"
            datadirmc = "/scratch/hh/current/cms/user/kirschen/v6_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2_smeared_Matthias_with_METcorr_nominal"
    
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
        if(MC_type=="Hpp"):
            datadirmc = "/scratch/hh/current/cms/user/stadie/QCD_Pt-15To3000_Tune23_Flat_7TeV-herwigpp_Fall10-START38_V12-v1Amerged"
    #    if(MC_type=="Z2wPU"):
    return datadirmc



def importDatatypesNewTrigger():
    return ["PrRe62pb","42X_corr","42X_PrRe","42X_combPrRe_ReRe","2fb_ReRe_PrRe","May10_pl_v4","Aug05_pl_v6","May10","PrReV4","Aug05","PrReV6","11BPrV1","Full2011","42XFull2011","Z2wPUsmeared_DMC","Z2wPU_DMC","Z2wPUSu11_DMC","11AReRe","11BReRe","TEST","2012_193336"]
