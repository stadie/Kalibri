#include "computeTriggerEfficiencies.C"


void runComputeTriggerEfficiencies(){
  //  run(1000000);
  //  run(-1);

  std::vector <TString> inputlist;
  inputlist.push_back("/afs/naf.desy.de/user/k/kirschen/scratch/Kalibri2/L2andJERScripts/filelist_Full2011_L1FastJet_AK5PF");
  inputlist.push_back("/afs/naf.desy.de/user/k/kirschen/scratch/Kalibri2/L2andJERScripts/filelist_Full2011_L1Offset_AK5JPT");
  inputlist.push_back("/afs/naf.desy.de/user/k/kirschen/scratch/Kalibri2/L2andJERScripts/filelist_Full2011_L1Offset_AK5Calo");

  std::vector <double> etacuts;
  etacuts.push_back(0.0);
  etacuts.push_back(1.3);
  etacuts.push_back(3.0);

  std::vector <TString> triggersel;
  triggersel.push_back("DiJetAve");
  triggersel.push_back("Jet");





 
  int nEvts = -1;
  double lumi = 4965;
  TriggerYearLabel="2011";


  for(size_t input_i=0;input_i<inputlist.size();input_i++){
    for(size_t eta_i=0;eta_i<etacuts.size();eta_i++){
      for(size_t triggersel_i=0;triggersel_i<triggersel.size();triggersel_i++){
	single_dijet=triggersel.at(triggersel_i);

	if(single_dijet=="DiJetAve"){
	  single_dijet_ptVarLabel = "p^{ave}_{T}";
	}
	else if(single_dijet=="Jet"){
	  single_dijet_ptVarLabel = "p^{lead}_{T}";
	}
	std::cout << "TESTING" << inputlist.at(input_i) << " " << lumi << " " <<etacuts.at(eta_i) << " " << triggersel.at(triggersel_i) << std::endl;
	std::cout << single_dijet_ptVarLabel << std::endl;
	computeTriggerEfficiences(inputlist.at(input_i),lumi,nEvts,etacuts.at(eta_i));

      }
    }
  }


}
