#define base_corr_cxx
#include "THelpers.h"
#include "base_corr.h"
#include <vector>
#include <TH2.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include <iomanip>
#include <cmath>


void base_corr::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L base_corr.C
//      Root > base_corr t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   cout << "Funzt..." << endl;
}

TH2D* base_corr::define_pt_histo(TString Name, TString Title, std::vector< std::pair <Double_t,Double_t> > pt_bins_, Int_t pt_i,Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup)
{
  std::stringstream name;
  std::string name_s;

  name <<Name << "pt_" << setfill('0') <<setw(4)<< pt_bins_[pt_i].first << "_to_" << setfill('0') <<setw(4)<<pt_bins_[pt_i].second;
  name_s=name.str();

TH2D *temp = new TH2D(name_s.c_str(), name_s.c_str()+Title,nbinsx,xlow,xup,nbinsy,ylow,yup);

 return temp;

}

TH1D* base_corr::define_pt_histo(TString Name, TString Title, std::vector< std::pair <Double_t,Double_t> > pt_bins_, Int_t pt_i,Int_t nbinsx, Double_t xlow, Double_t xup)
{
  std::stringstream name;
  std::string name_s;

  name <<Name << "pt_" << setfill('0') <<setw(4)<< pt_bins_[pt_i].first << "_to_" << setfill('0') <<setw(4)<<pt_bins_[pt_i].second;
  name_s=name.str();

TH1D *temp = new TH1D(name_s.c_str(), name_s.c_str()+Title,nbinsx,xlow,xup);

 return temp;

}

TProfile* base_corr::define_pt_histo(TString Name, TString Title, std::vector< std::pair <Double_t,Double_t> > pt_bins_, Int_t pt_i,Int_t nbinsx, Double_t xlow, Double_t xup, Double_t ylow, Double_t yup)
{
  std::stringstream name;
  std::string name_s;

  name <<Name << "pt_" << setfill('0') <<setw(4)<< pt_bins_[pt_i].first << "_to_" << setfill('0') <<setw(4)<<pt_bins_[pt_i].second;
  name_s=name.str();

TProfile *temp = new TProfile(name_s.c_str(), name_s.c_str()+Title,nbinsx,xlow,xup,ylow,yup);

 return temp;

}

TString base_corr::define_pt_histo_name(TString Name, std::vector< std::pair <Double_t,Double_t> > pt_bins_, Int_t pt_i)
{
  std::stringstream name;
  std::string name_s;

  name <<Name << "pt_" << setfill('0') <<setw(4)<< pt_bins_[pt_i].first << "_to_" << setfill('0') <<setw(4)<<pt_bins_[pt_i].second;
  name_s=name.str();

  TString temp = (TString) name_s.c_str();

 return temp;

}

   std::vector < Double_t > base_corr::get_EMF_vars_(Int_t genjet_i, Int_t match)
{
		     std::vector<Double_t> JetE_;
		     JetE_.push_back(JetE[match]);	    
		     std::vector<Double_t> L2L3JetE_;
		     L2L3JetE_.push_back(JetCorrL2L3[match] * JetE[match]);	    
		     std::vector<Double_t> GenJetColE_;
		     GenJetColE_.push_back(GenJetColE[genjet_i]);	    
		     std::vector<Double_t> Tow_E_;	    
		     std::vector<Double_t> Tow_Em_;	    
		     std::vector<Double_t> Tow_Had_;	    
		     std::vector<Double_t> Tow_EMFt_;	    
		     std::vector<Double_t> Tow_relETOW_;	    
		     std::vector<Double_t> Tow_relETOW_gen_;	    
		     std::vector<Double_t> Tow_relETOW_L2L3_;	    
		     std::vector<Double_t> Tow_Phi_;	    
		     std::vector<Double_t> Tow_Delta_Phi_;	    

		     std::vector <Double_t > rel_sum_class_;
		     rel_sum_class_.push_back(0.);
		     rel_sum_class_.push_back(0.);
		     rel_sum_class_.push_back(0.);
		     rel_sum_class_.push_back(0.);
		     rel_sum_class_.push_back(0.);
		     rel_sum_class_.push_back(0.);
		     rel_sum_class_.push_back(0.);
		     rel_sum_class_.push_back(0.);
		     rel_sum_class_.push_back(0.);
		     rel_sum_class_.push_back(0.);
		     rel_sum_class_.push_back(0.);
		     //Tow_EMFt<=0,Tow_EMFt>0<1,Tow_EMFt>=1,rest,Tow_EMFt>=1_inner,Tow_EMFt>=1_outer, 
		     //Tow_EMFt>=1_inner  * Tow_EMFt>=1_outer
		     //--------------------------------------
		     //           Tow_EMFt>=1,
		     //Tow_EMFt>=1 GenJetE, Tow_EMFt>=1_inner0.2, Tow_EMFt>=1_outer0.2
		     //Tow_EMFt>=1 L2L3JetE,
	     

		     Double_t jet_sigma_phi = JetEtWeightedSigmaPhi[match];
		     Float_t Jet_phi = JetPhi[match];
		        

                    for(Int_t tow_i=0;tow_i<NobjTow;tow_i++) // Ueber alle Tower loopen und alle EM und Had-Tower zum jeweiligen GenJet gehoeren zusammenzaehlen
                      {
                        if(Tow_jetidx[tow_i]==GenJetColJetIdx[genjet_i])
                          {
			    Tow_Phi_.push_back(TowPhi[tow_i]);
			    Tow_Delta_Phi_.push_back(deltaPhi (TowPhi[tow_i], Jet_phi));

       			    Tow_E_.push_back(TowE[tow_i]);
			    Tow_Em_.push_back(TowEm[tow_i]);
			    Tow_Had_.push_back(TowHad[tow_i]);
			    if(Tow_E_.back()!=0)Tow_EMFt_.push_back(Tow_Em_.back()/(Tow_Em_.back()+Tow_Had_.back()));	
			    Tow_relETOW_.push_back(Tow_E_.back()/JetE_.back());	
			    Tow_relETOW_gen_.push_back(Tow_E_.back()/GenJetColE_.back());
			    Tow_relETOW_L2L3_.push_back(Tow_E_.back()/L2L3JetE_.back());
		    //			    cout << Tow_EMFt_.back() << "and" <<Tow_Em_.back()/Tow_E_.back()<<endl;
			    if(Tow_EMFt_.back()<=0)rel_sum_class_[0]+=Tow_relETOW_.back();
			    else if(Tow_EMFt_.back()>0&&Tow_EMFt_.back()<1)rel_sum_class_[1]+=Tow_relETOW_.back();
			    else if(Tow_EMFt_.back()>=1)
			      {
				rel_sum_class_[2]+=Tow_relETOW_.back();
				rel_sum_class_[7]+=Tow_relETOW_gen_.back();
				rel_sum_class_[10]+=Tow_relETOW_L2L3_.back();
			    if(Tow_Delta_Phi_.back()<(jet_sigma_phi/2))
				rel_sum_class_[4]+=Tow_relETOW_.back();
			    if(Tow_Delta_Phi_.back()>(jet_sigma_phi/2))
				rel_sum_class_[5]+=Tow_relETOW_.back();

			    if(Tow_Delta_Phi_.back()<(0.2))
				rel_sum_class_[8]+=Tow_relETOW_.back();
			    if(Tow_Delta_Phi_.back()>(0.2))
				rel_sum_class_[9]+=Tow_relETOW_.back();



			      }
			    else rel_sum_class_[3]+=Tow_relETOW_.back();

			    
                          }
                      }

		    rel_sum_class_[6]=(rel_sum_class_[4] * rel_sum_class_[5]) / (rel_sum_class_[2] );

		    return rel_sum_class_;
  


}


   std::vector < Double_t > base_corr::get_new_sigma_phi_vars_(Int_t genjet_i, Int_t match)
{
  std::vector < Double_t > sigma_vars_;
  std::vector<Double_t> JetE_;
  JetE_.push_back(JetE[match]);	    
  std::vector<Double_t> Tow_E_;	    
  std::vector<Double_t> Tow_Et_;	    
  std::vector<Double_t> Tow_Et_over_E_;	    
  std::vector<Double_t> Tow_Em_;	    
  std::vector<Double_t> Tow_Had_;	    
  std::vector<Double_t> Tow_EMFt_;	    
  std::vector<Double_t> Tow_Phi_;	    
  std::vector<Double_t> Tow_Delta_Phi_;	    
  std::vector<Double_t> Tow_relETOW_;	    

  std::vector <Double_t > et_sum_;
  et_sum_.push_back(0.); //ECAL_SigmaPhi	 
  et_sum_.push_back(0.); //HCAL_SigmaPhi    
  et_sum_.push_back(0.); //E_SigmaPhi - Standard	 
  et_sum_.push_back(0.); //phi times towemf	 
  et_sum_.push_back(0.); //towemf weighted sigma phi	 

  std::vector <Double_t > et_weigh_sum_one_;
  et_weigh_sum_one_.push_back(0.); //ECAL_SigmaPhi	 
  et_weigh_sum_one_.push_back(0.); //HCAL_SigmaPhi    
  et_weigh_sum_one_.push_back(0.); //E_SigmaPhi - Standard    
  et_weigh_sum_one_.push_back(0.); //phi times towemf    
  et_weigh_sum_one_.push_back(0.); //towemf weighted sigma phi	 

  std::vector <Double_t > et_weigh_sum_squared_;
  et_weigh_sum_squared_.push_back(0.);  //ECAL_SigmaPhi	 
  et_weigh_sum_squared_.push_back(0.);  //HCAL_SigmaPhi    
  et_weigh_sum_squared_.push_back(0.);  //E_SigmaPhi - Standard
  et_weigh_sum_squared_.push_back(0.);  //phi times towemf
  et_weigh_sum_squared_.push_back(0.); //towemf weighted sigma phi	 
		        
  Float_t Jet_phi = JetPhi[match];

  for(Int_t tow_i=0;tow_i<NobjTow;tow_i++) // Ueber alle Tower loopen und alle EM und Had-Tower zum jeweiligen GenJet gehoeren zusammenzaehlen
    {
      if(Tow_jetidx[tow_i]==GenJetColJetIdx[genjet_i])
	{		
	  Tow_E_.push_back(TowE[tow_i]);
	  Tow_Et_.push_back(TowEt[tow_i]);
	  if(Tow_E_.back()>0.)Tow_Et_over_E_.push_back(Tow_Et_.back()/Tow_E_.back());
	  else Tow_Et_over_E_.push_back(0.);
	  Tow_Em_.push_back(TowEm[tow_i]);
	  Tow_Had_.push_back(TowHad[tow_i]);
	  Tow_Phi_.push_back(TowPhi[tow_i]);
	  Tow_Delta_Phi_.push_back(deltaPhi (TowPhi[tow_i], Jet_phi));
	  //	  cout << Tow_Delta_Phi_.back() << endl;
	  if(Tow_E_.back()!=0)Tow_EMFt_.push_back(Tow_Em_.back()/Tow_E_.back());	
	  Tow_relETOW_.push_back(Tow_E_.back()/JetE_.back());

	  et_sum_[0]+=Tow_Em_.back()  * Tow_Et_over_E_.back();
	  et_sum_[1]+=Tow_Had_.back() * Tow_Et_over_E_.back();
	  et_sum_[2]+=Tow_E_.back() * Tow_Et_over_E_.back();
	  et_sum_[3]+=Tow_E_.back() * Tow_Et_over_E_.back();
	  et_sum_[4]+=Tow_EMFt_.back();
			    
	  et_weigh_sum_one_[0]+= Tow_Em_.back()  * Tow_Et_over_E_.back() * Tow_Delta_Phi_.back();
	  et_weigh_sum_one_[1]+= Tow_Had_.back()  * Tow_Et_over_E_.back() * Tow_Delta_Phi_.back();
	  et_weigh_sum_one_[2]+= Tow_E_.back()  * Tow_Et_over_E_.back() * Tow_Delta_Phi_.back();
	  et_weigh_sum_one_[3]+= Tow_E_.back()  * Tow_Et_over_E_.back() * Tow_Delta_Phi_.back() * Tow_EMFt_.back();
	  et_weigh_sum_one_[4]+= Tow_EMFt_.back() * Tow_Delta_Phi_.back();

	  et_weigh_sum_squared_[0]+= Tow_Em_.back()   * Tow_Et_over_E_.back() * Tow_Delta_Phi_.back() * Tow_Delta_Phi_.back();
	  et_weigh_sum_squared_[1]+= Tow_Had_.back()  * Tow_Et_over_E_.back() * Tow_Delta_Phi_.back() * Tow_Delta_Phi_.back();
	  et_weigh_sum_squared_[2]+= Tow_E_.back()  * Tow_Et_over_E_.back() * Tow_Delta_Phi_.back() * Tow_Delta_Phi_.back();
	  et_weigh_sum_squared_[3]+= Tow_E_.back()  * Tow_Et_over_E_.back() * Tow_Delta_Phi_.back() * Tow_Delta_Phi_.back() *
	     Tow_EMFt_.back() * Tow_EMFt_.back();
	  et_weigh_sum_squared_[4]+= Tow_EMFt_.back() * Tow_Delta_Phi_.back() * Tow_Delta_Phi_.back();
	}
    }


  for(unsigned int i=0;i<et_sum_.size();i++)
    {
      Double_t stand_dev=TMath::Sqrt( (et_weigh_sum_squared_[i]/et_sum_[i]) - TMath::Power(et_weigh_sum_one_[i]/et_sum_[i],2) );
      sigma_vars_.push_back(stand_dev);
      //      cout << "standard deviation("<< i  <<"): " << stand_dev << "; SigmaPhi: " << JetEtWeightedSigmaPhi[match] << endl;
    }

  return sigma_vars_;
}


void base_corr::Fill_obvious_vars(Int_t genjet_i, Int_t match)
{

       //Retrieve data for matching jets...
       _GenJetColE=GenJetColE[genjet_i];   
       _GenJetColEt=GenJetColEt[genjet_i];
       _GenJetColPt=GenJetColPt[genjet_i];         
       _L2L3JetPt=JetPt[match]*JetCorrL2L3[match];	    
       _L2L3JetE=JetE[match]*JetCorrL2L3[match];
       _L2L3JetResponse=_L2L3JetPt/_GenJetColPt;	      
       _JetResponse=_JetPt/_GenJetColPt;	    
       _JetE=JetE[match];	      
       _JetPt=JetPt[match];	    
       _JetEt=JetEt[match];	    
       _JetEtWeightedSigmaPhi=JetEtWeightedSigmaPhi[match];
       _JetEtWeightedSigmaEta=JetEtWeightedSigmaEta[match];
       _JetEMF=JetEMF[match];         
       _JetEmE=_JetE*_JetEMF;	    
       _JetCorrEmE=_JetEmE-GenJetColEmE[genjet_i];
       _JetEMFCorr=_JetCorrEmE/_JetE; 	    


}

std::vector < Double_t > base_corr::get_X_Vars_(Int_t genjet_i, Int_t match)
{
  std::vector < Double_t > X_Var_;
  
  X_Var_.push_back(_GenJetColPt);
  X_Var_.push_back(_GenJetColEt);
  X_Var_.push_back(_GenJetColE);
  X_Var_.push_back(_JetPt);
  X_Var_.push_back(_JetEt);
  X_Var_.push_back(_JetE);
 
  return X_Var_;
}

std::vector < Double_t > base_corr::get_Correction_Vars_(Int_t genjet_i, Int_t match)
{

	     //Tower EMF Variables
	     std::vector < Double_t > Tow_EMF_vars_=get_EMF_vars_(genjet_i, match);

	     //new Sigma_Phi_variables
	     std::vector < Double_t > Tow_sigma_vars_=get_new_sigma_phi_vars_(genjet_i, match);

	     std::vector < Double_t > Correction_Var_;
 	     Correction_Var_.push_back(1);                                 
	     Correction_Var_.push_back(1);
	     Correction_Var_.push_back(_JetEMF);
	     Correction_Var_.push_back(_JetEMFCorr);
	     Correction_Var_.push_back(_JetEtWeightedSigmaPhi);
	     Correction_Var_.push_back(_JetEtWeightedSigmaEta);
	     Correction_Var_.push_back(_JetEtWeightedSigmaEta+_JetEtWeightedSigmaPhi);
	     Correction_Var_.push_back(_JetEtWeightedSigmaEta+Tow_EMF_vars_[2]);
	     Correction_Var_.push_back(_JetEtWeightedSigmaEta/Tow_sigma_vars_[0]);
	     Correction_Var_.push_back(_JetEtWeightedSigmaEta/Tow_sigma_vars_[4]);
	     Correction_Var_.push_back(Tow_EMF_vars_[0]);
	     Correction_Var_.push_back(Tow_EMF_vars_[1]);
	     Correction_Var_.push_back(Tow_EMF_vars_[2]);
	     Correction_Var_.push_back(Tow_EMF_vars_[4]);
	     Correction_Var_.push_back(Tow_EMF_vars_[5]);
	     Correction_Var_.push_back(Tow_EMF_vars_[6]);
	     Correction_Var_.push_back(Tow_EMF_vars_[7]);
	     Correction_Var_.push_back(Tow_EMF_vars_[8]);
	     Correction_Var_.push_back(Tow_EMF_vars_[9]);
	     Correction_Var_.push_back(Tow_EMF_vars_[10]);
	     Correction_Var_.push_back(Tow_sigma_vars_[0]);
	     Correction_Var_.push_back(Tow_sigma_vars_[1]);
	     Correction_Var_.push_back(Tow_sigma_vars_[2]);
	     Correction_Var_.push_back(Tow_sigma_vars_[3]);
	     Correction_Var_.push_back(Tow_sigma_vars_[4]);

	     return Correction_Var_;
 
}



void base_corr::Declare_Labels()
{

img_extension << ".gif";

 pt_bins_.push_back(make_pair(0,15));
 pt_bins_.push_back(make_pair(15,20 	  ));
 pt_bins_.push_back(make_pair(20,30 	  ));
 pt_bins_.push_back(make_pair(30,50 	  ));
 pt_bins_.push_back(make_pair(50,80 	  ));
 pt_bins_.push_back(make_pair(80,120	  ));
 pt_bins_.push_back(make_pair(120,170  ));
 pt_bins_.push_back(make_pair(170,230  ));
 pt_bins_.push_back(make_pair(230,300  ));
 pt_bins_.push_back(make_pair(300,380  ));
 pt_bins_.push_back(make_pair(380,470  ));
 pt_bins_.push_back(make_pair(470,600  ));
 pt_bins_.push_back(make_pair(600,800  ));
 pt_bins_.push_back(make_pair(800,1000 ));
 pt_bins_.push_back(make_pair(1000,1400));
 pt_bins_.push_back(make_pair(1400,1800));
 pt_bins_.push_back(make_pair(1800,2200));
 pt_bins_.push_back(make_pair(2200,2600));
 pt_bins_.push_back(make_pair(2600,3000));
 pt_bins_.push_back(make_pair(3000,3500));
 pt_bins_.push_back(make_pair(3500,9999));
 

 cout << "pt_bins: " << endl;
 for (unsigned int i=0;i<pt_bins_.size();i++)
   {
     cout << "Low, High: " << pt_bins_[i].first <<  " , " << pt_bins_[i].second << endl;
   }

 no_pt_bins_=pt_bins_.size();

   testout.setf(ios::fixed);
  testout.precision(0);
  testout.fill('0');
  testout.width(4);

  X_labels_.push_back("P_{T,Gen}");
  X_labels_.push_back("E_{T,Gen}");
  X_labels_.push_back("E_{Gen}");
  X_labels_.push_back("P_{T,Calo}");
  X_labels_.push_back("E_{T,Calo}");
  X_labels_.push_back("E_{Calo}");
  no_X_labels_=X_labels_.size();

  Corr_labels_.push_back(make_pair("no",""));
  Corr_labels_.push_back(make_pair("L2L3",""));
  Corr_labels_.push_back(make_pair("EMF",""));
  Corr_labels_.push_back(make_pair("CorrEMF",""));
  Corr_labels_.push_back(make_pair("Sigma_Phi",""));
  Corr_labels_.push_back(make_pair("Sigma_Eta",""));
  Corr_labels_.push_back(make_pair("Sigma_Phi_plus_Sigma_Eta",""));			  
  Corr_labels_.push_back(make_pair("Sigma_Phi_plus_Tow_EMF1_Sum",""));			  
  Corr_labels_.push_back(make_pair("Sigma_Phi_divided_by_Tow_Sigma_ECAL",""));			  
  Corr_labels_.push_back(make_pair("Sigma_Phi_divided_by_Tow_EMFt_weighted_Sigma_Phi",""));
  Corr_labels_.push_back(make_pair("Tow_EMF0_Sum",""));
  Corr_labels_.push_back(make_pair("Tow_EMF01_Sum",""));
  Corr_labels_.push_back(make_pair("Tow_EMF1_Sum",""));
  Corr_labels_.push_back(make_pair("Tow_inner_EMF1_Sum",""));
  Corr_labels_.push_back(make_pair("Tow_outer_EMF1_Sum",""));
  Corr_labels_.push_back(make_pair("Tow_inner_outer_comb_EMF1_Sum",""));
  Corr_labels_.push_back(make_pair("Tow_EMF1_Sum_GENJET",""));
  Corr_labels_.push_back(make_pair("Tow_inner0.2_EMF1_Sum",""));		  
  Corr_labels_.push_back(make_pair("Tow_outer0.2_EMF1_Sum",""));		  
  Corr_labels_.push_back(make_pair("Tow_EMF1_Sum_L2L3JET",""));		  
  Corr_labels_.push_back(make_pair("Tow_Sigma_ECAL",""));
  Corr_labels_.push_back(make_pair("Tow_Sigma_HCAL",""));
  Corr_labels_.push_back(make_pair("Tow_Sigma_class",""));
  Corr_labels_.push_back(make_pair("Tow_Sigma_times_EMF",""));
  Corr_labels_.push_back(make_pair("Tow_EMFt_weighted_Sigma_Phi",""));            


  no_Corr_labels_=Corr_labels_.size();
  Corr_labels_.push_back(make_pair("GMP_no",""));
  Corr_labels_.push_back(make_pair("GMP_L2L3",""));
  Corr_labels_.push_back(make_pair("GMP_EMF",""));
  Corr_labels_.push_back(make_pair("GMP_CorrEMF",""));
  Corr_labels_.push_back(make_pair("GMP_Sigma_Phi",""));
  Corr_labels_.push_back(make_pair("GMP_Sigma_Eta",""));
  Corr_labels_.push_back(make_pair("GMP_Sigma_Phi_plus_Sigma_Eta",""));			  
  Corr_labels_.push_back(make_pair("GMP_Sigma_Phi_plus_Tow_EMF1_Sum",""));			  
  Corr_labels_.push_back(make_pair("GMP_Sigma_Phi_divided_by_Tow_Sigma_ECAL",""));			  
  Corr_labels_.push_back(make_pair("GMP_Sigma_Phi_divided_by_Tow_EMFt_weighted_Sigma_Phi",""));
  Corr_labels_.push_back(make_pair("GMP_Tow_EMF0_Sum",""));
  Corr_labels_.push_back(make_pair("GMP_Tow_EMF01_Sum",""));
  Corr_labels_.push_back(make_pair("GMP_Tow_EMF1_Sum",""));
  Corr_labels_.push_back(make_pair("GMP_Tow_inner_EMF1_Sum",""));
  Corr_labels_.push_back(make_pair("GMP_Tow_outer_EMF1_Sum",""));
  Corr_labels_.push_back(make_pair("GMP_Tow_inner_outer_comb_EMF1_Sum",""));
  Corr_labels_.push_back(make_pair("GMP_Tow_EMF1_Sum_GENJET",""));
  Corr_labels_.push_back(make_pair("GMP_Tow_inner0.2_EMF1_Sum",""));		  
  Corr_labels_.push_back(make_pair("GMP_Tow_outer0.2_EMF1_Sum",""));		  
  Corr_labels_.push_back(make_pair("GMP_Tow_EMF1_Sum_L2L3JET",""));		  
  Corr_labels_.push_back(make_pair("GMP_Tow_Sigma_ECAL",""));
  Corr_labels_.push_back(make_pair("GMP_Tow_Sigma_HCAL",""));
  Corr_labels_.push_back(make_pair("GMP_Tow_Sigma_class",""));
  Corr_labels_.push_back(make_pair("GMP_Tow_Sigma_times_EMF",""));
  Corr_labels_.push_back(make_pair("GMP_Tow_EMFt_weighted_Sigma_Phi",""));            


  no_all_Corr_labels_=Corr_labels_.size();





  double_gauss_labels_.push_back("Mean");
  double_gauss_labels_.push_back("Error_{Mean}");
  double_gauss_labels_.push_back("#sigma_{Rel}");
  double_gauss_labels_.push_back("Error_{#sigma_{Rel}}");




}

std::vector<Double_t> base_corr::doublefit_gaus(TH1D *histo)
{
  histo->Fit("gaus","q");//,"same");

Double_t mean = histo->GetFunction("gaus")->GetParameter(1);
Double_t sigma = histo->GetFunction("gaus")->GetParameter(2);

Double_t spread=1.5;
Double_t xlow = mean-spread*sigma;
Double_t xhigh = mean + spread*sigma;
 histo->Fit("gaus","q","",xlow,xhigh);//,"same",xlow,xhigh);

 std::vector<Double_t> parameters;
 parameters.push_back(histo->GetFunction("gaus")->GetParameter(1)); //mean
 parameters.push_back(histo->GetFunction("gaus")->GetParError(1));  //mean_e
 parameters.push_back(histo->GetFunction("gaus")->GetParameter(2)/histo->GetFunction("gaus")->GetParameter(1)); //sigma/mean
 parameters.push_back(histo->GetFunction("gaus")->GetParError(2)/histo->GetFunction("gaus")->GetParameter(1));  //sigma/mean_e

 return parameters;
 ///////////////////////////////////////////////////////ADATP HERE
}


