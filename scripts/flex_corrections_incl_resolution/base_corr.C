#define base_corr_cxx
#include "THelpers.h"
#include "base_corr.h"
#include <vector>
#include <TH2.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>
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

  name <<X_labels_[bin_choice] <<Name << "pt_" << setfill('0') <<setw(4)<< pt_bins_[pt_i].first << "_to_" << setfill('0') <<setw(4)<<pt_bins_[pt_i].second;
  name_s=name.str();

TH2D *temp = new TH2D(name_s.c_str(), name_s.c_str()+Title,nbinsx,xlow,xup,nbinsy,ylow,yup);
 temp->Sumw2();
 return temp;

}

TH1D* base_corr::define_pt_histo(TString Name, TString Title, std::vector< std::pair <Double_t,Double_t> > pt_bins_, Int_t pt_i,Int_t nbinsx, Double_t xlow, Double_t xup)
{
  std::stringstream name;
  std::string name_s;

  name <<X_labels_[bin_choice]  <<Name << "pt_" << setfill('0') <<setw(4)<< pt_bins_[pt_i].first << "_to_" << setfill('0') <<setw(4)<<pt_bins_[pt_i].second;
  name_s=name.str();

TH1D *temp = new TH1D(name_s.c_str(), name_s.c_str()+Title,nbinsx,xlow,xup);
 temp->Sumw2();

 return temp;

}

TProfile* base_corr::define_pt_histo(TString Name, TString Title, std::vector< std::pair <Double_t,Double_t> > pt_bins_, Int_t pt_i,Int_t nbinsx, Double_t xlow, Double_t xup, Double_t ylow, Double_t yup)
{
  std::stringstream name;
  std::string name_s;

  name <<X_labels_[bin_choice]  <<Name << "pt_" << setfill('0') <<setw(4)<< pt_bins_[pt_i].first << "_to_" << setfill('0') <<setw(4)<<pt_bins_[pt_i].second;
  name_s=name.str();

TProfile *temp = new TProfile(name_s.c_str(), name_s.c_str()+Title,nbinsx,xlow,xup,ylow,yup);
 temp->Sumw2();

 return temp;

}

THStack* base_corr::define_pt_stack(TString Name, TString Title, std::vector< std::pair <Double_t,Double_t> > pt_bins_, Int_t pt_i)
{
  std::stringstream name;
  std::string name_s;
  //  if(isstack)cout <<"isstack"<<endl;
  name <<X_labels_[bin_choice]  <<Name << "pt_" << setfill('0') <<setw(4)<< pt_bins_[pt_i].first << "_to_" << setfill('0') <<setw(4)<<pt_bins_[pt_i].second;
  name_s=name.str();
  
THStack *temp = new THStack(name_s.c_str(), Title);

 return temp;

}

TString base_corr::define_pt_histo_name(TString Name, std::vector< std::pair <Double_t,Double_t> > pt_bins_, Int_t pt_i)
{
  std::stringstream name;
  std::string name_s;

  name <<X_labels_[bin_choice]  <<Name << "pt_" << setfill('0') <<setw(4)<< pt_bins_[pt_i].first << "_to_" << setfill('0') <<setw(4)<<pt_bins_[pt_i].second;
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
		     rel_sum_class_.push_back(0.);
		     //Tow_EMFt<=0,Tow_EMFt>0<1,Tow_EMFt>=1,rest,Tow_EMFt>=1_inner,Tow_EMFt>=1_outer, 
		     //Tow_EMFt>=1_inner  * Tow_EMFt>=1_outer
		     //--------------------------------------
		     //           Tow_EMFt>=1,
		     //Tow_EMFt>=1 GenJetE, Tow_EMFt>=1_inner0.2, Tow_EMFt>=1_outer0.2
		     //Tow_EMFt>=1 L2L3JetE,
		     //Tow_E(EMF=1-tower with highest energy)	     

		     Double_t jet_sigma_phi = JetEtWeightedSigmaPhi[match];
		     Float_t Jet_phi = JetPhi[match];
		        
		     
		     Double_t highest_energy_EMF1_tower=0.;
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

			    
			    if(Tow_EMFt_.back()==1)
			      {
				if(Tow_E_.back()> highest_energy_EMF1_tower) highest_energy_EMF1_tower=Tow_E_.back();
			      }

                          }
                      }

		    rel_sum_class_[6]=(rel_sum_class_[4] * rel_sum_class_[5]) / (rel_sum_class_[2] );
		    rel_sum_class_[11]=highest_energy_EMF1_tower/L2L3JetE_.back();

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
  std::vector<Double_t> Tow_Eta_;	    
  std::vector<Double_t> Tow_Delta_Phi_;	    
  std::vector<Double_t> Tow_Delta_R_;	    
  std::vector<Double_t> Tow_relETOW_;	    

  std::vector <Double_t > et_sum_;
  std::vector <Double_t > et_weigh_sum_one_;
  std::vector <Double_t > et_weigh_sum_squared_;

  Int_t no_of_variables=7;
  for(Int_t vars_i=0;vars_i<no_of_variables;vars_i++) 
    {
      et_sum_.push_back(0.);
      et_weigh_sum_one_.push_back(0.);
      et_weigh_sum_squared_.push_back(0.);
    }

		        
  Float_t Jet_phi = JetPhi[match];
  Float_t Jet_eta = JetEta[match];

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
	  Tow_Delta_R_.push_back(deltaR (TowEta[tow_i], TowPhi[tow_i], Jet_eta, Jet_phi));
	  //	  cout << Tow_Delta_Phi_.back() << endl;
	  if(Tow_E_.back()!=0)Tow_EMFt_.push_back(Tow_Em_.back()/Tow_E_.back());	
	  Tow_relETOW_.push_back(Tow_E_.back()/JetE_.back());

	  //	  area_norm


	  std::vector <Double_t > value_;
	  value_.push_back(Tow_Delta_Phi_.back());                        //ECAL_SigmaPhi	 		 
	  value_.push_back(Tow_Delta_Phi_.back());			  //HCAL_SigmaPhi    			 
	  value_.push_back(Tow_Delta_Phi_.back());			  //E_SigmaPhi - Standard	 	 
	  value_.push_back(Tow_Delta_Phi_.back() * Tow_EMFt_.back());	  //phi times towemf	 		 
	  value_.push_back(Tow_Delta_Phi_.back());			  //towemf weighted sigma phi	 	 
	  value_.push_back(Tow_Delta_R_.back());			  //ET weighted DeltaR->SigmaR	 
	  value_.push_back(Tow_Delta_R_.back());			  //areanormalized ET weighted DeltaR->SigmaR	 

	  std::vector <Double_t > weight_;
	  weight_.push_back(Tow_Em_.back()  * Tow_Et_over_E_.back());    //ECAL_SigmaPhi	 		 
	  weight_.push_back(Tow_Had_.back()  * Tow_Et_over_E_.back());	 //HCAL_SigmaPhi    			 
	  weight_.push_back(Tow_E_.back()  * Tow_Et_over_E_.back());	 //E_SigmaPhi - Standard	 	 
	  weight_.push_back(Tow_E_.back()  * Tow_Et_over_E_.back());	 //phi times towemf	 		 
	  weight_.push_back(Tow_EMFt_.back());				 //towemf weighted sigma phi	 	 
	  weight_.push_back(Tow_E_.back() * Tow_Et_over_E_.back());	 //ET weighted DeltaR->SigmaR	 
	  weight_.push_back(Tow_E_.back() * Tow_Et_over_E_.back()/deltar_area_norm->GetBinContent(deltar_area_norm->FindBin(Tow_Delta_R_.back())));	 //ET weighted DeltaR->SigmaR	 


	  for(Int_t vars_i=0;vars_i<no_of_variables;vars_i++) 
	    {
	      et_sum_[vars_i]+=weight_[vars_i];
	      et_weigh_sum_one_[vars_i]+=value_[vars_i]*weight_[vars_i];
	      et_weigh_sum_squared_[vars_i]+=value_[vars_i]*value_[vars_i]*weight_[vars_i];
	    }


	}
    }

  for(Int_t i=0;i<no_of_variables;i++)
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
       _L2JetPt=JetPt[match]*JetCorrL2[match];	    
       _L2L3JetPt=JetPt[match]*JetCorrL2L3[match];	    
       _L2L3JetE=JetE[match]*JetCorrL2L3[match];
       _L2JetResponse=_L2JetPt/_GenJetColPt;	      
       _L2L3JetResponse=_L2L3JetPt/_GenJetColPt;	      
       _JetEta=JetEta[match];	    
       _JetPt=JetPt[match];	    
       _JetResponse=_JetPt/_GenJetColPt;	    
       _JetE=JetE[match];	      
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
  X_Var_.push_back(_L2L3JetPt);
 
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
	     Correction_Var_.push_back(1);
	     Correction_Var_.push_back(_JetEta);
	     Correction_Var_.push_back(_JetEMF);
	     Correction_Var_.push_back(_JetEMFCorr);
	     Correction_Var_.push_back(_JetEtWeightedSigmaPhi);
	     Correction_Var_.push_back(_JetEtWeightedSigmaEta);
	     Correction_Var_.push_back((_JetEtWeightedSigmaPhi+_JetEtWeightedSigmaEta)/2.);
	     Correction_Var_.push_back(-0.7485*_JetEtWeightedSigmaPhi+(-0.6631)*_JetEtWeightedSigmaEta);
	     Correction_Var_.push_back(-0.6631*_JetEtWeightedSigmaPhi+(+0.7485)*_JetEtWeightedSigmaEta);
// 2x2 matrix is as follows

//      |      1    |      2    |
// -------------------------------
//    1 |    -0.7485     -0.6631
//    2 |    -0.6631      0.7485
// 	     Correction_Var_.push_back(_JetEtWeightedSigmaPhi*_JetEtWeightedSigmaPhi+_JetEtWeightedSigmaEta*_JetEtWeightedSigmaEta);
 	     Correction_Var_.push_back(TMath::Sqrt(_JetEtWeightedSigmaPhi*_JetEtWeightedSigmaPhi+_JetEtWeightedSigmaEta*_JetEtWeightedSigmaEta));
// 	     Correction_Var_.push_back((_JetEtWeightedSigmaPhi*_JetEtWeightedSigmaPhi+_JetEtWeightedSigmaEta*_JetEtWeightedSigmaEta)/_JetEtWeightedSigmaPhi);
// 	     Correction_Var_.push_back((_JetEtWeightedSigmaPhi+_JetEtWeightedSigmaEta)/_JetEtWeightedSigmaPhi);
// 	     Correction_Var_.push_back(_JetEtWeightedSigmaEta+Tow_EMF_vars_[2]);
// 	     Correction_Var_.push_back(_JetEtWeightedSigmaEta/Tow_sigma_vars_[0]);
// 	     Correction_Var_.push_back(_JetEtWeightedSigmaEta/Tow_sigma_vars_[4]);
////	     Correction_Var_.push_back(Tow_EMF_vars_[11]);
	     Correction_Var_.push_back(Tow_EMF_vars_[7]);
	     Correction_Var_.push_back(Tow_EMF_vars_[10]);
	     Correction_Var_.push_back(Tow_EMF_vars_[2]);
	     Correction_Var_.push_back(Tow_EMF_vars_[0]);
	     Correction_Var_.push_back(Tow_EMF_vars_[1]);
////	     Correction_Var_.push_back(Tow_EMF_vars_[0]+Tow_EMF_vars_[1]);
////	     //	     Correction_Var_.push_back(Tow_EMF_vars_[4]);
// 	     //	     Correction_Var_.push_back(Tow_EMF_vars_[5]);
// 	     //	     Correction_Var_.push_back(Tow_EMF_vars_[6]);
// 	     Correction_Var_.push_back(Tow_EMF_vars_[8]);
// 	     Correction_Var_.push_back(Tow_EMF_vars_[9]);
 	     Correction_Var_.push_back(Tow_sigma_vars_[0]);
 	     Correction_Var_.push_back(Tow_sigma_vars_[1]);
// 	     //	     Correction_Var_.push_back(Tow_sigma_vars_[2]);
// 	     Correction_Var_.push_back(Tow_sigma_vars_[3]);
// 	     //	     Correction_Var_.push_back(Tow_sigma_vars_[4]);
// 	     Correction_Var_.push_back(Tow_sigma_vars_[5]);
// 	     Correction_Var_.push_back(Tow_sigma_vars_[6]);
	     Correction_Var_.push_back((_JetEtWeightedSigmaPhi+_JetEtWeightedSigmaEta)/2.);

	     return Correction_Var_;
 
}



void base_corr::Declare_Labels()
{
  if(statbox==0){
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0000);
  }

   TColor *color2 = gROOT->GetColor(2);
   TColor *color3 = gROOT->GetColor(3);
   TColor *color5 = gROOT->GetColor(5);
   TColor *color20 = gROOT->GetColor(20);
   color20->SetRGB(1,1,0);
   //   TColor *color = gROOT->GetColor(10);

   color5->SetRGB(0.3,0.3,0.3);
   color3->SetRGB(0.35,0.83,0.33);

   //Color:8  Red=0.350000 Green=0.830000 Blue=0.330000 Name=Color8

   //   color5->SetRGB(1.,0.,0.);

   //Color:9  Red=0.350000 Green=0.330000 Blue=0.850000 Name=Color9
   //root [16]

   good_markers.push_back(20);
   good_markers.push_back(24);
   good_markers.push_back(21);
   good_markers.push_back(25);
   good_markers.push_back(22);
   good_markers.push_back(26);
   good_markers.push_back(29);
   good_markers.push_back(30);
   good_markers.push_back(23);
   good_markers.push_back(28);
   good_markers.push_back(34);

img_extension << ".png";

//5 10 12 15 18 22 26 30 35 40 45 51 57 64 72 80 90
// pt_bins_.push_back(make_pair(0,15));
 pt_bins_.push_back(make_pair(15,18 	  ));
 pt_bins_.push_back(make_pair(18,22 	  ));
 pt_bins_.push_back(make_pair(22,26 	  ));
 pt_bins_.push_back(make_pair(26,30 	  ));
 pt_bins_.push_back(make_pair(30,35 	  ));
 pt_bins_.push_back(make_pair(35,40 	  ));
 pt_bins_.push_back(make_pair(40,45 	  ));
 pt_bins_.push_back(make_pair(45,55 	  ));
 pt_bins_.push_back(make_pair(55,64 	  ));
 pt_bins_.push_back(make_pair(64,80 	  ));
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


//traditionally this way
//  pt_bins_.push_back(make_pair(0,15));
//  pt_bins_.push_back(make_pair(15,20 	  ));
//  pt_bins_.push_back(make_pair(20,30 	  ));
//  pt_bins_.push_back(make_pair(30,50 	  ));
//  pt_bins_.push_back(make_pair(50,80 	  ));
//  pt_bins_.push_back(make_pair(80,120	  ));
//  pt_bins_.push_back(make_pair(120,170  ));
//  pt_bins_.push_back(make_pair(170,230  ));
//  pt_bins_.push_back(make_pair(230,300  ));
//  pt_bins_.push_back(make_pair(300,380  ));
//  pt_bins_.push_back(make_pair(380,470  ));
//  pt_bins_.push_back(make_pair(470,600  ));
//  pt_bins_.push_back(make_pair(600,800  ));
//  pt_bins_.push_back(make_pair(800,1000 ));
//  pt_bins_.push_back(make_pair(1000,1400));
//  pt_bins_.push_back(make_pair(1400,1800));
//  pt_bins_.push_back(make_pair(1800,2200));
//  pt_bins_.push_back(make_pair(2200,2600));
//  pt_bins_.push_back(make_pair(2600,3000));
//  pt_bins_.push_back(make_pair(3000,3500));
//  pt_bins_.push_back(make_pair(3500,9999));

 no_small_pt_bins_=pt_bins_.size();
 pt_bins_.push_back(make_pair(50,52));
 pt_bins_.push_back(make_pair(60,62));
 pt_bins_.push_back(make_pair(70,72));
 

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

 L2L3_response_label="p_{T}^{L2L3}/p_{T}^{gen}";
 L2L3_raw_response_label="p_{T}^{raw}/p_{T}^{gen}";
 L2L3_L2_response_label="p_{T}^{L2}/p_{T}^{gen}";


  PDG_labels_.push_back("L_Quarks");
  PDG_labels_.push_back("Gluons");
  PDG_labels_.push_back("C_Quarks");
  PDG_labels_.push_back("B_Quarks");
  PDG_labels_.push_back("other_or_not_matched");
  PDG_labels_.push_back("All_Flav");

  no_PDG_labels_= PDG_labels_.size();
  //  no_PDG_labels_nm_= PDG_labels_.size();


  param_fit_labels_.push_back("Const");
  param_fit_labels_.push_back("X0");
  param_fit_labels_.push_back("B_X-X0_");
  param_fit_labels_.push_back("C_X-X0_2");
  param_fit_labels_.push_back("Sigma");

  param_fit_y_edges_.push_back(make_pair(-0.1,1.1)); 
  param_fit_y_edges_.push_back(make_pair(-0.2,0.4)); 
  param_fit_y_edges_.push_back(make_pair(-2.0,1.0)); 
  param_fit_y_edges_.push_back(make_pair(-10.0,15.0)); 
  param_fit_y_edges_.push_back(make_pair(-.05,0.25)); 


  eta_region_labels_.push_back("barrel_0_0_1_305");
  eta_region_labels_.push_back("tr_endcap_1_305_2_65");
  eta_region_labels_.push_back("HE_2_65_2_964_");
  eta_region_labels_.push_back("Cherenkov_2_964_5_191");
  eta_region_labels_.push_back("All_0_0_5_191");

  eta_region_.push_back(make_pair(0.0,1.305)); 
  eta_region_.push_back(make_pair(1.305,2.65)); 
  eta_region_.push_back(make_pair(2.65,2.964)); 
  eta_region_.push_back(make_pair(2.964,5.191)); 
  eta_region_.push_back(make_pair(0.0,5.191)); 
  no_eta_region_=eta_region_.size();


  X_labels_.push_back("p_{T}^{gen}");
  X_labels_.push_back("E_{T}^{gen}");
  X_labels_.push_back("E^{Gen}");
  X_labels_.push_back("p_{T}^{calo}");
  X_labels_.push_back("E_{T}^{calo}");
  X_labels_.push_back("E^{Calo}");
  X_labels_.push_back("p_{T}^{L2L3}");
  no_X_labels_=X_labels_.size();

  Corr_labels_.push_back(make_pair("no","no"));							
  Corr_labels_.push_back(make_pair("L2","L2"));							
  Corr_labels_.push_back(make_pair("L2L3","L2L3"));							
  Corr_labels_.push_back(make_pair("Eta","#eta"));							
  Corr_labels_.push_back(make_pair("EMF","EMF"));							
  Corr_labels_.push_back(make_pair("CorrEMF","EMF_{Corr}"));						
  Corr_labels_.push_back(make_pair("Phi_Phi-mom","#sigma_{#phi#phi}"));						
  Corr_labels_.push_back(make_pair("Eta_Eta-mom","#sigma_{#eta#eta}"));						
  Corr_labels_.push_back(make_pair("JW_raw","jet width"));			  	
  Corr_labels_.push_back(make_pair("A_decorr_80_120_Sigma_Phi_plus_Sigma_Eta","A_{decorr}"));			  
  Corr_labels_.push_back(make_pair("B_decorr_80_120_Sigma_Phi_plus_Sigma_Eta","B_{decorr}"));			  
  //   Corr_labels_.push_back(make_pair("Sigma_Phi_2_plus_Sigma_Eta_2",""));			  	
   Corr_labels_.push_back(make_pair("SQRT_Sigma_Phi_2_plus_Sigma_Eta_2",""));			  
//   Corr_labels_.push_back(make_pair("Sigma_Phi_2_plus_Sigma_Eta_2_div_by_Sigma_Phi",""));			  
//   Corr_labels_.push_back(make_pair("Sigma_Phi_plus_Sigma_Eta_div_by_Sigma_Phi",""));			  
//   Corr_labels_.push_back(make_pair("Sigma_Phi_plus_Tow_EMF1_Sum",""));			  	
//   Corr_labels_.push_back(make_pair("Sigma_Phi_divided_by_Tow_Sigma_ECAL",""));			  
//   Corr_labels_.push_back(make_pair("Sigma_Phi_divided_by_Tow_EMFt_weighted_Sigma_Phi",""));	
////  Corr_labels_.push_back(make_pair("Tow_EMF1_highest_energy_L2L3",""));					
  Corr_labels_.push_back(make_pair("Tow_EMF1_Sum_GENJET","E_{rel,gen}^{EMF_{tow}=1}"));					
  Corr_labels_.push_back(make_pair("Tow_EMF1_Sum_L2L3JET","E_{rel,L2L3}^{EMF_{tow}=1}"));		  			
  Corr_labels_.push_back(make_pair("Tow_EMF1_Sum","E_{rel}^{EMF_{tow}=1}"));						
  Corr_labels_.push_back(make_pair("Tow_EMF0_Sum","E_{rel}^{EMF_{tow}=0}"));						
  Corr_labels_.push_back(make_pair("Tow_EMF01_Sum","E_{rel}^{EMF_{tow}=0to1}"));					
////  Corr_labels_.push_back(make_pair("Tow_EMF0_and_01_Sum",""));					
////  //  Corr_labels_.push_back(make_pair("Tow_inner_EMF1_Sum",""));				
  //  Corr_labels_.push_back(make_pair("Tow_outer_EMF1_Sum",""));				
  //  Corr_labels_.push_back(make_pair("Tow_inner_outer_comb_EMF1_Sum",""));			
//   Corr_labels_.push_back(make_pair("Tow_inner0.2_EMF1_Sum",""));		  		
//   Corr_labels_.push_back(make_pair("Tow_outer0.2_EMF1_Sum",""));		  		
   Corr_labels_.push_back(make_pair("Tow_Sigma_ECAL","#sigma_{#phi#phi}^{ECAL}"));					
   Corr_labels_.push_back(make_pair("Tow_Sigma_HCAL","#sigma_{#phi#phi}^{HCAL}"));					
//   //  Corr_labels_.push_back(make_pair("Tow_Sigma_class",""));					
//   Corr_labels_.push_back(make_pair("Tow_Sigma_times_EMF",""));					
//   //  Corr_labels_.push_back(make_pair("Tow_EMFt_weighted_Sigma_Phi",""));            		
//   Corr_labels_.push_back(make_pair("Tow_Sigma_R",""));            				
//   Corr_labels_.push_back(make_pair("Tow_Sigma_R_area_normalized_weight",""));                   
  Corr_labels_.push_back(make_pair("JW_Fully_parametrized","JW"));                   

  no_Corr_labels_=Corr_labels_.size();

  for (Int_t Corr_i=0;Corr_i<no_Corr_labels_;Corr_i++)
    {
  Corr_labels_.push_back(make_pair("GMP_" + Corr_labels_[Corr_i].first,"GMP_" + Corr_labels_[Corr_i].second));
    }


  no_all_Corr_labels_=Corr_labels_.size();


  for (Int_t Corr_i=0;Corr_i<no_Corr_labels_;Corr_i++)
    {
  Corr_labels_.push_back(make_pair("Parametrisation_" + Corr_labels_[Corr_i].first,"Parametrisation_" + Corr_labels_[Corr_i].second));
    }

  no_all_Corr_labels_plus_parametrisations_=Corr_labels_.size();

  corr_var_x_edges_.push_back(make_pair(-0.1,1.05));    //  ("no",""));							
  corr_var_x_edges_.push_back(make_pair(-0.1,1.05));    //  ("L2",""));							
  corr_var_x_edges_.push_back(make_pair(-0.1,1.05));    //  ("L2L3",""));							
  corr_var_x_edges_.push_back(make_pair(-5.2,5.2));//  ("Eta",""));							
  corr_var_x_edges_.push_back(make_pair(-0.1,1.05));//  ("EMF",""));							
  corr_var_x_edges_.push_back(make_pair(-1.05,1.05)); //     ("CorrEMF",""));						
  corr_var_x_edges_.push_back(make_pair(-0.1,.35));     // ("Sigma_Phi",""));						
  corr_var_x_edges_.push_back(make_pair(-0.1,.35));     // ("Sigma_Eta",""));						
  corr_var_x_edges_.push_back(make_pair(-0.1,.35));     // ("Sigma_Phi_plus_Sigma_Eta_div_by_2",""));			  	
  corr_var_x_edges_.push_back(make_pair(-0.45,.05));     // ("A_decorr_80_120_Sigma_Phi_plus_Sigma_Eta",""));		
  corr_var_x_edges_.push_back(make_pair(-.2,.2));    //  ("B_decorr_80_120_Sigma_Phi_plus_Sigma_Eta",""));		
//   corr_var_x_edges_.push_back(make_pair(-0.1,.15));     // ("Sigma_Phi_2_plus_Sigma_Eta_2",""));			  	
   corr_var_x_edges_.push_back(make_pair(-0.1,.35));     // ("SQRT_Sigma_Phi_2_plus_Sigma_Eta_2",""));			
//   corr_var_x_edges_.push_back(make_pair(-0.1,0.65));   //   ("Sigma_Phi_2_plus_Sigma_Eta_2_div_by_Sigma_Phi",""));	
//   corr_var_x_edges_.push_back(make_pair(-0.55,2.05));   //   ("Sigma_Phi_plus_Sigma_Eta_div_by_Sigma_Phi",""));		
//   corr_var_x_edges_.push_back(make_pair(-0.1,1.55));   //   ("Sigma_Phi_plus_Tow_EMF1_Sum",""));			  	
  //  corr_var_x_edges_.push_back(make_pair(-1.05,1.05));   // //air("Sigma_Phi_divided_by_Tow_Sigma_ECAL",""));		
  //  corr_var_x_edges_.push_back(make_pair(-1.05,1.05));   // //air("Sigma_Phi_divided_by_Tow_EMFt_weighted_Sigma_Phi",""));	
////  corr_var_x_edges_.push_back(make_pair(-0.1,1.05));    // ("Tow_EMF1_highest_energy_L2L3",""));	
  corr_var_x_edges_.push_back(make_pair(-0.1,1.05));    //  ("Tow_EMF1_Sum_GENJET",""));					
  corr_var_x_edges_.push_back(make_pair(-0.1,1.05));    //  ("Tow_EMF1_Sum_L2L3JET",""));		  			
  corr_var_x_edges_.push_back(make_pair(-0.1,1.05));    //  ("Tow_EMF1_Sum",""));							
  corr_var_x_edges_.push_back(make_pair(-0.1,1.05));    //  ("Tow_EMF0_Sum",""));						
  corr_var_x_edges_.push_back(make_pair(-0.1,1.05));    //  ("Tow_EMF01_Sum",""));					
////  corr_var_x_edges_.push_back(make_pair(-0.1,1.05));    //  ("Tow_EMF0_and_01_Sum",""));
//////   corr_var_x_edges_.push_back(make_pair(-0.1,1.05));    //  pair("Tow_inner_EMF1_Sum",""));				
//   corr_var_x_edges_.push_back(make_pair(-0.1,1.05));    //  pair("Tow_outer_EMF1_Sum",""));				
//   corr_var_x_edges_.push_back(make_pair(-0.1,1.05));    //  pair("Tow_inner_outer_comb_EMF1_Sum",""));			
//   corr_var_x_edges_.push_back(make_pair(-0.1,1.05));    //  ("Tow_inner0.2_EMF1_Sum",""));		  		
//   corr_var_x_edges_.push_back(make_pair(-0.1,1.05));    //  ("Tow_outer0.2_EMF1_Sum",""));		  		
   corr_var_x_edges_.push_back(make_pair(-0.1,.35));     // ("Tow_Sigma_ECAL",""));					
   corr_var_x_edges_.push_back(make_pair(-0.1,.35));     // ("Tow_Sigma_HCAL",""));					
// //  corr_var_x_edges_.push_back(make_pair(-0.1,.5));     // pair("Tow_Sigma_class",""));					
//   corr_var_x_edges_.push_back(make_pair(-1.,0.35));   //	      ("Tow_Sigma_times_EMF",""));					
// //  corr_var_x_edges_.push_back(make_pair(-1.,1.05));   //	      pair("Tow_EMFt_weighted_Sigma_Phi",""));            		
//   corr_var_x_edges_.push_back(make_pair(-0.1,.2));	    //  ("Tow_Sigma_R",""));            				
//   corr_var_x_edges_.push_back(make_pair(-0.1,.2));	    //  ("Tow_Sigma_R_area_normalized_weight",""));                   
  corr_var_x_edges_.push_back(make_pair(-0.1,.35));     // ("Fully parametrized Sigma Phi Eta",""));			  	

  for (Int_t Corr_i=0;Corr_i<no_Corr_labels_;Corr_i++)
    {
  corr_var_x_edges_.push_back(corr_var_x_edges_[Corr_i]);
    }



  double_gauss_labels_.push_back("Mean");
  double_gauss_labels_.push_back("Error_{Mean}");
  double_gauss_labels_.push_back("#sigma_{Rel}");
  double_gauss_labels_.push_back("Error_{#sigma_{Rel}}");
  //  double_gauss_labels_.push_back("#frac{#sigma}{Mean} / #frac{#sigma_{L2L3}}{Mean_{L2L3}}");


  deltar_area_norm= new TH1D("deltar_area_norm", "deltar_area_norm",  s_phi_bins_x,s_phi_xlow,s_phi_xhigh);
  Double_t area, binlowedge, binwidth;

  for(Int_t bin_i=1;bin_i<s_phi_bins_x;bin_i++)
    {
      binlowedge =  deltar_area_norm->GetBinLowEdge(bin_i);
      binwidth =  deltar_area_norm->GetBinWidth(bin_i);
  
      area = TMath::Pi() * (2*binlowedge* binwidth + binwidth * binwidth);
      deltar_area_norm->SetBinContent(bin_i,area);
    }


}

std::vector<Double_t> base_corr::doublefit_gaus(TH1D *histo)
{
 std::vector<Double_t> parameters;

if(histo->Fit("gaus","q")==0)
{//;//,"same");

Double_t mean = histo->GetFunction("gaus")->GetParameter(1);
Double_t sigma = histo->GetFunction("gaus")->GetParameter(2);

Double_t spread=1.5;
Double_t xlow = mean-spread*sigma;
Double_t xhigh = mean + spread*sigma;
 histo->Fit("gaus","q","",xlow,xhigh);//,"same",xlow,xhigh);

 parameters.push_back(histo->GetFunction("gaus")->GetParameter(1)); //mean
 parameters.push_back(histo->GetFunction("gaus")->GetParError(1));  //mean_e
 parameters.push_back(histo->GetFunction("gaus")->GetParameter(2)/histo->GetFunction("gaus")->GetParameter(1)); //sigma/mean
 parameters.push_back(histo->GetFunction("gaus")->GetParError(2)/histo->GetFunction("gaus")->GetParameter(1));  //sigma/mean_e

//      parameters.push_back(0.1);
//      parameters.push_back(0.1);
//      parameters.push_back(0.1);
//      parameters.push_back(0.1);
//      cout << "fit nicht erfolgreich" << endl;
}

else
  {
 parameters.push_back(0.5);
 parameters.push_back(0.5);
 parameters.push_back(0.5);
 parameters.push_back(0.5);
  }

 return parameters;
 ///////////////////////////////////////////////////////ADATP HERE
}


void base_corr::draw_graphs(std::vector <TGraphErrors*> graphs_, Double_t ylow, Double_t yhigh, TLegend *legend, TString PDF_PNG_name, Double_t xlow, Double_t xhigh)
{

  if(graphs_.size()>=1)
    {
graphs_[0]->Draw("ALP");
graphs_[0]->GetXaxis()->SetTitle(graphs_[0]->GetXaxis()->GetTitle() + (TString) "[GeV]");
graphs_[0]->GetYaxis()->SetRangeUser(ylow,yhigh);
graphs_[0]->GetXaxis()->SetRangeUser(xlow,xhigh);
graphs_[0]->Draw("ALP");
    }
  int markers[11]={20,24,21,25,22,26,29,30,23,28,34};

 for(unsigned int a_i=0; a_i<graphs_.size();a_i++)
   {
      graphs_[a_i]->SetLineColor(a_i+1);
      graphs_[a_i]->SetMarkerColor(a_i+1);
      graphs_[a_i]->SetMarkerSize(1.8);
      graphs_[a_i]->SetLineStyle(a_i+1);
      graphs_[a_i]->SetLineWidth(2);
      if(a_i<11)graphs_[a_i]->SetMarkerStyle(markers[a_i]);
      graphs_[a_i]->Draw("LP same"); 
   }
 legend->Draw();

 test->Print(PDF_PNG_name+".root");
 test->Print(PDF_PNG_name+".png");
 test->Print(PDF_PNG_name+".eps");

  test->SetLogx();

 test->Print(PDF_PNG_name+"_logX_" +".root");
 test->Print(PDF_PNG_name+"_logX_" +".png");
 test->Print(PDF_PNG_name+"_logX_" +".eps");

  test->SetLogx(0);


}



TGraphErrors* base_corr::make_graph(std::vector < std::vector < Double_t > >  Double_gauss_, Int_t y_para, TString title, Int_t X_par)
{
  //benoetigt: X_choice und double_gauss_labels_

  Int_t size = Double_gauss_.size();

  Double_t X_array[size];
  Double_t X_e_array[size];
  Double_t Y_array[size];
  Double_t Y_e_array[size];

  //X_choice
  if(X_par==-1)X_par=X_choice;

  for (Int_t a_i =0;a_i<size;a_i++)
    {
      X_array[a_i]=tlj_X_counts_all_[X_par].at(a_i)->GetMean(1);
      X_e_array[a_i]=tlj_X_counts_all_[X_par].at(a_i)->GetMean(11);
      Y_array[a_i]=Double_gauss_[a_i].at(y_para);
      Y_e_array[a_i]=Double_gauss_[a_i].at(y_para+1);
    }

  //vorsicht... alles gerade etwas instabil geworden...
  TGraphErrors* temp = new TGraphErrors(size, X_array, Y_array, X_e_array, Y_e_array);
  temp->SetName(title);
  if(title.Contains("igma_{Rel}_ov")&&(titleprint==0)){temp->SetTitle(";"+X_labels_[X_par]+";"+"#frac{#sigma}{#mu} / #frac{#sigma_{L2L3}}{#mu_{L2L3}}");
  temp->Draw();
  temp->SetLineWidth(2);


  temp->GetXaxis()->SetTitle(X_labels_[X_par]);
  //  temp->GetYaxis()->SetTitle("#frac{#sigma}{Mean} / #frac{#sigma_{L2L3}}{Mean_{L2L3}}");
  temp->GetYaxis()->SetTitleOffset(1.2);
  temp->GetYaxis()->SetTitleSize(0.03);


  }
  else {
    temp->SetTitle (title+";"+X_labels_[X_par]+";"+double_gauss_labels_[y_para]);
  temp->GetXaxis()->SetTitle(X_labels_[X_par]);
  temp->GetYaxis()->SetTitle(double_gauss_labels_[y_para]);

  }

  return temp;
}



TGraphErrors* base_corr::make_graph(std::vector < TString > params_labels_, std::vector < std::vector < Double_t > >  Double_params_, Int_t y_para, TString title, Int_t X_par, Int_t x_offset)
{
  //benoetigt: X_choice und double_gauss_labels_

  Int_t size = Double_params_.size();

  Double_t X_array[size];
  Double_t X_e_array[size];
  Double_t Y_array[size];
  Double_t Y_e_array[size];


  //X_choice
  if(X_par==-1)X_par=X_choice;

  for (Int_t a_i =0;a_i<size;a_i++)
    {
      X_array[a_i]=tlj_X_counts_all_[X_par].at(a_i+x_offset)->GetMean(1);
      X_e_array[a_i]=tlj_X_counts_all_[X_par].at(a_i+x_offset)->GetMean(11);
      Y_array[a_i]=Double_params_[a_i].at(2*y_para);
      Y_e_array[a_i]=Double_params_[a_i].at(2*y_para+1);
    }

  TGraphErrors* temp = new TGraphErrors(size, X_array, Y_array, X_e_array, Y_e_array);
  temp->SetName(title);
  if(titleprint==0)title="";
  temp->Draw();
  temp->SetLineWidth(2);

  temp->SetTitle (title+";"+X_labels_[X_par]+";"+params_labels_[y_para]);
  temp->GetXaxis()->SetTitle(X_labels_[X_par]);
  temp->GetYaxis()->SetTitle(params_labels_[y_para]);


  return temp;
}


