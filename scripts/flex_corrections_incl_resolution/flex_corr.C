#define flex_corr_cxx

#include "THelpers.h"
#include "base_corr.C"
#include "flex_corr.h"


void flex_corr::Book_Histos()
{

  for (Int_t X_i=0;X_i<no_X_labels_;X_i++)
    {
      std::vector < TH1D* >  tlj_X_counts_;
       for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	 {
	   //       	   tlj_X_counts_.push_back(define_pt_histo( "tlj_counts_"+X_labels_[X_i], ";pt;counts", pt_bins_, pt_i,100,0.75*pt_bins_[pt_i].first,1.25*pt_bins_[pt_i].second));
	   tlj_X_counts_.push_back(define_pt_histo( "tlj_counts_"+(X_labels_[X_i]+"_"), ";"+(X_labels_[X_i]+";counts"), pt_bins_, pt_i,100,0.75*pt_bins_[pt_i].first,1.25*pt_bins_[pt_i].second));
	 }
       tlj_X_counts_all_.push_back(tlj_X_counts_);
    }
  
  for (Int_t Corr_i=0;Corr_i<no_Corr_labels_;Corr_i++)
    {
      std::vector < TH2D* >  tlj_X_response_2D_;
      std::vector < TProfile* >  tlj_X_response_prof_;
       for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	 {
	   tlj_X_response_2D_.push_back(define_pt_histo( "tlj_X_response_2D_"+Corr_labels_[Corr_i].first+"_", ";"+Corr_labels_[Corr_i].first+";L2L3-Response", pt_bins_, pt_i, s_phi_bins_x,s_phi_xlow,s_phi_xhigh,response_bins,response_low, response_high));
	   tlj_X_response_prof_.push_back(define_pt_histo( "tlj_X_response_prof_"+Corr_labels_[Corr_i].first+"_", ";"+Corr_labels_[Corr_i].first+";L2L3-Response", pt_bins_, pt_i, s_phi_bins_x,s_phi_xlow,s_phi_xhigh,response_low, response_high));
	 }
              tlj_X_response_2D_all_.push_back(tlj_X_response_2D_);
              tlj_X_response_prof_all_.push_back(tlj_X_response_prof_);
    }

}

void flex_corr::Write_Histos()
{

  TString root_name_binning= "histos_" + X_labels_[bin_choice] + ".root";

   TFile *outf = new TFile(root_name_binning,"RECREATE");


  for (Int_t Corr_i=0;Corr_i<no_Corr_labels_;Corr_i++)
    {
      std::vector < TH1D* >  tlj_X_response_GMP_mean_;
       for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	 {
	   tlj_X_response_2D_all_[Corr_i].at(pt_i)->FitSlicesY(0,tlj_X_response_2D_all_[Corr_i].at(pt_i)->GetXaxis()->FindBin(s_phi_xlow),tlj_X_response_2D_all_[Corr_i].at(pt_i)->GetXaxis()->FindBin(s_phi_xhigh));
	   TString temp =  (TString) tlj_X_response_2D_all_[Corr_i].at(pt_i)->GetName() + "_1";
	   tlj_X_response_GMP_mean_.push_back((TH1D*)gDirectory->Get(temp.Data()));
	   tlj_X_response_GMP_mean_.back()->GetYaxis()->SetRangeUser(-0.5,2.0);

	 }
              tlj_X_response_GMP_mean_all_.push_back(tlj_X_response_GMP_mean_);
    }



  for (Int_t Corr_i=0;Corr_i<no_Corr_labels_;Corr_i++)
    {
      draw_TH2D_save_PS(img_choice, tlj_X_response_2D_all_[Corr_i],tlj_X_response_2D_all_[Corr_i].at(0)->GetName(),"nice","colz","z1");
      draw_TH1D_save_PS( img_choice,tlj_X_response_GMP_mean_all_[Corr_i],tlj_X_response_GMP_mean_all_[Corr_i].at(0)->GetName() + (TString) "_MEAN_GMP", "nice");
      draw_Colz_Prof_GMP_save_PS(img_choice,tlj_X_response_2D_all_[Corr_i], tlj_X_response_prof_all_[Corr_i], tlj_X_response_GMP_mean_all_[Corr_i], tlj_X_response_GMP_mean_all_[Corr_i].at(0)->GetName() + (TString) "_ALL_Colz_Prof_GMP", "nice", "colz", "z1");

      for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	{
	  tlj_X_response_2D_all_[Corr_i].at(pt_i)->Write();
	  tlj_X_response_prof_all_[Corr_i].at(pt_i)->Write();
	  tlj_X_response_GMP_mean_all_[Corr_i].at(pt_i)->Write();
	}
    }

  for (Int_t X_i=0;X_i<no_X_labels_;X_i++)
    {
      draw_TH1D_save_PS(img_choice, tlj_X_counts_all_[X_i],tlj_X_counts_all_[X_i].at(0)->GetName(), "nice");
       for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	 {
	   tlj_X_counts_all_[X_i].at(pt_i)->Write();
	 }
    }

   outf->Close();


}


void flex_corr::Loop()
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



   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!!
   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!!
   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!!
   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!!

   Declare_Labels();
   Book_Histos();

   
   cout << "Bitte Binningvariable w�hlen (ACHTUNG: I.A. sollte die Binningvariable hier die gleiche sein wie beim Rausschreiben der Korrekturplots (bisher Genpt)!!! Eta-Selektion bisher weiter in GenJetColEta):" <<endl;
   for(Int_t i=0;i<no_X_labels_;i++)  cout << i << ": " << X_labels_[i] << endl;
   cin >>bin_choice;
   fflush(stdin);
   cout << "Es wurde " << bin_choice << ", also " << X_labels_[bin_choice]<< ", ausgewaehlt. Danke..." << endl;

   cout << "Bitte Exportformat f�r Bilder w�hlen ('A' f�r keinen Export)" <<endl;
   cin >>img_choice;
   fflush(stdin);
   cout << "Es wurde " << img_choice << " ausgewaehlt. Danke..." << endl;


   Long64_t nentries = fChain->GetEntriesFast();

   //         nentries=100000;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     if(jentry % 1000 == 0) printf("event %d\n", jentry); //Alle tausend Events was sagen
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;



      //Jet quality variables...
      Bool_t towers_gt_1=0;
      Double_t genjet_jet_deltar=0;


      //loop over two leading jets and 
      //1. calculate jet quality variables
      //2. fill histos with correction variables of good jets in (pt)-bins

          for(Int_t genjet_i =0; genjet_i<2;genjet_i++)//�ber erste zwei colgenjet
	{
      //DeltaR ausrechnen (Gen -Jet Nummerierung und passender CaloJet)
	  genjet_jet_deltar= delta_r(GenJetColEta[genjet_i],JetEta[GenJetColJetIdx[genjet_i]],GenJetColPhi[genjet_i],JetPhi[GenJetColJetIdx[genjet_i]]); 

	  Bool_t isolated=1;//Isolation ueberpruefen
	  for(Int_t caljet_i=0;caljet_i<NobjJet;caljet_i++)
	    {
		  if(caljet_i!=GenJetColJetIdx[genjet_i]&&delta_r(GenJetColEta[genjet_i],JetEta[caljet_i],GenJetColPhi[genjet_i],JetPhi[caljet_i])<iso_max)
		    {
		      isolated=0;
		    }
	    }

	

	  
	  towers_gt_1=0;
	  Int_t no_towers=0;

	  for(Int_t a_i=0;a_i<NobjTow;a_i++) // Ueber alle Tower loopen und alle EM und Had-Tower zum jeweiligen GenJet gehoeren zusammenzaehlen
	    {
	      if(Tow_jetidx[a_i]==GenJetColJetIdx[genjet_i])
		{
		  no_towers++;
		}
	    }
	  if(no_towers>1)towers_gt_1=1;// ueberpruefe explizit, dass Jet mehr als einen Tower hat...


	  if(std::abs(GenJetColEta[genjet_i])<1.1&&genjet_jet_deltar<0.25&&isolated&&towers_gt_1)//demand Eta-selection, deltar (genjet,jet) <.25, isolated and more than one tower per jet
	    {
	      Int_t match=GenJetColJetIdx[genjet_i];


	      //Get variables needed for filling...
       	      base_corr::Fill_obvious_vars(genjet_i, match);
	      std::vector < Double_t > X_Var_=get_X_Vars_(genjet_i,match);
	      std::vector < Double_t > Correction_Var_=get_Correction_Vars_(genjet_i,match);


	      //Fill histos now...
	     for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	       {//Jeder der zwei l. Jet wird dem pt-bin zugeordnet, der zu ihm passt
		 if(X_Var_[bin_choice]>pt_bins_[pt_i].first&&X_Var_[bin_choice]<pt_bins_[pt_i].second)
		   {

		     for (Int_t X_i=0;X_i<no_X_labels_;X_i++)
		       {
			 tlj_X_counts_all_[X_i].at(pt_i)->Fill(X_Var_[X_i]);
		       }
		     for(Int_t Corr_i=0;Corr_i<no_Corr_labels_;Corr_i++)
		       {
			 tlj_X_response_2D_all_[Corr_i].at(pt_i)->Fill(Correction_Var_[Corr_i],_L2L3JetResponse);
			 tlj_X_response_prof_all_[Corr_i].at(pt_i)->Fill(Correction_Var_[Corr_i],_L2L3JetResponse);
		       }
		   }
	       }



	    }



	}



   }


   std::cout << "Funzt...auch hier :)" << std::endl;

   gStyle->SetOptStat(1);
   gStyle->SetOptFit(1);
   gStyle->SetPalette(1);

   test->Range(-1.5,-0.625,13.5,5.625);
   test->SetBorderSize(2);
   test->SetFrameFillColor(0);
//   test->SetGridy();
   test->SetSelected(test);

   Write_Histos();

}

