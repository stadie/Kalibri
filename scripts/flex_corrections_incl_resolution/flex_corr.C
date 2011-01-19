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

	   if(titleprint==0){
	     tlj_X_counts_.back()->SetTitle("");
	   }

	 }
       tlj_X_counts_all_.push_back(tlj_X_counts_);
    }
  
  for (Int_t Corr_i=0;Corr_i<no_Corr_labels_;Corr_i++)
    {
      std::vector < TH1D* >  tlj_Corr_Vars_counts_;
      std::vector < TH2D* >  tlj_X_response_2D_;
      std::vector < TProfile* >  tlj_X_response_prof_;
      std::vector < TH2D* >  tlj_X_raw_response_2D_;
      std::vector < TProfile* >  tlj_X_raw_response_prof_;
      std::vector < TH2D* >  tlj_X_L2_response_2D_;
      std::vector < TProfile* >  tlj_X_L2_response_prof_;
       for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	 {
	   tlj_Corr_Vars_counts_.push_back(define_pt_histo("tlj_Corr_Vars_counts_"+Corr_labels_[Corr_i].first+"_", ";"+Corr_labels_[Corr_i].second+";counts",pt_bins_,pt_i, s_phi_bins_x,corr_var_x_edges_[Corr_i].first,corr_var_x_edges_[Corr_i].second));
	   tlj_X_response_2D_.push_back(define_pt_histo( "tlj_X_response_2D_"+Corr_labels_[Corr_i].first+"_", ";"+Corr_labels_[Corr_i].second+";"+L2L3_response_label, pt_bins_, pt_i, s_phi_bins_x,corr_var_x_edges_[Corr_i].first,corr_var_x_edges_[Corr_i].second,response_bins,response_low, response_high));
	   tlj_X_response_prof_.push_back(define_pt_histo( "tlj_X_response_prof_"+Corr_labels_[Corr_i].first+"_", ";"+Corr_labels_[Corr_i].second+";"+L2L3_response_label, pt_bins_, pt_i, s_phi_bins_x,corr_var_x_edges_[Corr_i].first,corr_var_x_edges_[Corr_i].second,response_low, response_high));
	   tlj_X_raw_response_2D_.push_back(define_pt_histo( "tlj_X_raw_response_2D_"+Corr_labels_[Corr_i].first+"_", ";"+Corr_labels_[Corr_i].second+";"+L2L3_raw_response_label, pt_bins_, pt_i, s_phi_bins_x,corr_var_x_edges_[Corr_i].first,corr_var_x_edges_[Corr_i].second,response_bins,0, response_high));
	   tlj_X_raw_response_prof_.push_back(define_pt_histo( "tlj_X_raw_response_prof_"+Corr_labels_[Corr_i].first+"_", ";"+Corr_labels_[Corr_i].second+";"+L2L3_raw_response_label, pt_bins_, pt_i, s_phi_bins_x,corr_var_x_edges_[Corr_i].first,corr_var_x_edges_[Corr_i].second,0, response_high));
	   tlj_X_L2_response_2D_.push_back(define_pt_histo( "tlj_X_L2_response_2D_"+Corr_labels_[Corr_i].first+"_", ";"+Corr_labels_[Corr_i].second+";"+L2L3_L2_response_label, pt_bins_, pt_i, s_phi_bins_x,corr_var_x_edges_[Corr_i].first,corr_var_x_edges_[Corr_i].second,response_bins,0, response_high));
	   tlj_X_L2_response_prof_.push_back(define_pt_histo( "tlj_X_L2_response_prof_"+Corr_labels_[Corr_i].first+"_", ";"+Corr_labels_[Corr_i].second+";"+L2L3_L2_response_label, pt_bins_, pt_i, s_phi_bins_x,corr_var_x_edges_[Corr_i].first,corr_var_x_edges_[Corr_i].second,0, response_high));
	   if(titleprint==0){
	     tlj_Corr_Vars_counts_.back()->SetTitle("");
	     tlj_X_response_2D_.back()->SetTitle("");
	     tlj_X_response_prof_.back()->SetTitle("");
	     tlj_X_raw_response_2D_.back()->SetTitle("");
	     tlj_X_raw_response_prof_.back()->SetTitle("");
	     tlj_X_L2_response_2D_.back()->SetTitle("");
	     tlj_X_L2_response_prof_.back()->SetTitle("");
	   }


	 }
       tlj_Corr_Vars_counts_all_.push_back(tlj_Corr_Vars_counts_);
       tlj_X_response_2D_all_.push_back(tlj_X_response_2D_);
       tlj_X_response_prof_all_.push_back(tlj_X_response_prof_);
       tlj_X_raw_response_2D_all_.push_back(tlj_X_raw_response_2D_);
       tlj_X_raw_response_prof_all_.push_back(tlj_X_raw_response_prof_);
       tlj_X_L2_response_2D_all_.push_back(tlj_X_L2_response_2D_);
       tlj_X_L2_response_prof_all_.push_back(tlj_X_L2_response_prof_);
    }

}

void flex_corr::Write_Histos()
{

  TString root_name_binning= "histos_" + eta_region_labels_[eta_choice] + X_labels_[bin_choice] + ".root";

   TFile *outf = new TFile(root_name_binning,"RECREATE");


  for (Int_t Corr_i=0;Corr_i<no_Corr_labels_;Corr_i++)
    {
      std::vector < TH1D* >  tlj_X_response_GMP_mean_;
      std::vector < TH1D* >  tlj_X_response_GMP_sigma_;
      std::vector < TH1D* >  tlj_X_response_GMP_rel_sigma_;
      std::vector < TH1D* >  tlj_X_raw_response_GMP_mean_;
      std::vector < TH1D* >  tlj_X_L2_response_GMP_mean_;
       for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	 {
	   tlj_X_response_2D_all_[Corr_i].at(pt_i)->FitSlicesY(0,0,-1,cut_fitslices);//tlj_X_response_2D_all_[Corr_i].at(pt_i)->GetXaxis()->FindBin(corr_var_x_edges_[Corr_i].first),tlj_X_response_2D_all_[Corr_i].at(pt_i)->GetXaxis()->FindBin(corr_var_x_edges_[Corr_i].second));
	   TString temp =  (TString) tlj_X_response_2D_all_[Corr_i].at(pt_i)->GetName() + "_1";
	   tlj_X_response_GMP_mean_.push_back((TH1D*)gDirectory->Get(temp.Data()));
	   tlj_X_response_GMP_mean_.back()->GetYaxis()->SetRangeUser(-0.5,2.0);


	   temp =  (TString) tlj_X_response_2D_all_[Corr_i].at(pt_i)->GetName() + "_2";
	   tlj_X_response_GMP_sigma_.push_back((TH1D*)gDirectory->Get(temp.Data()));
	   tlj_X_response_GMP_sigma_.back()->GetYaxis()->SetRangeUser(-0.1,0.5);

	   tlj_X_response_GMP_rel_sigma_.push_back((TH1D*)tlj_X_response_GMP_sigma_.back()->Clone());
	   tlj_X_response_GMP_rel_sigma_.back()->Divide(tlj_X_response_GMP_sigma_.back(),tlj_X_response_GMP_mean_.back());
	   tlj_X_response_GMP_rel_sigma_.back()->GetXaxis()->SetRangeUser(0,0.25);
	   tlj_X_response_GMP_rel_sigma_.back()->GetYaxis()->SetRangeUser(-0.1,0.5);
	   tlj_X_response_GMP_rel_sigma_.back()->GetYaxis()->SetTitle("#sigma/#mu");
	   tlj_X_response_GMP_rel_sigma_.back()->SetStats(0);
	   tlj_X_response_GMP_rel_sigma_.back()->SetMarkerStyle(20);
	   tlj_X_response_GMP_rel_sigma_.back()->SetMarkerColor(2);
	   tlj_X_response_GMP_rel_sigma_.back()->SetLineColor(2);
	   tlj_X_response_GMP_rel_sigma_.back()->SetTitle("");
	   //	   TPaveStats *st = (TPaveStats*)tlj_X_response_GMP_rel_sigma_.back()->GetListOfFunctions()->FindObject("stats")


	   tlj_X_raw_response_2D_all_[Corr_i].at(pt_i)->FitSlicesY(0,0,-1,cut_fitslices);//tlj_X_response_2D_all_[Corr_i].at(pt_i)->GetXaxis()->FindBin(corr_var_x_edges_[Corr_i].first),tlj_X_response_2D_all_[Corr_i].at(pt_i)->GetXaxis()->FindBin(corr_var_x_edges_[Corr_i].second));
	   temp =  (TString) tlj_X_raw_response_2D_all_[Corr_i].at(pt_i)->GetName() + "_1";
	   tlj_X_raw_response_GMP_mean_.push_back((TH1D*)gDirectory->Get(temp.Data()));
	   tlj_X_raw_response_GMP_mean_.back()->GetYaxis()->SetRangeUser(-0.5,2.0);

	   tlj_X_L2_response_2D_all_[Corr_i].at(pt_i)->FitSlicesY(0,0,-1,cut_fitslices);//tlj_X_response_2D_all_[Corr_i].at(pt_i)->GetXaxis()->FindBin(corr_var_x_edges_[Corr_i].first),tlj_X_response_2D_all_[Corr_i].at(pt_i)->GetXaxis()->FindBin(corr_var_x_edges_[Corr_i].second));
	   temp =  (TString) tlj_X_L2_response_2D_all_[Corr_i].at(pt_i)->GetName() + "_1";
	   tlj_X_L2_response_GMP_mean_.push_back((TH1D*)gDirectory->Get(temp.Data()));
	   tlj_X_L2_response_GMP_mean_.back()->GetYaxis()->SetRangeUser(-0.5,2.0);

	 }
              tlj_X_response_GMP_mean_all_.push_back(tlj_X_response_GMP_mean_);
              tlj_X_response_GMP_sigma_all_.push_back(tlj_X_response_GMP_sigma_);
              tlj_X_response_GMP_rel_sigma_all_.push_back(tlj_X_response_GMP_rel_sigma_);


              tlj_X_raw_response_GMP_mean_all_.push_back(tlj_X_raw_response_GMP_mean_);
              tlj_X_L2_response_GMP_mean_all_.push_back(tlj_X_L2_response_GMP_mean_);
    }


    TString saveFolder = "Corr_Plots";
    if(chdir(saveFolder) != 0){
      mkdir(saveFolder, S_IRWXU|S_IRWXG|S_IRWXO);
      chdir(saveFolder);
    }
 


  for (Int_t Corr_i=0;Corr_i<no_Corr_labels_;Corr_i++)
    {
      draw_TH2D_save_PS("PS_Corr_Vars_vs_Response",false,img_choice, tlj_X_response_2D_all_[Corr_i], eta_region_labels_[eta_choice] + X_labels_[bin_choice] +tlj_X_response_2D_all_[Corr_i].at(0)->GetName(),"nice","colz","z1");
      draw_TH2D_save_PS("PS_Corr_Vars_vs_Response",false,img_choice, tlj_X_response_2D_all_[Corr_i], eta_region_labels_[eta_choice] + X_labels_[bin_choice] +tlj_X_response_2D_all_[Corr_i].at(0)->GetName(),"nice","colz","z0");
      draw_TH1D_save_PS("PS_prof_and_GMP_all_Corrs", img_choice,tlj_X_response_GMP_mean_all_[Corr_i], eta_region_labels_[eta_choice] + X_labels_[bin_choice] + tlj_X_response_GMP_mean_all_[Corr_i].at(0)->GetName() + (TString) "_MEAN_GMP", "nice");
      draw_TH1D_save_PS("PS_GMP_prof_sigma_all_Corrs", img_choice,tlj_X_response_GMP_sigma_all_[Corr_i], eta_region_labels_[eta_choice] + X_labels_[bin_choice] + tlj_X_response_GMP_sigma_all_[Corr_i].at(0)->GetName() + (TString) "_MEAN_GMP_sigma", "nice");
      draw_TH1D_save_PS("PS_GMP_prof_rel_sigma_all_Corrs", img_choice,tlj_X_response_GMP_rel_sigma_all_[Corr_i], eta_region_labels_[eta_choice] + X_labels_[bin_choice] + tlj_X_response_GMP_rel_sigma_all_[Corr_i].at(0)->GetName() + (TString) "_MEAN_GMP_rel_sigma", "nice");


      draw_Colz_Prof_GMP_save_PS("PS_Corr_Vars_vs_Response",img_choice,tlj_X_response_2D_all_[Corr_i], tlj_X_response_prof_all_[Corr_i], tlj_X_response_GMP_mean_all_[Corr_i],  eta_region_labels_[eta_choice] + X_labels_[bin_choice] + tlj_X_response_GMP_mean_all_[Corr_i].at(0)->GetName() + (TString) "_ALL_Colz_Prof_GMP", "nice", "colz", "z1");
      draw_Colz_Prof_GMP_save_PS("PS_Corr_Vars_vs_Response",img_choice,tlj_X_response_2D_all_[Corr_i], tlj_X_response_prof_all_[Corr_i], tlj_X_response_GMP_mean_all_[Corr_i],  eta_region_labels_[eta_choice] + X_labels_[bin_choice] + tlj_X_response_GMP_mean_all_[Corr_i].at(0)->GetName() + (TString) "_ALL_Colz_Prof_GMP", "nice", "colz", "z0");

      draw_TH1D_save_PS("PS_Corr_Var_Counts",img_choice, tlj_Corr_Vars_counts_all_[Corr_i], eta_region_labels_[eta_choice] + X_labels_[bin_choice] +tlj_Corr_Vars_counts_all_[Corr_i].at(0)->GetName(), "nice","hist","y1");
      draw_TH1D_save_PS("PS_Corr_Var_Counts",img_choice, tlj_Corr_Vars_counts_all_[Corr_i], eta_region_labels_[eta_choice] + X_labels_[bin_choice] +tlj_Corr_Vars_counts_all_[Corr_i].at(0)->GetName(), "nice","hist","y0");


      draw_TH2D_save_PS("PS_Corr_Vars_vs_RAW_Response",false,img_choice, tlj_X_raw_response_2D_all_[Corr_i], eta_region_labels_[eta_choice] + X_labels_[bin_choice] +tlj_X_raw_response_2D_all_[Corr_i].at(0)->GetName(),"nice","colz","z0",0,1.2);
      draw_Colz_Prof_GMP_save_PS("PS_Corr_Vars_vs_RAW_Response",img_choice,tlj_X_raw_response_2D_all_[Corr_i], tlj_X_raw_response_prof_all_[Corr_i], tlj_X_raw_response_GMP_mean_all_[Corr_i],  eta_region_labels_[eta_choice] + X_labels_[bin_choice] + tlj_X_raw_response_GMP_mean_all_[Corr_i].at(0)->GetName() + (TString) "_ALL_Colz_Prof_GMP", "nice", "colz", "z0",0,-1,0,1.2);

      draw_TH2D_save_PS("PS_Corr_Vars_vs_L2_Response",false,img_choice, tlj_X_L2_response_2D_all_[Corr_i], eta_region_labels_[eta_choice] + X_labels_[bin_choice] +tlj_X_L2_response_2D_all_[Corr_i].at(0)->GetName(),"nice","colz","z0",0,1.2);
      draw_Colz_Prof_GMP_save_PS("PS_Corr_Vars_vs_L2_Response",img_choice,tlj_X_L2_response_2D_all_[Corr_i], tlj_X_L2_response_prof_all_[Corr_i], tlj_X_L2_response_GMP_mean_all_[Corr_i],  eta_region_labels_[eta_choice] + X_labels_[bin_choice] + tlj_X_L2_response_GMP_mean_all_[Corr_i].at(0)->GetName() + (TString) "_ALL_Colz_Prof_GMP", "nice", "colz", "z0",0,-1,0,1.2);


      for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	{
	  tlj_Corr_Vars_counts_all_[Corr_i].at(pt_i)->Write();
	  tlj_X_response_2D_all_[Corr_i].at(pt_i)->Write();
	  tlj_X_response_prof_all_[Corr_i].at(pt_i)->Write();
	  tlj_X_response_GMP_mean_all_[Corr_i].at(pt_i)->Write();
	}
    }

  for (Int_t X_i=0;X_i<no_X_labels_;X_i++)
    {
      draw_TH1D_save_PS("PS_X_Counts",img_choice, tlj_X_counts_all_[X_i],  eta_region_labels_[eta_choice] + X_labels_[bin_choice] +tlj_X_counts_all_[X_i].at(0)->GetName(), "nice");
       for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	 {
	   tlj_X_counts_all_[X_i].at(pt_i)->Write();
	 }
    }

  chdir("..");



   outf->Close();


}


void flex_corr::Loop(Bool_t setvalues, Int_t par_bin_choice, Int_t par_eta_choice, TString par_img_choice)
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

   if(setvalues){
 bin_choice	= par_bin_choice;	
 eta_choice	= par_eta_choice;	
 img_choice	= par_img_choice;
   }


   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!!
   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!!
   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!!
   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!!

   Declare_Labels();

   
   cout << "Bitte Binningvariable wählen (ACHTUNG: I.A. sollte die Binningvariable hier die gleiche sein wie beim Rausschreiben der Korrekturplots (bisher Genpt)!!! Eta-Selektion bisher weiter in GenJetColEta):" <<endl;
   if(!setvalues){
   for(Int_t i=0;i<no_X_labels_;i++)  cout << i << ": " << X_labels_[i] << endl;
   cin >>bin_choice;
   fflush(stdin);
   }
   cout << "Es wurde " << bin_choice << ", also " << X_labels_[bin_choice]<< ", ausgewaehlt. Danke..." << endl;

   cout << "Bitte Pseudorapiditäts-Region auswählen:" <<endl;
   if(!setvalues){
   for(Int_t i=0;i<no_eta_region_;i++)  cout << i << ": " << eta_region_labels_[i] << endl;
   cin >>eta_choice;
   fflush(stdin);
   }
   cout << "Es wurde " << eta_choice << ", also " << eta_region_labels_[eta_choice]<< ", ausgewaehlt. Danke..." << endl;

   cout << "Bitte Exportformat für Bilder wählen ('A' für keinen Export)" <<endl;
   if(!setvalues){
   cin >>img_choice;
   fflush(stdin);
   }
   cout << "Es wurde " << img_choice << " ausgewaehlt. Danke..." << endl;

   Book_Histos();



   Long64_t nentries = fChain->GetEntriesFast();

   if(entries_to_run>0)nentries=entries_to_run;

   //nentries=100000

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

          for(Int_t genjet_i =0; genjet_i<2;genjet_i++)//über erste zwei colgenjet
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


	  if(std::abs(GenJetColEta[genjet_i])>eta_region_[eta_choice].first&&std::abs(GenJetColEta[genjet_i])<eta_region_[eta_choice].second&&genjet_jet_deltar<0.25)
// before it was like this.... &&isolated&&towers_gt_1)//demand Eta-selection, deltar (genjet,jet) <.25, isolated and more than one tower per jet
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
			 tlj_Corr_Vars_counts_all_[Corr_i].at(pt_i)->Fill(Correction_Var_[Corr_i]);
			 tlj_X_response_2D_all_[Corr_i].at(pt_i)->Fill(Correction_Var_[Corr_i],_L2L3JetResponse);
			 tlj_X_response_prof_all_[Corr_i].at(pt_i)->Fill(Correction_Var_[Corr_i],_L2L3JetResponse);
			 tlj_X_raw_response_2D_all_[Corr_i].at(pt_i)->Fill(Correction_Var_[Corr_i],_JetResponse);
			 tlj_X_raw_response_prof_all_[Corr_i].at(pt_i)->Fill(Correction_Var_[Corr_i],_JetResponse);
			 tlj_X_L2_response_2D_all_[Corr_i].at(pt_i)->Fill(Correction_Var_[Corr_i],_L2JetResponse);
			 tlj_X_L2_response_prof_all_[Corr_i].at(pt_i)->Fill(Correction_Var_[Corr_i],_L2JetResponse);
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


