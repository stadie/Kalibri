#define flex_flav_distr_cxx

#include "THelpers.h"
#include "base_corr.C"
#include "flex_flav_distr.h"







void flex_flav_distr::Book_Histos()
{
  //  std::vector < std::vector < std::vector <TH1D*> > > tlj_Corr_Vars_counts_all_pt_all_PDG_;
  //  std::vector < std::vector < TH1D* > > tlj_L2L3_Response_all_PDG_;

  for (Int_t Corr_i=0;Corr_i<no_Corr_labels_;Corr_i++)
    {
      std::vector < std::vector < TH1D* > >  tlj_Corr_Vars_counts_all_PDG_;
      for(unsigned int PDG_i=0;PDG_i<no_PDG_labels_;PDG_i++)  
	{
	  std::vector < TH1D* >  tlj_Corr_Vars_counts_all_pt_;
	  for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	    {
	      tlj_Corr_Vars_counts_all_pt_.push_back(define_pt_histo( "tlj_counts_"+(X_labels_[bin_choice]+"_")+Corr_labels_[Corr_i].first+"_"+ PDG_labels_[PDG_i]+"_", ";"+(Corr_labels_[Corr_i].second+";counts"), pt_bins_, pt_i,s_phi_bins_x,corr_var_x_edges_[Corr_i].first,corr_var_x_edges_[Corr_i].second));
	      if(titleprint==0){
		tlj_Corr_Vars_counts_all_pt_.back()->SetTitle("");
	      }
	    }
	  tlj_Corr_Vars_counts_all_PDG_.push_back(tlj_Corr_Vars_counts_all_pt_);
	}
      tlj_Corr_Vars_counts_all_pt_all_PDG_.push_back(tlj_Corr_Vars_counts_all_PDG_);
    }

     for(unsigned int PDG_i=0;PDG_i<no_PDG_labels_;PDG_i++)  
       {
	 std::vector < TH1D* >  tlj_L2L3_Response_all_pt_;
	   for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	     {
	        tlj_L2L3_Response_all_pt_.push_back(define_pt_histo( "tlj_counts_"+(X_labels_[bin_choice]+"_"+ PDG_labels_[PDG_i]+"_"), ";"+(L2L3_response_label+";counts"), pt_bins_, pt_i,100,0,2.0));
	      if(titleprint==0){
		tlj_L2L3_Response_all_pt_.back()->SetTitle("");
	      }
	     }
	 tlj_L2L3_Response_all_PDG_.push_back(tlj_L2L3_Response_all_pt_);
       }
}

void flex_flav_distr::Write_Histos()
{

  TString root_name_binning= "histos_flav_distr_" + eta_region_labels_[eta_choice] + X_labels_[bin_choice] + ".root";

   TFile *outf = new TFile(root_name_binning,"RECREATE");



    TString saveFolder = "Flav_distri";
    if(chdir(saveFolder) != 0){
      mkdir(saveFolder, S_IRWXU|S_IRWXG|S_IRWXO);
      chdir(saveFolder);
    }
      for(unsigned int PDG_i=0;PDG_i<no_PDG_labels_;PDG_i++)  
       {
	       draw_TH1D_save_PS("PS_Response_PDG",img_choice, tlj_L2L3_Response_all_PDG_[PDG_i], eta_region_labels_[eta_choice] + X_labels_[bin_choice] +tlj_L2L3_Response_all_PDG_[PDG_i].at(0)->GetName(), "nice","","y0");
      for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	{
	  tlj_L2L3_Response_all_PDG_[PDG_i].at(pt_i)->Write();
	}

       }



  for (Int_t Corr_i=0;Corr_i<no_Corr_labels_;Corr_i++)
    {
      std::vector < std::vector < TH1D* > >  tlj_Corr_Vars_counts_all_PDG_;
      for(unsigned int PDG_i=0;PDG_i<no_PDG_labels_;PDG_i++)  
	{
	  tlj_Corr_Vars_counts_all_pt_all_PDG_[Corr_i].at(PDG_i);
	       draw_TH1D_save_PS("PS_PDG_Corr_counts_"+Corr_labels_[Corr_i].first,img_choice, tlj_Corr_Vars_counts_all_pt_all_PDG_[Corr_i].at(PDG_i), eta_region_labels_[eta_choice] + X_labels_[bin_choice] +Corr_labels_[Corr_i].first + tlj_Corr_Vars_counts_all_pt_all_PDG_[Corr_i].at(PDG_i).at(0)->GetName(), "nice","","y0");
	       for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
		 {
		   tlj_Corr_Vars_counts_all_pt_all_PDG_[Corr_i].at(PDG_i).at(pt_i)->Write();
		 }

	}
    }

  ////BUILD THSTACKS....
  //  std::vector < std::vector  <THStack*>  > tlj_Corr_Vars_counts_all_pt_all_Corr_stacked_;

  //	      tlj_Corr_Vars_counts_all_pt_.push_back(define_pt_histo( "tlj_counts_"+(X_labels_[bin_choice]+"_")+Corr_labels_[Corr_i].first+"_"+ PDG_labels_[PDG_i]+"_", ";"+(Corr_labels_[Corr_i].second+";counts"), pt_bins_, pt_i,s_phi_bins_x,corr_var_x_edges_[Corr_i].first,corr_var_x_edges_[Corr_i].second));



      int nEntries =3;
      TLegend *leg = new TLegend(0.2,0.85-nEntries*0.07,0.4,0.85);
      leg->SetBorderSize(0);
      leg->SetFillColor(0);
      leg->SetTextFont(42);
      leg->SetTextSize(0.04);

      leg->AddEntry(tlj_Corr_Vars_counts_all_pt_all_PDG_[0].at(0).at(0),"Light Quarks","EPL");
      leg->AddEntry(tlj_Corr_Vars_counts_all_pt_all_PDG_[0].at(1).at(0),"Gluons","EPL");
      leg->AddEntry(tlj_Corr_Vars_counts_all_pt_all_PDG_[0].at(2).at(0),"c","EPL");
      leg->AddEntry(tlj_Corr_Vars_counts_all_pt_all_PDG_[0].at(3).at(0),"b","EPL");
      leg->AddEntry(tlj_Corr_Vars_counts_all_pt_all_PDG_[0].at(4).at(0),"not matched","EPL");

      TLegend *leg3 = new TLegend(0.2,0.85-nEntries*0.07,0.4,0.85);
      leg3->SetBorderSize(0);
      leg3->SetFillColor(0);
      leg3->SetTextFont(42);
      leg3->SetTextSize(0.04);
	 
      leg3->AddEntry(tlj_Corr_Vars_counts_all_pt_all_PDG_[0].at(0).at(0),"Light Quarks","EPL");
      leg3->AddEntry(tlj_Corr_Vars_counts_all_pt_all_PDG_[0].at(1).at(0),"Gluons","EPL");
      leg3->AddEntry(tlj_Corr_Vars_counts_all_pt_all_PDG_[0].at(2).at(0),"c","EPL");
      leg3->AddEntry(tlj_Corr_Vars_counts_all_pt_all_PDG_[0].at(3).at(0),"b","EPL");
      //      leg3->AddEntry(tlj_Corr_Vars_counts_all_pt_all_PDG_[0].at(4).at(0),"not matched","EPL");


      TLegend *leg2 = new TLegend(0.2,0.85-nEntries*0.07,0.4,0.85);
      leg2->SetBorderSize(0);
      leg2->SetFillColor(0);
      leg2->SetTextFont(42);
      leg2->SetTextSize(0.04);
	 
      leg2->AddEntry(tlj_Corr_Vars_counts_all_pt_all_PDG_[0].at(0).at(0),"Light Quarks","F");
      leg2->AddEntry(tlj_Corr_Vars_counts_all_pt_all_PDG_[0].at(1).at(0),"Gluons","F");
      leg2->AddEntry(tlj_Corr_Vars_counts_all_pt_all_PDG_[0].at(2).at(0),"c","F");
      leg2->AddEntry(tlj_Corr_Vars_counts_all_pt_all_PDG_[0].at(3).at(0),"b","F");
      //      leg->AddEntry(tlj_Corr_Vars_counts_all_pt_all_PDG_[0].at(4).at(0),"not matched","F");

      TLegend *leg4 = new TLegend(0.2,0.85-nEntries*0.07,0.4,0.85);
      leg4->SetBorderSize(0);
      leg4->SetFillColor(0);
      leg4->SetTextFont(42);
      leg4->SetTextSize(0.04);
	 
      leg4->AddEntry(tlj_Corr_Vars_counts_all_pt_all_PDG_[0].at(0).at(0),"Light Quarks","F");
      leg4->AddEntry(tlj_Corr_Vars_counts_all_pt_all_PDG_[0].at(1).at(0),"Gluons","F");
      leg4->AddEntry(tlj_Corr_Vars_counts_all_pt_all_PDG_[0].at(2).at(0),"c","F");
      leg4->AddEntry(tlj_Corr_Vars_counts_all_pt_all_PDG_[0].at(3).at(0),"b","F");
      leg4->AddEntry(tlj_Corr_Vars_counts_all_pt_all_PDG_[0].at(4).at(0),"not matched","F");




  for (Int_t Corr_i=0;Corr_i<no_Corr_labels_;Corr_i++)
    {
      std::vector  <THStack*> tlj_Corr_Vars_counts_all_pt_stacked_;
      for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	{
	  tlj_Corr_Vars_counts_all_pt_stacked_.push_back(define_pt_stack("tlj_counts_"+(X_labels_[bin_choice]+"_")+Corr_labels_[Corr_i].first+"_", ";"+(Corr_labels_[Corr_i].second+";counts"), pt_bins_, pt_i));
	  //	  if(titleprint==0){
	  //	    tlj_Corr_Vars_counts_all_pt_stacked_.back()->SetTitle("");
	  //	  }

	  //THStack *stack_test= new THStack(TString("hs_")+tlj_Corr_Vars_counts_all_pt_all_PDG_[Corr_i].at(0).at(pt_i)->GetName(),TString("hstitle_")+tlj_Corr_Vars_counts_all_pt_all_PDG_[Corr_i].at(0).at(pt_i)->GetName());
	  for(unsigned int PDG_i=0;PDG_i<no_PDG_labels_-1;PDG_i++)  
	    {
	      tlj_Corr_Vars_counts_all_pt_all_PDG_[Corr_i].at(PDG_i).at(pt_i)->SetFillColor(PDG_i+1);
	      //	      tlj_Corr_Vars_counts_all_pt_all_PDG_[Corr_i].at(PDG_i).at(pt_i)->SetFillStyle(3001+PDG_i);
       	      tlj_Corr_Vars_counts_all_pt_all_PDG_[Corr_i].at(PDG_i).at(pt_i)->SetMarkerColor(PDG_i+1);
	      tlj_Corr_Vars_counts_all_pt_all_PDG_[Corr_i].at(PDG_i).at(pt_i)->SetMarkerStyle(good_markers[PDG_i]);
	      tlj_Corr_Vars_counts_all_pt_stacked_.back()->Add(tlj_Corr_Vars_counts_all_pt_all_PDG_[Corr_i].at(PDG_i).at(pt_i));
	    }
	  //	   tlj_Corr_Vars_counts_all_pt_stacked_.push_back(stack_test);
	    //tlj_Corr_Vars_counts_all_pt_all_PDG_[Corr_i].at(PDG_i).at(pt_i)->Write();
	}

       draw_THStack_save_PS(leg,tlj_Corr_Vars_counts_all_pt_stacked_,TString("PS_")+tlj_Corr_Vars_counts_all_pt_stacked_.at(0)->GetName(), "nice", "hist_lego_nostack_legend", "x0_y0_z0", 0, -1,0,-1,img_choice, "PS_corr_var_stacks_full");
    }



  for (Int_t Corr_i=0;Corr_i<no_Corr_labels_;Corr_i++)
    {
      std::vector  <THStack*> tlj_Corr_Vars_counts_all_pt_stacked_rel_;
      for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	{
	  tlj_Corr_Vars_counts_all_pt_stacked_rel_.push_back(define_pt_stack("tlj_counts_rel_"+(X_labels_[bin_choice]+"_")+Corr_labels_[Corr_i].first+"_", ";"+(Corr_labels_[Corr_i].second+";fraction"), pt_bins_, pt_i));
	  //	  if(titleprint==0){
	  //	    tlj_Corr_Vars_counts_all_pt_stacked_.back()->SetTitle("");
	  //	  }

	  //THStack *stack_test= new THStack(TString("hs_")+tlj_Corr_Vars_counts_all_pt_all_PDG_[Corr_i].at(0).at(pt_i)->GetName(),TString("hstitle_")+tlj_Corr_Vars_counts_all_pt_all_PDG_[Corr_i].at(0).at(pt_i)->GetName());
	  for(unsigned int PDG_i=0;PDG_i<no_PDG_labels_-1;PDG_i++)  
	    {
	      tlj_Corr_Vars_counts_all_pt_all_PDG_[Corr_i].at(PDG_i).at(pt_i)->Divide(tlj_Corr_Vars_counts_all_pt_all_PDG_[Corr_i].at(PDG_i).at(pt_i),tlj_Corr_Vars_counts_all_pt_all_PDG_[Corr_i].at(no_PDG_labels_-1).at(pt_i));
	      tlj_Corr_Vars_counts_all_pt_all_PDG_[Corr_i].at(PDG_i).at(pt_i)->SetFillColor(PDG_i+1);
	      tlj_Corr_Vars_counts_all_pt_all_PDG_[Corr_i].at(PDG_i).at(pt_i)->SetMarkerColor(PDG_i+1);
	      tlj_Corr_Vars_counts_all_pt_all_PDG_[Corr_i].at(PDG_i).at(pt_i)->SetMarkerStyle(good_markers[PDG_i]);
	      tlj_Corr_Vars_counts_all_pt_stacked_rel_.back()->Add(tlj_Corr_Vars_counts_all_pt_all_PDG_[Corr_i].at(PDG_i).at(pt_i));
	    }
	  //	   tlj_Corr_Vars_counts_all_pt_stacked_.push_back(stack_test);
	    //tlj_Corr_Vars_counts_all_pt_all_PDG_[Corr_i].at(PDG_i).at(pt_i)->Write();
	}


     draw_THStack_save_PS(leg4,tlj_Corr_Vars_counts_all_pt_stacked_rel_,TString("PS_")+tlj_Corr_Vars_counts_all_pt_stacked_rel_.at(0)->GetName(), "nice", "hist_lego_nostack_legend_norm", "x0_y0_z0", 0, -1,-0.05,1.05,img_choice, "PS_corr_var_stacks_full");
    }




  std::vector  <THStack*> tlj_L2L3_Response_flav_;
      for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	{
	  tlj_L2L3_Response_flav_.push_back(define_pt_stack("tlj_L2L3_response_"+(X_labels_[bin_choice]+"_"), ";"+(L2L3_response_label+";counts"), pt_bins_, pt_i));
	  for(unsigned int PDG_i=0;PDG_i<no_PDG_labels_-1;PDG_i++)  
	    {	
	      //	      tlj_L2L3_Response_all_PDG_[PDG_i].at(pt_i)->Divide(tlj_L2L3_Response_all_PDG_[PDG_i].at(pt_i),tlj_L2L3_Response_all_PDG_[no_PDG_labels_-1].at(pt_i));
	      tlj_L2L3_Response_all_PDG_[PDG_i].at(pt_i)->SetFillColor(PDG_i+1);
	      tlj_L2L3_Response_all_PDG_[PDG_i].at(pt_i)->SetMarkerColor(PDG_i+1);
	      tlj_L2L3_Response_all_PDG_[PDG_i].at(pt_i)->SetMarkerStyle(good_markers[PDG_i]);
	      tlj_L2L3_Response_flav_.back()->Add(tlj_L2L3_Response_all_PDG_[PDG_i].at(pt_i));

	   }
	}
      draw_THStack_save_PS(leg,tlj_L2L3_Response_flav_,TString("PS_")+tlj_L2L3_Response_all_PDG_[0].at(0)->GetName(), "nice", "hist_lego_nostack_legend_norm", "x0_y0_z0", 0, -1,0,-1,img_choice, "PS_response_stacks_full");

  std::vector  <THStack*> new_tlj_L2L3_Response_flav_;
      for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	{
	  new_tlj_L2L3_Response_flav_.push_back(define_pt_stack("new_tlj_L2L3_response_"+(X_labels_[bin_choice]+"_"), ";"+(L2L3_response_label+";counts"), pt_bins_, pt_i));
	  for(unsigned int PDG_i=0;PDG_i<no_PDG_labels_-2;PDG_i++)  
	    {	
	      //	      tlj_L2L3_Response_all_PDG_[PDG_i].at(pt_i)->Divide(tlj_L2L3_Response_all_PDG_[PDG_i].at(pt_i),tlj_L2L3_Response_all_PDG_[no_PDG_labels_-1].at(pt_i));
	      tlj_L2L3_Response_all_PDG_[PDG_i].at(pt_i)->SetFillColor(PDG_i+1);
	      tlj_L2L3_Response_all_PDG_[PDG_i].at(pt_i)->SetMarkerColor(PDG_i+1);
	      tlj_L2L3_Response_all_PDG_[PDG_i].at(pt_i)->SetMarkerStyle(good_markers[PDG_i]);
	      new_tlj_L2L3_Response_flav_.back()->Add(tlj_L2L3_Response_all_PDG_[PDG_i].at(pt_i));

	   }
	}
      draw_THStack_save_PS(leg3,new_tlj_L2L3_Response_flav_,TString("PS_few_")+tlj_L2L3_Response_all_PDG_[0].at(0)->GetName(), "nice", "hist_lego_nostack_legend_norm", "x0_y0_z0", 0.5, 1.5,0,-1,img_choice, "PS_response_stacks_full");



  std::vector  <THStack*> tlj_L2L3_Response_flav_rel_;
      for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	{
	  tlj_L2L3_Response_flav_rel_.push_back(define_pt_stack("tlj_L2L3_response_rel_"+(X_labels_[bin_choice]+"_"), ";"+(L2L3_response_label+";counts"), pt_bins_, pt_i));
	  for(unsigned int PDG_i=0;PDG_i<no_PDG_labels_-1;PDG_i++)  
	    {	
       	      tlj_L2L3_Response_all_PDG_[PDG_i].at(pt_i)->Divide(tlj_L2L3_Response_all_PDG_[PDG_i].at(pt_i),tlj_L2L3_Response_all_PDG_[no_PDG_labels_-1].at(pt_i));
	      tlj_L2L3_Response_all_PDG_[PDG_i].at(pt_i)->SetFillColor(PDG_i+1);
	      tlj_L2L3_Response_all_PDG_[PDG_i].at(pt_i)->SetMarkerColor(PDG_i+1);
	      tlj_L2L3_Response_all_PDG_[PDG_i].at(pt_i)->SetMarkerStyle(good_markers[PDG_i]);
	      tlj_L2L3_Response_flav_rel_.back()->Add(tlj_L2L3_Response_all_PDG_[PDG_i].at(pt_i));

	   }
	}
      draw_THStack_save_PS(leg,tlj_L2L3_Response_flav_rel_,TString("PS_rel_")+tlj_L2L3_Response_all_PDG_[0].at(0)->GetName(), "nice", "hist_lego_nostack_legend_norm", "x0_y0_z0", 0, -1,-0.05,1.05,img_choice, "PS_response_stacks_full");

    

  

  chdir("..");



   outf->Close();


}


void flex_flav_distr::Loop()
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

   
   cout << "Bitte Binningvariable wählen (ACHTUNG: I.A. sollte die Binningvariable hier die gleiche sein wie beim Rausschreiben der Korrekturplots (bisher Genpt)!!! Eta-Selektion bisher weiter in GenJetColEta):" <<endl;
   for(Int_t i=0;i<no_X_labels_;i++)  cout << i << ": " << X_labels_[i] << endl;
   cin >>bin_choice;
   fflush(stdin);
   cout << "Es wurde " << bin_choice << ", also " << X_labels_[bin_choice]<< ", ausgewaehlt. Danke..." << endl;

   cout << "Bitte Pseudorapiditäts-Region auswählen:" <<endl;
   for(Int_t i=0;i<no_eta_region_;i++)  cout << i << ": " << eta_region_labels_[i] << endl;
   cin >>eta_choice;
   fflush(stdin);
   cout << "Es wurde " << eta_choice << ", also " << eta_region_labels_[eta_choice]<< ", ausgewaehlt. Danke..." << endl;

   cout << "Bitte Exportformat für Bilder wählen ('A' für keinen Export)" <<endl;
   cin >>img_choice;
   fflush(stdin);
   cout << "Es wurde " << img_choice << " ausgewaehlt. Danke..." << endl;

   Book_Histos();



   Long64_t nentries = fChain->GetEntriesFast();
   //      nentries =100000;

   if(entries_to_run>0)nentries=entries_to_run;

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
	      Int_t PDG_selector;
	      Int_t PDG_id=GenPartId_algo[match];
	      //	      cout << PDG_id << endl;
	      if(std::abs(PDG_id)==3||std::abs(PDG_id)==2||std::abs(PDG_id)==1)PDG_selector=0;
	      else if(std::abs(PDG_id)==21)PDG_selector=1;
	      else if(std::abs(PDG_id)==4)PDG_selector=2;
	      else if(std::abs(PDG_id)==5)PDG_selector=3;
	      else {PDG_selector=4;}


	      //Get variables needed for filling...
       	      base_corr::Fill_obvious_vars(genjet_i, match);
	      std::vector < Double_t > X_Var_=get_X_Vars_(genjet_i,match);
	      std::vector < Double_t > Correction_Var_=get_Correction_Vars_(genjet_i,match);


	      //Fill histos now...
	     for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	       {//Jeder der zwei l. Jet wird dem pt-bin zugeordnet, der zu ihm passt
		 if(X_Var_[bin_choice]>pt_bins_[pt_i].first&&X_Var_[bin_choice]<pt_bins_[pt_i].second)
		   {
		     //		     cout << "test " << endl;
		     //HIER HIN
		     //HIER HIN
		     //HIER HIN
		     //HIER HIN

		     tlj_L2L3_Response_all_PDG_[PDG_selector].at(pt_i)->Fill(_L2L3JetResponse);
		     tlj_L2L3_Response_all_PDG_[no_PDG_labels_-1].at(pt_i)->Fill(_L2L3JetResponse);

		     //		     cout << "test2 " << endl;

		     for (Int_t Corr_i=0;Corr_i<no_Corr_labels_;Corr_i++)
		       {
			 tlj_Corr_Vars_counts_all_pt_all_PDG_[Corr_i].at(PDG_selector).at(pt_i)->Fill(Correction_Var_[Corr_i]);
			 tlj_Corr_Vars_counts_all_pt_all_PDG_[Corr_i].at(no_PDG_labels_-1).at(pt_i)->Fill(Correction_Var_[Corr_i]);
		       }

		     //HIER HIN
		     //HIER HIN
		     //HIER HIN
		     //HIER HIN
		     //HIER HIN


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


