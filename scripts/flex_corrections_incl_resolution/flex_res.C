#define flex_res_cxx
#include "THelpers.h"
#include "base_corr.C"
#include "flex_res.h"

void flex_res::Import_Histos()
{

  TString root_name_binning= "histos_" + X_labels_[bin_choice] + ".root";

  cout << "DEBUG: " << root_name_binning << endl;

  TFile* _file=new TFile(root_name_binning,"OPEN");

  //    TFile f("file.root");
    if (_file->IsZombie()) {
       cout << "Error opening file" << endl;
       exit(-1);
    }
    else cout << "ROOT-Datei erfolgreich geladen. " << endl;



  //  std::vector < std::vector < TH1D* > > tlj_X_response_all_corrections_all_;


  for (Int_t X_i=0;X_i<no_X_labels_;X_i++)
    {
      std::vector < TH1D* >  tlj_X_counts_;
       for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	 {
	   //	   cout << define_pt_histo_name("tlj_counts_"+(X_labels_[X_i]+"_"), pt_bins_, pt_i).Data() << endl;
      	   tlj_X_counts_.push_back((TH1D*)_file->Get(define_pt_histo_name("tlj_counts_"+(X_labels_[X_i]+"_"), pt_bins_, pt_i).Data()));
	 }
       tlj_X_counts_all_.push_back(tlj_X_counts_);
    }

  for (Int_t Corr_i=0;Corr_i<no_Corr_labels_;Corr_i++)
    {
      std::vector < TH2D* >  tlj_X_response_2D_;
      std::vector < TProfile* >  tlj_X_response_prof_;
      std::vector < TH1D* >  tlj_X_response_GMP_mean_;
       for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	 {
	   tlj_X_response_2D_.push_back((TH2D*)_file->Get(define_pt_histo_name( "tlj_X_response_2D_"+Corr_labels_[Corr_i].first+"_", pt_bins_, pt_i)));
	   tlj_X_response_prof_.push_back((TProfile*)_file->Get(define_pt_histo_name( "tlj_X_response_prof_"+Corr_labels_[Corr_i].first+"_", pt_bins_, pt_i)));
	   TString temp =  (TString) tlj_X_response_2D_.back()->GetName() + "_1";
	   tlj_X_response_GMP_mean_.push_back((TProfile*)_file->Get(temp.Data()));

	 }

              tlj_X_response_2D_all_.push_back(tlj_X_response_2D_);
              tlj_X_response_prof_all_.push_back(tlj_X_response_prof_);
              tlj_X_response_GMP_mean_all_.push_back(tlj_X_response_GMP_mean_);

    }

  for (Int_t Corr_i=0;Corr_i<no_all_Corr_labels_;Corr_i++)
    {
      std::vector < TH1D* > tlj_X_response_all_corrections_;
       for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	 {
	   if(Corr_i<no_Corr_labels_) tlj_X_response_all_corrections_.push_back((TH1D*)tlj_X_response_prof_all_[Corr_i].at(pt_i));
	   else tlj_X_response_all_corrections_.push_back(tlj_X_response_GMP_mean_all_[Corr_i-no_Corr_labels_].at(pt_i));
	   //	   cout << "ACHTUNG!!!!!!!!!!!!!!" << tlj_X_response_all_corrections_.back()->GetName() << endl;

	 }
       tlj_X_response_all_corrections_all_.push_back(tlj_X_response_all_corrections_);
    }


  //  tlj_X_response_all_corrections_all_[9].at(0)->Draw();

}


void flex_res::Book_Histos()
{

  for (Int_t Corr_i=0;Corr_i<no_all_Corr_labels_;Corr_i++)
    {
      std::vector < TH1D* >  tlj_corrected_response_barrel_;
       for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	 {
	   tlj_corrected_response_barrel_.push_back(define_pt_histo( "tlj_corrected_response_barrel_"+Corr_labels_[Corr_i].first+"_", ";corrected Response", pt_bins_, pt_i, response_bins, response_low, response_high));
	   //	   cout << tlj_corrected_response_barrel_.back()->GetName() << endl;
	   //	   cout << "Corr_i: " <<Corr_i <<" pt_i: " <<pt_i << endl;
	 }
       tlj_corrected_response_barrel_all_.push_back(tlj_corrected_response_barrel_);
    }


 for(unsigned int Sel_Corr_i=0;Sel_Corr_i<Corr_selected_labels_.size();Sel_Corr_i++)
   {
  std::vector < std::vector < TH2D* > > tlj_Sel_Correlations_2D_one_corr_;
   for(unsigned int Sel_Corr_j=0;Sel_Corr_j<Corr_selected_labels_.size();Sel_Corr_j++)
     {
      std::vector < TH2D* >  tlj_Sel_Correlations_2D_;
       for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	 {
	   tlj_Sel_Correlations_2D_.push_back(define_pt_histo( "tlj_Sel_Correlations_2D_"+
	   Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first+"_vs_" + Corr_labels_[Corr_selected_labels_[Sel_Corr_j]].first + "_", ";"+Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first+";"+  Corr_labels_[Corr_selected_labels_[Sel_Corr_j]].first, pt_bins_, pt_i, s_phi_bins_x,s_phi_xlow,s_phi_xhigh, s_phi_bins_x,s_phi_xlow,s_phi_xhigh));
	 }
       tlj_Sel_Correlations_2D_one_corr_.push_back(tlj_Sel_Correlations_2D_);
     }
   tlj_Sel_Correlations_2D_all_.push_back(tlj_Sel_Correlations_2D_one_corr_);
   }


}

void flex_res::Write_Histos()
{


 root_resol_name_binning= "histos_resolution_" + X_labels_[bin_choice] + ".root";

   TFile *outf = new TFile(root_resol_name_binning,"RECREATE");

     

  for (Int_t Corr_i=0;Corr_i<no_all_Corr_labels_;Corr_i++)
    {
      draw_TH1D_save_with_Gauss_Fit_PS( tlj_corrected_response_barrel_all_[Corr_i],tlj_corrected_response_barrel_all_[Corr_i].at(0)->GetName(),"nice");
      for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	{
	  tlj_corrected_response_barrel_all_[Corr_i].at(pt_i)->Write();
	  //	  tlj_X_response_all_ tlj_X_response_all_corrections_.back()->Draw();
	}
      
    }
   

  for (Int_t Corr_i=0;Corr_i<no_all_Corr_labels_;Corr_i++)
    {
      draw_TH1D_save_PS(tlj_X_response_all_corrections_all_[Corr_i],"tlj_X_response_all_corrections_all_"+Corr_labels_[Corr_i].first,"nice");  
      for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	{
	  tlj_X_response_all_corrections_all_[Corr_i].at(pt_i)->Write();
	  //	  tlj_X_response_all_ tlj_X_response_all_corrections_.back()->Draw();
	}


    }

  for(unsigned int Sel_Corr_i=0;Sel_Corr_i<Corr_selected_labels_.size();Sel_Corr_i++)
    {
      //			 std::vector < std::vector < TH2D* > > tlj_Sel_Correlations_2D_one_corr_;
      for(unsigned int Sel_Corr_j=0;Sel_Corr_j<Corr_selected_labels_.size();Sel_Corr_j++)
	{
	  //			     std::vector < TH2D* >  tlj_Sel_Correlations_2D_;
	  
	        draw_TH2D_save_PS(true,img_choice,tlj_Sel_Correlations_2D_all_[Sel_Corr_i].at(Sel_Corr_j)," tlj_Sel_Correlations_"+Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first + "_vs_" + Corr_labels_[Corr_selected_labels_[Sel_Corr_j]].first ,"nice","colz");
	  for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	    {
	      tlj_Sel_Correlations_2D_all_[Sel_Corr_i].at(Sel_Corr_j).at(pt_i)->Write();
	    }
	}
    }


   outf->Close();


}

void flex_res::Write_TGraphErrors()
{

   TFile *outf = new TFile(root_resol_name_binning,"UPDATE");

      draw_TGraphErrors_save_PS(img_choice,tlj_Response_Graphs_mean_all_,tlj_Response_Graphs_mean_all_.at(0)->GetName(),"nice","ALP","no","x1");
      draw_TGraphErrors_save_PS(img_choice,tlj_Response_Graphs_mean_all_,tlj_Response_Graphs_mean_all_.at(0)->GetName(),"nice");
      draw_TGraphErrors_save_PS(img_choice,tlj_Response_Graphs_rel_sigma_all_,tlj_Response_Graphs_rel_sigma_all_.at(0)->GetName(),"nice","ALP", "yes","x1",0,-1,0,0.5);
      draw_TGraphErrors_save_PS(img_choice,tlj_Response_Graphs_rel_sigma_all_,tlj_Response_Graphs_rel_sigma_all_.at(0)->GetName(),"nice","ALP", "yes","x0",0,-1,0,0.5);
      draw_TGraphErrors_save_PS(img_choice,tlj_Response_Graphs_rel_resol_all_,tlj_Response_Graphs_rel_resol_all_.at(0)->GetName(),"nice","ALP", "no","x1");
      draw_TGraphErrors_save_PS(img_choice,tlj_Response_Graphs_rel_resol_all_,tlj_Response_Graphs_rel_resol_all_.at(0)->GetName(),"nice","ALP", "no","x0");

  for (Int_t Corr_i=0;Corr_i<no_all_Corr_labels_;Corr_i++)
    {
      tlj_Response_Graphs_mean_all_.at(Corr_i)->Write();
      tlj_Response_Graphs_rel_sigma_all_.at(Corr_i)->Write();
      tlj_Response_Graphs_rel_resol_all_.at(Corr_i)->Write();
    }


  if(Corr_selected_labels_.size()>0)
{
      draw_TGraphErrors_save_PS(img_choice,tlj_Response_Graphs_mean_selec_,tlj_Response_Graphs_mean_selec_.at(0)->GetName(),"nice","ALP","no","x1");
      draw_TGraphErrors_save_PS(img_choice,tlj_Response_Graphs_mean_selec_,tlj_Response_Graphs_mean_selec_.at(0)->GetName(),"nice","ALP","no","x0");
      draw_TGraphErrors_save_PS(img_choice,tlj_Response_Graphs_rel_sigma_selec_,tlj_Response_Graphs_rel_sigma_selec_.at(0)->GetName(),"nice","ALP","yes","x0",0,-1,0,0.5);
      draw_TGraphErrors_save_PS(img_choice,tlj_Response_Graphs_rel_sigma_selec_,tlj_Response_Graphs_rel_sigma_selec_.at(0)->GetName(),"nice","ALP","yes","x1",0,-1,0,0.5);
      draw_TGraphErrors_save_PS(img_choice,tlj_Response_Graphs_rel_resol_selec_,tlj_Response_Graphs_rel_resol_selec_.at(0)->GetName(),"nice","ALP","no","x0");
      draw_TGraphErrors_save_PS(img_choice,tlj_Response_Graphs_rel_resol_selec_,tlj_Response_Graphs_rel_resol_selec_.at(0)->GetName(),"nice","ALP","no","x1");
}


   for(unsigned int Sel_Corr_i=0;Sel_Corr_i<Corr_selected_labels_.size();Sel_Corr_i++)
     {
       tlj_Response_Graphs_mean_selec_.at(Sel_Corr_i)->Write();
       tlj_Response_Graphs_rel_sigma_selec_.at(Sel_Corr_i)->Write();
       tlj_Response_Graphs_rel_resol_selec_.at(Sel_Corr_i)->Write();
     }



 TLegend *leg_all;
 leg_all = new TLegend(0.1,0.1,0.48,0.3);
 //   leg->SetHeader("Legende");
 for (Int_t Corr_i=0;Corr_i<no_all_Corr_labels_;Corr_i++)
   {
   leg_all->AddEntry(tlj_Response_Graphs_mean_all_[Corr_i],Corr_labels_[Corr_i].first,"lep");
   }

 TLegend *leg_selected;
 leg_selected = new TLegend(0.1,0.1,0.48,0.3);
 //   leg->SetHeader("Legende");
 for(unsigned int Sel_Corr_i=0;Sel_Corr_i<Corr_selected_labels_.size();Sel_Corr_i++)
   {
   leg_selected->AddEntry(tlj_Response_Graphs_mean_selec_[Sel_Corr_i],Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first,"lep");
   }



 draw_graphs( tlj_Response_Graphs_mean_all_, 0., 1.3, leg_all, "RESPONSE_res_All_response_graph_mean_w_mean_e");
 draw_graphs( tlj_Response_Graphs_rel_sigma_all_, -0.1, 0.5, leg_all, "RESPONSE_res_All_response_graph_rel_sigma");
 draw_graphs( tlj_Response_Graphs_rel_resol_all_, 0.6, 1.1, leg_all, "RESPONSE_res_All_response_graph_rel_sigma_relative_to_L2L3");


 TString selection;
   for(unsigned int Sel_Corr_i=0;Sel_Corr_i<Corr_selected_labels_.size();Sel_Corr_i++)
     {
       selection+=Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first +"-";
     }
   cout << selection << endl;


  if(Corr_selected_labels_.size()>0)
{
 draw_graphs( tlj_Response_Graphs_mean_selec_, 0.85, 1.15, leg_selected, "RESPONSE_res_selected_" + selection + "_response_graph_mean");
 draw_graphs( tlj_Response_Graphs_rel_sigma_selec_, -0.1, 0.5, leg_selected, "RESPONSE_res_selected_" + selection + "_response_graph_rel_sigma");
 draw_graphs( tlj_Response_Graphs_rel_resol_selec_, 0.6, 1.1, leg_selected, "RESPONSE_res_selected_" + selection + "_response_graph_rel_sigma_relative_to_L2L3");
}



   outf->Close();

}

void flex_res::draw_graphs(std::vector <TGraphErrors*> graphs_, Double_t ylow, Double_t yhigh, TLegend *legend, TString PDF_PNG_name)
{

  if(graphs_.size()>=1)
    {
graphs_[0]->GetYaxis()->SetRangeUser(ylow,yhigh);
graphs_[0]->Draw("ALP");
    }

 for(unsigned int a_i=0; a_i<graphs_.size();a_i++)
   {
      graphs_[a_i]->SetLineColor(a_i+1);
      graphs_[a_i]->Draw("same"); 
   }
 legend->Draw();

 test->Print(PDF_PNG_name+".pdf");
 test->Print(PDF_PNG_name+".png");

  test->SetLogx();

 test->Print(PDF_PNG_name+"_logX_" +".pdf");
 test->Print(PDF_PNG_name+"_logX_" +".png");

  test->SetLogx(0);


}


TGraphErrors* flex_res::make_graph(std::vector < std::vector < Double_t > >  Double_gauss_, Int_t y_para, TString title, Int_t X_par)
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

  TGraphErrors* temp = new TGraphErrors(size, X_array, Y_array, X_e_array, Y_e_array);
  temp->SetTitle (title+";"+X_labels_[X_par]+";"+double_gauss_labels_[y_para]);
  temp->SetName(title);

  return temp;
}


void flex_res::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L flex_res.C
//      Root > flex_res t
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

   cout << "Bitte Binningvariable wählen (ACHTUNG: I.A. sollte die Binningvariable hier die gleiche sein wie beim Rausschreiben der Korrekturplots (bisher Genpt)!!!):" <<endl;
   for(Int_t i=0;i<no_X_labels_;i++)  cout << i << ": " << X_labels_[i] << endl;
   cin >>bin_choice;
   fflush(stdin);
   cout << "Es wurde " << bin_choice << ", also " << X_labels_[bin_choice]<< ", ausgewaehlt. Danke..." << endl;

   

   cout << "Bitte Parametrisierung wählen:" <<endl;
   for(Int_t i=0;i<no_X_labels_;i++)  cout << i << ": " << X_labels_[i] << endl;
   cin >>X_choice;
   fflush(stdin);
   cout << "Es wurde " << X_choice << ", also " << X_labels_[X_choice]<< ", ausgewaehlt. Danke..." << endl;



   cout << "Bitte Korrekturen zum Vergleich auswaehlen:" <<endl;
   //   cout << "Test: " << (char)66 << endl;

   //   cout << Corr_labels_[2].first << endl;

   for(Int_t i=0;i<no_all_Corr_labels_;i++)  
     {
       cout << (char)(i+65) << ": " << Corr_labels_[i].first << endl;
     }
   cin >> Corr_choice;
   fflush(stdin);
   cout << "Es wurde " << Corr_choice << " ausgewaehlt. Danke... Gewaehlt wurde also: " << endl;


   for(Int_t i=0;i<no_all_Corr_labels_;i++)  
     {
       if(Corr_choice.Contains((char)(i+65)))
	 {
	   cout << (char)(i+65) << ": " << Corr_labels_[i].first << endl;
	   Corr_selected_labels_.push_back(i);
	 }
     }
   if(Corr_selected_labels_.size()==0)cout << "Es wurde nichts ausgewaehlt!!!" << endl;

   cout << "Bitte Exportformat für Bilder wählen ('A' für keinen Export)" <<endl;
   cin >>img_choice;
   fflush(stdin);
   cout << "Es wurde " << img_choice << " ausgewaehlt. Danke..." << endl;

   Import_Histos();
   Book_Histos();





   
   Long64_t nentries = fChain->GetEntriesFast();


   //          nentries=20000;

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


	  if(std::abs(GenJetColEta[genjet_i])<1.1&&genjet_jet_deltar<0.25&&isolated&&towers_gt_1)//demand Eta-selection, deltar (genjet,jet) <.25, isolated and more than one tower per jet
	    {
	      Int_t match=GenJetColJetIdx[genjet_i];


	      //Get variables needed for filling...
       	      base_corr::Fill_obvious_vars(genjet_i, match);
	      std::vector < Double_t > X_Var_=get_X_Vars_(genjet_i,match);
	      std::vector < Double_t > Correction_Var_=get_Correction_Vars_(genjet_i,match);

	     for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	       {//Jeder der zwei l. Jet wird dem pt-bin zugeordnet, der zu ihm passt
		 if(X_Var_[bin_choice]>pt_bins_[pt_i].first&&X_Var_[bin_choice]<pt_bins_[pt_i].second)
		   {
		     std::vector < Double_t > Corrected_Response_;
		     Corrected_Response_.push_back(_JetResponse);				
		     Corrected_Response_.push_back(_L2L3JetResponse);
		     for(Int_t Corr_i=2;Corr_i<no_Corr_labels_;Corr_i++)
		       {
			 Corrected_Response_.push_back(_L2L3JetResponse/tlj_X_response_prof_all_[Corr_i].at(pt_i)->GetBinContent(tlj_X_response_prof_all_[Corr_i].at(pt_i)->GetXaxis()->FindBin(Correction_Var_[Corr_i])));				
		       }
		     Corrected_Response_.push_back(_JetResponse);	
		     Corrected_Response_.push_back(_L2L3JetResponse);
		     for(Int_t Corr_i=2+no_Corr_labels_;Corr_i<no_all_Corr_labels_;Corr_i++)
		       {
			 Corrected_Response_.push_back(_L2L3JetResponse/tlj_X_response_all_corrections_all_[Corr_i].at(pt_i)->GetBinContent(tlj_X_response_all_corrections_all_[Corr_i].at(pt_i)->GetXaxis()->FindBin(Correction_Var_[Corr_i-no_Corr_labels_])));				
		       }



		     //  cout << "test" << endl;
		     for(Int_t Corr_i=0;Corr_i<no_all_Corr_labels_;Corr_i++)
		       {
			 tlj_corrected_response_barrel_all_[Corr_i].at(pt_i)->Fill(Corrected_Response_[Corr_i]);
			 //	   tlj_X_response_2D_all_[Corr_i].at(pt_i)->Fill(Correction_Var_[Corr_i],L2L3JetResponse_.back());
			 //			 tlj_X_response_prof_all_[Corr_i].at(pt_i)->Fill(Correction_Var_[Corr_i],L2L3JetResponse_.back());
			 //  cout << "test" << endl;
		       }


		     for(unsigned int Sel_Corr_i=0;Sel_Corr_i<Corr_selected_labels_.size();Sel_Corr_i++)
		       {
			 //			 std::vector < std::vector < TH2D* > > tlj_Sel_Correlations_2D_one_corr_;
			 for(unsigned int Sel_Corr_j=0;Sel_Corr_j<Corr_selected_labels_.size();Sel_Corr_j++)
			   {
			     //			     std::vector < TH2D* >  tlj_Sel_Correlations_2D_;

			     tlj_Sel_Correlations_2D_all_[Sel_Corr_i].at(Sel_Corr_j).at(pt_i)->Fill(Correction_Var_[Corr_selected_labels_[Sel_Corr_i]],Correction_Var_[Corr_selected_labels_[Sel_Corr_j]]);
			   }
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





   std::vector < std::vector < std::vector < Double_t > > > Double_gauss_all_;
   std::vector < std::vector < std::vector < Double_t > > > Rel_resol_all_;

   for(Int_t Corr_i=0;Corr_i<no_all_Corr_labels_;Corr_i++)
     {//Jeder der zwei l. Jet wird dem pt-bin zugeordnet, der zu ihm passt

       std::vector < std::vector < Double_t > > Double_gauss_;
       std::vector < std::vector < Double_t > > Rel_resol_;
       for(Int_t pt_i=0;pt_i<no_pt_bins_;pt_i++)
	 {
	   cout << "Anzahl der Einträge: " << tlj_corrected_response_barrel_all_[Corr_i].at(pt_i)->GetEntries() << " Corr_Nummer: " << Corr_i << "  Pt_i " << pt_i << endl;

	   //	   tlj_corrected_response_barrel_all_[Corr_i].at(pt_i)->Draw();
	   if(tlj_corrected_response_barrel_all_[Corr_i].at(pt_i)->GetEntries()>100)
	     {
	       cout << "test1" << endl;
	   cout << "Anzahl der Einträge: " << tlj_corrected_response_barrel_all_[Corr_i].at(pt_i)->GetEntries() << " Corr_Nummer: " << Corr_i << "  Pt_i " << pt_i << endl;
	       Double_gauss_.push_back(doublefit_gaus(tlj_corrected_response_barrel_all_[Corr_i].at(pt_i)));
	       std::vector < Double_t > rel_resol_pair_;
	       cout << "test2" << endl;
	       rel_resol_pair_.push_back(doublefit_gaus(tlj_corrected_response_barrel_all_[Corr_i].at(pt_i)).at(2)/doublefit_gaus(tlj_corrected_response_barrel_all_[1].at(pt_i)).at(2));
	       cout << "test3" << endl;
	       rel_resol_pair_.push_back(doublefit_gaus(tlj_corrected_response_barrel_all_[Corr_i].at(pt_i)).at(3)/doublefit_gaus(tlj_corrected_response_barrel_all_[1].at(pt_i)).at(2));
	       cout << "test4" << endl;

	       Rel_resol_.push_back(rel_resol_pair_);
	     }
	   else break;
	 }
       Double_gauss_all_.push_back(Double_gauss_);
       Rel_resol_all_.push_back(Rel_resol_);
     }



   for(Int_t Corr_i=0;Corr_i<no_all_Corr_labels_;Corr_i++)
     {
       //Mean, default X
       tlj_Response_Graphs_mean_all_.push_back(make_graph(Double_gauss_all_[Corr_i], 0, X_labels_[X_choice] +"vs.Mean_" + Corr_labels_[Corr_i].first, X_choice));
       //Rel. Sigma, default X
       tlj_Response_Graphs_rel_sigma_all_.push_back(make_graph(Double_gauss_all_[Corr_i], 2, X_labels_[X_choice] +"vs.#sigma_{Rel}_" + Corr_labels_[Corr_i].first, X_choice));
       //Rel. Sigma relative to L2L3, default X
       tlj_Response_Graphs_rel_resol_all_.push_back(make_graph(Rel_resol_all_[Corr_i], 0, X_labels_[X_choice] +"vs.#sigma_{Rel}_over_#sigma_{Rel,L2L3}_" + Corr_labels_[Corr_i].first, X_choice));
       if(Corr_i>=no_Corr_labels_)
	 {
	   tlj_Response_Graphs_mean_all_.back()->SetLineStyle(2);
	   tlj_Response_Graphs_rel_sigma_all_.back()->SetLineStyle(2);
	   tlj_Response_Graphs_rel_resol_all_.back()->SetLineStyle(2);
	 }
     }
   
   for(unsigned int Sel_Corr_i=0;Sel_Corr_i<Corr_selected_labels_.size();Sel_Corr_i++)
     {
       //Personal selection Mean, default X
       tlj_Response_Graphs_mean_selec_.push_back(make_graph(Double_gauss_all_[Corr_selected_labels_[Sel_Corr_i]], 0, "Selected_" + X_labels_[X_choice] +"vs.Mean_" + Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first, X_choice));
       //Rel. Sigma, default X
       tlj_Response_Graphs_rel_sigma_selec_.push_back(make_graph(Double_gauss_all_[Corr_selected_labels_[Sel_Corr_i]], 2, "Selected_" + X_labels_[X_choice] +"vs.#sigma_{Rel}_" + Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first, X_choice));
       //Rel. Sigma relative to L2L3, default X
       tlj_Response_Graphs_rel_resol_selec_.push_back(make_graph(Rel_resol_all_[Corr_selected_labels_[Sel_Corr_i]], 0, "Selected" + X_labels_[X_choice] +"vs.#sigma_{Rel}_over_#sigma_{Rel,L2L3}_" + Corr_labels_[Corr_selected_labels_[Sel_Corr_i]].first, X_choice));

       if(Corr_selected_labels_[Sel_Corr_i]>=no_Corr_labels_)
	 {
	   tlj_Response_Graphs_mean_selec_.back()->SetLineStyle(2);
	   tlj_Response_Graphs_rel_sigma_selec_.back()->SetLineStyle(2);
	   tlj_Response_Graphs_rel_resol_selec_.back()->SetLineStyle(2);
	 }

     }
   



   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!!
   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!!
   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!!
   /////////////////////////////////HISTO BOOKING!!!!!!!!!!!!!!!!!!!!!!!!
         Write_Histos();
         Write_TGraphErrors();







}
