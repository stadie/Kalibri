#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TChain.h"
#include "cmath" 
#include "TPad.h"
#include "TStyle.h"
#include "TSpectrum.h"
#include "TString.h"
#include "TVector2.h"
#include "TFile.h"
#include "TVirtualFitter.h"
#include "TLorentzVector.h"
void Summarizer() //type 0: ptave; 1: leadpt; 2: barrelpt
{
  bool verbose=false; //activate for debug mode
  if(verbose)std::cout<<"start the summarizer"<<endl;
  //some general settings
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(11111111);
  gStyle->SetPalette(1);
  //define alpha bins  
  if(verbose)std::cout<<"setting binning"<<endl;
  int alphaBin[4]={10, 15, 20, 30}; //special JEC variant
  //int alphaBin[4]={10, 20, 30, 40}; //normal dijet variant  
  int etaBin[8]={0, 8, 13, 19, 25, 30, 32, 52};
  //loading the files  
  if(verbose)std::cout<<"defining files"<<endl;
  TFile *f;				      
  TFile *save;
  for(int choice=0; choice<3; choice++)
    {
      if(verbose)std::cout<<"load choice="<<choice<<endl;
      for(int type=0; type<2; type++)
	{
	  if(verbose)std::cout<<"load type="<<type<<endl;  
	  //  if(type==0)f=TFile::Open("../20122012AB_195396_CORR2012SQLV7_AK5_MC_Su12Z2Star_PUS6S7_kostas_SQLV7AlternativeExtrapolAK5/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_JEC_Mikko/plots/KalibriPlots.root");
	  //  else if(type==1)f=TFile::Open("../20122012AB_195396_CORR2012SQLV7_AK5_MC_Su12Z2Star_PUS6S7_kostas_SQLV7AlternativeExtrapolAK5/dijetsFall10_TuneZ2_AK5PFCHS_weighted_residuals_JEC_Mikko/plots/KalibriPlots.root");
	  //if(type==0)f=TFile::Open("../20122012AB_196531_CORR2012SQLV7_AK5_MC_Su12Z2Star_PUS6S7_kostas_NewTry/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_JEC_Mikko/plots/KalibriPlots.root");
	  //else if(type==1)f=TFile::Open("../20122012AB_196531_CORR2012SQLV7_AK5_MC_Su12Z2Star_PUS6S7_kostas_NewTry/dijetsFall10_TuneZ2_AK5PFCHS_weighted_residuals_JEC_Mikko/plots/KalibriPlots.root");
	  //if(type==0)f=TFile::Open("../20122012ABC_199429_CORR2012SQLV7_AK5_MC_Su12Z253_kostas_62fbSingleJetMC53/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_JEC_Mikko/plots/KalibriPlots.root");
	  //else if(type==1)f=TFile::Open("../20122012ABC_199429_CORR2012SQLV7_AK5_MC_Su12Z253_kostas_62fbSingleJetMC53/dijetsFall10_TuneZ2_AK5PFCHS_weighted_residuals_JEC_Mikko/plots/KalibriPlots.root");
	  /*if(type==0)f=TFile::Open("../20122012ABC_200601_CORR2012SQLV7_AK5_MC_Su12Z253_kostas_85fbDiJetMC53/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_JEC_Mikko/plots/KalibriPlots.root");
	    else
	    
	    if(type==1)f=TFile::Open("../20122012ABC_200601_CORR2012SQLV7_AK5_MC_Su12Z253_kostas_85fbDiJetMC53/dijetsFall10_TuneZ2_AK5PFCHS_weighted_residuals_JEC_Mikko/plots/KalibriPlots.root");
	  */
	  //if(type==0)f=TFile::Open("../20122012ABC_202305_CORR2012SQLV7_AK5_MC_Su12Z253_Smear_kostas_MPF_2012ABC-SJ-53S-FullwLR/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_JEC_Mikko/plots/KalibriPlots.root");
	  //else if(type==1)f=TFile::Open("../20122012ABC_202305_CORR2012SQLV7_AK5_MC_Su12Z253_Smear_kostas_MPF_2012ABC-SJ-53S-FullwLR/dijetsFall10_TuneZ2_AK5PFCHS_weighted_residuals_JEC_Mikko/plots/KalibriPlots.root");
	  //if(type==0)f=TFile::Open("../20122012ABC_203002_CORR2012FallV1_AK5_MC_Su12Z253_Smear_kostas_12-1fb-DJ-Fall2012V1/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_JEC_Mikko/plots/KalibriPlots.root");
	  //else if(type==1)f=TFile::Open("../20122012ABC_203002_CORR2012FallV1_AK5_MC_Su12Z253_Smear_kostas_12-1fb-DJ-Fall2012V1/dijetsFall10_TuneZ2_AK5PFCHS_weighted_residuals_JEC_Mikko/plots/KalibriPlots.root");	  
	  	 
          if(type==0)f=TFile::Open("../20122012ABCD_ReReco_CORR2012FallV5_AK5_MC_Su12Z253_V11_T1T2_kostas_SJ-ABCD-ReReco/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_JEC_Mikko/plots/KalibriPlots.root");
	  else if(type==1)f=TFile::Open("../20122012ABCD_ReReco_CORR2012FallV5_AK5_MC_Su12Z253_V11_T1T2_kostas_SJ-ABCD-ReReco/dijetsFall10_TuneZ2_AK5PFCHS_weighted_residuals_JEC_Mikko/plots/KalibriPlots.root");	
	  if(verbose)std::cout<<"generate/update Summary.root"<<endl;   
	  //define write-target
	  if(type==0 && choice==0)save=new TFile("Winter12V1_SJ.root", "recreate", "plots", 1);
	  if(type>0 || choice>0)save=new TFile("Winter12V1_SJ.root", "update", "plots", 1);
	  //loop over RR/MPF method
	  for(int method=0; method<2; method++)
	    {
	      if(verbose)std::cout<<"method loop="<<method<<endl; 
	      //loop over alpha bins
	      for(Int_t alpha=0; alpha<4; alpha++)
		{
		  if(verbose)std::cout<<"alpha bin"<<alphaBin[alpha]<<endl; 
		  //loop over eta bins
		  for(Int_t eta=0; eta<7; eta++)
		    {
		      if(verbose)std::cout<<"eta bin"<<etaBin[eta]<<endl; 		  
		      //readout and write histograms
		      TH1F *ReadData, *ReadMC, *WriteRatio;
		      //read histograms
		      if(method==0)
			{
			  if(choice==0) //ptave
			    {
			      ReadData=(TH1F*)f->Get(TString::Format("AbsAsymmetryVsPt%d/AbsAsymmetryVsPt%d_AsymmetryVsMeanPt_data_L2L3_AbsEta%d_RatioOfMeans", alphaBin[alpha], alphaBin[alpha], eta));
			      ReadData->GetYaxis()->SetRangeUser(0.95,1.05);
			      ReadMC=(TH1F*)f->Get(TString::Format("AbsAsymmetryVsPt%d/AbsAsymmetryVsPt%d_AsymmetryVsMeanPt_MC_L2L3_AbsEta%d_RatioOfMeans", alphaBin[alpha], alphaBin[alpha], eta));
			      ReadMC->GetYaxis()->SetRangeUser(0.95,1.05);
			    }
			  else if(choice==1) //leadpt
			    {
			      ReadData=(TH1F*)f->Get(TString::Format("AbsAsymmetryVsLeadPt%d/AbsAsymmetryVsLeadPt%d_AsymmetryVsJetLeadPt_data_L2L3_AbsEta%d_RatioOfMeans", alphaBin[alpha], alphaBin[alpha], eta));
			      ReadData->GetYaxis()->SetRangeUser(0.95,1.05);
			      ReadMC=(TH1F*)f->Get(TString::Format("AbsAsymmetryVsLeadPt%d/AbsAsymmetryVsLeadPt%d_AsymmetryVsJetLeadPt_MC_L2L3_AbsEta%d_RatioOfMeans", alphaBin[alpha], alphaBin[alpha], eta));
			      ReadMC->GetYaxis()->SetRangeUser(0.95,1.05);
			    }	
			  else if(choice==2) //barrelpt
			    {
			      ReadData=(TH1F*)f->Get(TString::Format("AbsAsymmetryVsVBarrPt%d/AbsAsymmetryVsVBarrPt%d_AsymmetryVsJet2Pt_data_L2L3_AbsEta%d_RatioOfMeans", alphaBin[alpha], alphaBin[alpha], eta));
			      ReadData->GetYaxis()->SetRangeUser(0.95,1.05);
			      ReadMC=(TH1F*)f->Get(TString::Format("AbsAsymmetryVsVBarrPt%d/AbsAsymmetryVsVBarrPt%d_AsymmetryVsJet2Pt_MC_L2L3_AbsEta%d_RatioOfMeans", alphaBin[alpha], alphaBin[alpha], eta));
			      ReadMC->GetYaxis()->SetRangeUser(0.95,1.05);
			    }	      
			}
		      else if(method==1)
			{
			  if(choice==0) //ptave
			    {		    
			      ReadData=(TH1F*)f->Get(TString::Format("AbsMPFVsPtAve%d/AbsMPFVsPtAve%d_MPFResponseVsMeanPt_data_L2L3_AbsEta%d_Mean", alphaBin[alpha], alphaBin[alpha], eta));
			      ReadMC=(TH1F*)f->Get(TString::Format("AbsMPFVsPtAve%d/AbsMPFVsPtAve%d_MPFResponseVsMeanPt_MC_L2L3_AbsEta%d_Mean", alphaBin[alpha], alphaBin[alpha], eta));
			    }
			  else if(choice==1) //leadpt
			    {
			      ReadData=(TH1F*)f->Get(TString::Format("AbsMPFVsLeadPt%d/AbsMPFVsLeadPt%d_MPFResponseVsJetLeadPt_data_L2L3_AbsEta%d_Mean", alphaBin[alpha], alphaBin[alpha], eta));
			      ReadMC=(TH1F*)f->Get(TString::Format("AbsMPFVsLeadPt%d/AbsMPFVsLeadPt%d_MPFResponseVsJetLeadPt_MC_L2L3_AbsEta%d_Mean", alphaBin[alpha], alphaBin[alpha], eta));
			    }			
			  else if(choice==2) //barrelpt
			    {
			      ReadData=(TH1F*)f->Get(TString::Format("AbsMPFVsPt%d/AbsMPFVsPt%d_MPFResponseVsJet2Pt_data_L2L3_AbsEta%d_Mean", alphaBin[alpha], alphaBin[alpha], eta));
			      ReadMC=(TH1F*)f->Get(TString::Format("AbsMPFVsPt%d/AbsMPFVsPt%d_MPFResponseVsJet2Pt_MC_L2L3_AbsEta%d_Mean", alphaBin[alpha], alphaBin[alpha], eta));
			    }			
			}
		      //calculate ratio
		      WriteRatio=(TH1F*)ReadData->Clone();
		      WriteRatio->SetDrawOption("P");
		      size_t bins=WriteRatio->GetNbinsX();
		      for(Int_t bin=0; bin<bins; bin++)
			{
			  float ratio;
			  float error;
			  if(ReadMC->GetBinContent(bin+1)>0) //check for filled bins
			    {
			      ratio=ReadData->GetBinContent(bin+1)/ReadMC->GetBinContent(bin+1);
			      error=TMath::Sqrt(TMath::Power(ReadData->GetBinError(bin+1)/ReadMC->GetBinContent(bin+1),2)+TMath::Power(ReadMC->GetBinError(bin+1)*ReadData->GetBinContent(bin+1),2)*TMath::Power(ReadMC->GetBinContent(bin+1),-4));
			      WriteRatio->SetBinContent(bin+1, ratio);
			      WriteRatio->SetBinError(bin+1, error);
			    }
			  else
			    {
			      WriteRatio->SetBinContent(bin+1, 0);
			      WriteRatio->SetBinError(bin+1, 0);
			    }
			}
		      //write histograms
		      if(method==0)
			{
			  if(type==0)
			    {
			      if(choice==0) //ptave
				{
				  ReadData->Write(TString::Format("PtBal_a%d_eta%d_%d_data", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				  ReadMC->Write(TString::Format("PtBal_a%d_eta%d_%d_MC", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));			  
				  WriteRatio->Write(TString::Format("PtBal_a%d_eta%d_%d_ratio", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				}
			      else if(choice==1) //leadpt
				{
				  ReadData->Write(TString::Format("PtBal_VsLeadPt_a%d_eta%d_%d_data", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				  ReadMC->Write(TString::Format("PtBal_VsLeadPt_a%d_eta%d_%d_MC", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));			  
				  WriteRatio->Write(TString::Format("PtBal_VsLeadPt_a%d_eta%d_%d_ratio", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				}
			      else if(choice==2) //barrpt
				{
				  ReadData->Write(TString::Format("PtBal_VsBarrPt_a%d_eta%d_%d_data", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				  ReadMC->Write(TString::Format("PtBal_VsBarrPt_a%d_eta%d_%d_MC", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));			  
				  WriteRatio->Write(TString::Format("PtBal_VsBarrPt_a%d_eta%d_%d_ratio", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				}			    
			    }
			  else if(type==1)
			    {
			      if(choice==0) //ptave
				{
				  ReadData->Write(TString::Format("PtBalchs_a%d_eta%d_%d_data", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				  ReadMC->Write(TString::Format("PtBalchs_a%d_eta%d_%d_MC", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				  WriteRatio->Write(TString::Format("PtBalchs_a%d_eta%d_%d_ratio", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				}
			      else if(choice==1) //leadpt
				{
				  ReadData->Write(TString::Format("PtBalchs_VsLeadPt_a%d_eta%d_%d_data", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				  ReadMC->Write(TString::Format("PtBalchs_VsLeadPt_a%d_eta%d_%d_MC", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				  WriteRatio->Write(TString::Format("PtBalchs_VsLeadPt_a%d_eta%d_%d_ratio", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				}
			      else if(choice==2) //barrpt
				{
				  ReadData->Write(TString::Format("PtBalchs_VsBarrPt_a%d_eta%d_%d_data", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				  ReadMC->Write(TString::Format("PtBalchs_VsBarrPt_a%d_eta%d_%d_MC", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				  WriteRatio->Write(TString::Format("PtBalchs_VsBarrPt_a%d_eta%d_%d_ratio", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				}
			    }		        
			}
		      else if(method==1)
			{
			  if(type==0)
			    {
			      if(choice==0) //ptave
				{
				  ReadData->Write(TString::Format("MPF_VsPtAve_a%d_eta%d_%d_data", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				  ReadMC->Write(TString::Format("MPF_VsPtAve_a%d_eta%d_%d_MC", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				  WriteRatio->Write(TString::Format("MPF_VsPtAve_a%d_eta%d_%d_ratio", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				}
			      else if(choice==1) //leadpt
				{
				  ReadData->Write(TString::Format("MPF_VsLeadPt_a%d_eta%d_%d_data", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				  ReadMC->Write(TString::Format("MPF_VsLeadPt_a%d_eta%d_%d_MC", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				  WriteRatio->Write(TString::Format("MPF_VsLeadPt_a%d_eta%d_%d_ratio", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				}
			      else if(choice==2) //barrpt
				{
				  ReadData->Write(TString::Format("MPF_a%d_eta%d_%d_data", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				  ReadMC->Write(TString::Format("MPF_a%d_eta%d_%d_MC", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				  WriteRatio->Write(TString::Format("MPF_a%d_eta%d_%d_ratio", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				}	
			    }
			  else if(type==1)
			    {
			      if(choice==0) //ptave
				{
				  ReadData->Write(TString::Format("MPFchs_VsPtAve_a%d_eta%d_%d_data", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				  ReadMC->Write(TString::Format("MPFchs_VsPtAve_a%d_eta%d_%d_MC", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				  WriteRatio->Write(TString::Format("MPFchs_VsPtAve_a%d_eta%d_%d_ratio", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				}
			      else if(choice==1) //leadpt
				{
				  ReadData->Write(TString::Format("MPFchs_VsLeadPt_a%d_eta%d_%d_data", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				  ReadMC->Write(TString::Format("MPFchs_VsLeadPt_a%d_eta%d_%d_MC", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				  WriteRatio->Write(TString::Format("MPFchs_VsLeadPt_a%d_eta%d_%d_ratio", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				}
			      else if(choice==2) //barrpt
				{
				  ReadData->Write(TString::Format("MPFchs_a%d_eta%d_%d_data", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				  ReadMC->Write(TString::Format("MPFchs_a%d_eta%d_%d_MC", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				  WriteRatio->Write(TString::Format("MPFchs_a%d_eta%d_%d_ratio", alphaBin[alpha], etaBin[eta], etaBin[eta+1]));
				}	
			    }
			}
		      if(verbose)std::cout<<"finished one iteration"<<endl; 
		    }
		}
	    }
	  save->Close();
	  f->Close();
	}
    }
}
