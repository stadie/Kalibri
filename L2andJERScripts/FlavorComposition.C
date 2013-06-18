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
#include "TGraphAsymmErrors.h"
void FlavorComposition() //type 0: ptave; 1: leadpt; 2: barrelpt
{
  bool verbose=false; //activate for debug mode
  if(verbose)std::cout<<"start the Flavor composition analyzer"<<endl;
  //some general settings
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  //define alpha bins  
  if(verbose)std::cout<<"setting binning"<<endl;
  //loading the files  
  if(verbose)std::cout<<"defining files"<<endl;	//change input files, here			      
  TFile *f=TFile::Open("../20132013ABCD_ReReco_CORR2013SummerV1_AK5_MC_Su12Z253_V11_T1T2_kostas_SummerV3-flavor-physical/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_JEC_Mikko/plots/KalibriPlots.root");
  //TFile *f=TFile::Open("../20132013ABCD_ReReco_CORR2013SummerV1_AK5_MC_Su12Hpp53_kostas_SummerV3-flavor-physical/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_JEC_Mikko/plots/KalibriPlots.root");
  //TFile *f=TFile::Open("../20132013ABCD_ReReco_CORR2013SummerV1_AK5_MC_Su12Z253_V11_T1T2_kostas_DJ-SummerV3-flavor/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_JEC_Mikko/plots/KalibriPlots.root");
  //TFile *f=TFile::Open("../20132013ABCD_ReReco_CORR2013SummerV1_AK5_MC_Su12Hpp53_kostas_DJ-SummerV3-flavor/dijetsFall10_TuneZ2_AK5PF_weighted_residuals_JEC_Mikko/plots/KalibriPlots.root");
  if(verbose)std::cout<<"generate/update ResComp.root"<<endl;   
  //define write-target
  TFile *save=new TFile("FlavorComparison_pythia_physical-extrapol.root", "recreate", "plots", 1);     
  TH1F *Readscheme=(TH1F*)f->Get("AbsFlavorVsPtAve20/AbsFlavorVsPtAve20_2D_AbsEta0_MC_L2L3_X"); //this is just a histogram to get the correct binning, no need to reset this, ever
  //here the alpha and eta bins are defined. Make sure to get these right from the very start
  int alphabins[4]={10,15,20,30};
  double alphaedges[6]={0.,0.075,0.125,0.175,0.225,0.375};
  double etabins[7]={0, 1.305, 1.93, 2.5, 2.964, 3.2, 5.191};
  double edges[Readscheme->GetNbinsX()+1]; //the pt-binning is read out, here
  for(int bin=1; bin<Readscheme->GetNbinsX()+1;bin++)
    {
      edges[bin-1]=Readscheme->GetXaxis()->GetBinLowEdge(bin);
      if(verbose)std::cout<<edges[bin-1]<<std::endl;
    }
  edges[Readscheme->GetNbinsX()]=Readscheme->GetXaxis()->GetBinUpEdge(Readscheme->GetNbinsX());
  if(verbose)std::cout<<edges[Readscheme->GetNbinsX()]<<std::endl;
  //defining the top work histograms
  TH1F *undefined, *gluon, *uds, *c, *b;
  TH3F *Undefinedalpha=new TH3F("Undefinedalpha","undefined flavour", 5, alphaedges, Readscheme->GetNbinsX(),edges, 6, etabins);
  TH3F *Gluonalpha=new TH3F("Gluonalpha","gluon flavour", 5, alphaedges, Readscheme->GetNbinsX(),edges, 6, etabins);
  TH3F *Lightalpha=new TH3F("Lightalpha","uds flavour", 5, alphaedges, Readscheme->GetNbinsX(),edges, 6, etabins);    
  TH3F *Charmalpha=new TH3F("Charmalpha","c flavour", 5, alphaedges, Readscheme->GetNbinsX(),edges, 6, etabins);    
  TH3F *Bottomalpha=new TH3F("Bottomalpha"," flavour", 5, alphaedges, Readscheme->GetNbinsX(),edges, 6, etabins);      
  for(int alpha=0; alpha<4; alpha++) //loop over all alpha bins
    {    
    for(int eta=0; eta<6; eta++) //loop over all eta bins
      {		  
      undefined=new TH1F(TString::Format("undefined%d-a%d",eta,alphabins[alpha]),"undefined flavour",Readscheme->GetNbinsX(),edges);
      undefined->GetXaxis()->SetTitle("#bar{p}_{T}");
      undefined->GetYaxis()->SetTitle("event fraction");
      gluon=new TH1F(TString::Format("gluon%d-a%d",eta,alphabins[alpha]),"gluon flavour",Readscheme->GetNbinsX(),edges);
      gluon->GetXaxis()->SetTitle("#bar{p}_{T}");
      gluon->GetYaxis()->SetTitle("event fraction");
      uds=new TH1F(TString::Format("uds%d-a%d",eta,alphabins[alpha]),"light flavour",Readscheme->GetNbinsX(),edges);
      uds->GetXaxis()->SetTitle("#bar{p}_{T}");
      uds->GetYaxis()->SetTitle("event fraction"); 
      c=new TH1F(TString::Format("c%d-a%d",eta,alphabins[alpha]),"charmed flavour",Readscheme->GetNbinsX(),edges);
      c->GetXaxis()->SetTitle("#bar{p}_{T}");
      c->GetYaxis()->SetTitle("event fraction");
      b=new TH1F(TString::Format("b%d-a%d",eta,alphabins[alpha]),"bottom flavour",Readscheme->GetNbinsX(),edges);
      b->GetXaxis()->SetTitle("#bar{p}_{T}");
      b->GetYaxis()->SetTitle("event fraction");
      undefined->Sumw2();
      gluon->Sumw2();
      uds->Sumw2();
      c->Sumw2();
      b->Sumw2();
      //read histograms
      TH2F *Readundefined, *Readgluon, *Readuds, *Readc, *Readb;
      Readundefined=(TH2F*)f->Get(TString::Format("AbsFlavorVsPtAve%d/AbsFlavorVsPtAve%d_FlavorVsMeanPt_MC_L2L3_AbsEta%d",alphabins[alpha],alphabins[alpha],eta));
      Readgluon=(TH2F*)f->Get(TString::Format("AbsFlavorVsPtAve%d/AbsFlavorVsPtAve%d_FlavorVsMeanPt_MC_L2L3_AbsEta%d",alphabins[alpha],alphabins[alpha],eta));
      Readuds=(TH2F*)f->Get(TString::Format("AbsFlavorVsPtAve%d/AbsFlavorVsPtAve%d_FlavorVsMeanPt_MC_L2L3_AbsEta%d",alphabins[alpha],alphabins[alpha],eta));
      Readc=(TH2F*)f->Get(TString::Format("AbsFlavorVsPtAve%d/AbsFlavorVsPtAve%d_FlavorVsMeanPt_MC_L2L3_AbsEta%d",alphabins[alpha],alphabins[alpha],eta)); 
      Readb=(TH2F*)f->Get(TString::Format("AbsFlavorVsPtAve%d/AbsFlavorVsPtAve%d_FlavorVsMeanPt_MC_L2L3_AbsEta%d",alphabins[alpha],alphabins[alpha],eta));
      TH1D *Pundefined=(TH1D*)Readundefined->ProjectionX(TString::Format("Pundefined%d",1),1,1);
      TH1D *Pgluon=(TH1D*)Readgluon->ProjectionX(TString::Format("Pgluon%d",2),2,2);
      TH1D *Puds=(TH1D*)Readuds->ProjectionX(TString::Format("Puds%d",3),3,3);
      TH1D *Pc=(TH1D*)Readc->ProjectionX(TString::Format("Pc%d",4),4,4);
      TH1D *Pb=(TH1D*)Readb->ProjectionX(TString::Format("Pb%d",5),5,5);
      for(int bins=1; bins<Pundefined->GetNbinsX()+1; bins++) //loop over all pt-bins -> set bins and their errors with normalized contents
        {
	  if(verbose)std::cout<<bins<<endl;
	  double scale=Pundefined->GetBinContent(bins)+Pgluon->GetBinContent(bins)+Puds->GetBinContent(bins)+Pc->GetBinContent(bins)+Pb->GetBinContent(bins); 
	  if(scale>0) //make sure to actually have events (division by zero is evil)
            {
              undefined->SetBinContent(bins, Pundefined->GetBinContent(bins)/scale);
              gluon->SetBinContent(bins, Pgluon->GetBinContent(bins)/scale);
              uds->SetBinContent(bins, Puds->GetBinContent(bins)/scale);
              c->SetBinContent(bins, Pc->GetBinContent(bins)/scale);
              b->SetBinContent(bins, Pb->GetBinContent(bins)/scale);
	      undefined->SetBinError(bins, Pundefined->GetBinError(bins)/scale);
	      gluon->SetBinError(bins, Pgluon->GetBinError(bins)/scale);
	      uds->SetBinError(bins, Puds->GetBinError(bins)/scale);
	      c->SetBinError(bins, Pc->GetBinError(bins)/scale);
	      b->SetBinError(bins, Pb->GetBinError(bins)/scale);
	      Undefinedalpha->SetBinContent(alpha+2, bins, eta+1, Pundefined->GetBinContent(bins)/scale);
	      Gluonalpha->SetBinContent(alpha+2, bins, eta+1, Pgluon->GetBinContent(bins)/scale);
	      Lightalpha->SetBinContent(alpha+2, bins, eta+1, Puds->GetBinContent(bins)/scale);
	      Charmalpha->SetBinContent(alpha+2, bins, eta+1, Pc->GetBinContent(bins)/scale);
	      Bottomalpha->SetBinContent(alpha+2, bins, eta+1, Pb->GetBinContent(bins)/scale);
	      Undefinedalpha->SetBinError(alpha+2, bins, eta+1, Pundefined->GetBinError(bins)/scale);
	      Gluonalpha->SetBinError(alpha+2, bins, eta+1, Pgluon->GetBinError(bins)/scale);
	      Lightalpha->SetBinError(alpha+2, bins, eta+1, Puds->GetBinError(bins)/scale);
	      Charmalpha->SetBinError(alpha+2, bins, eta+1, Pc->GetBinError(bins)/scale);
	      Bottomalpha->SetBinError(alpha+2, bins, eta+1, Pb->GetBinError(bins)/scale);	      
	    }
	}
      //getting rid of unnecessary stuff
      /*delete Readundefined;
      delete Readgluon;
      delete Readuds;
      delete Readc;
      delete Readb;
      delete Pundefined;
      delete Pgluon;
      delete Puds;
      delete Pc;
      delete Pb;*/
      //all the canvas and legends settings
      TCanvas *canvas=new TCanvas(TString::Format("canvas%d-a%d",eta,alphabins[alpha]),TString::Format("%.3f#leq|#eta|<%.3f",etabins[eta],etabins[eta+1]));      
      TLegend *leg=new TLegend(0.7,0.6,0.9,0.9,"flavor");
      undefined->SetLineColor(1);
      undefined->SetMarkerColor(1);
      undefined->SetMarkerSize(2);
      leg->AddEntry(undefined,"undefined", "lp");
      gluon->SetLineColor(2);
      gluon->SetMarkerColor(2);
      gluon->SetMarkerSize(2);
      leg->AddEntry(gluon,"gluon", "lp");     
      uds->SetLineColor(3);
      uds->SetMarkerColor(3);
      uds->SetMarkerSize(3);
      leg->AddEntry(uds,"uds", "lp");
      c->SetLineColor(4);
      c->SetMarkerColor(4);
      c->SetMarkerSize(4);
      leg->AddEntry(c,"c", "lp");
      b->SetLineColor(5);
      b->SetMarkerColor(5);
      b->SetMarkerSize(5);
      leg->AddEntry(b,"b", "lp");      
      leg->SetFillColor(0);
      undefined->SetMinimum(0);
      gluon->SetMinimum(0);  
      uds->SetMinimum(0);  
      c->SetMinimum(0);  
      b->SetMinimum(0);  
      undefined->SetMaximum(1);
      gluon->SetMaximum(1);  
      uds->SetMaximum(1);  
      c->SetMaximum(1);  
      b->SetMaximum(1);  
      //writing everything to file
      undefined->Write();
      gluon->Write();
      uds->Write();
      c->Write();
      b->Write();
      undefined->SetTitle(TString::Format("%.3f#leq|#eta|<%.3f flavor fraction summary",etabins[eta],etabins[eta+1])); 
      undefined->Draw("p");
      gluon->Draw("same,p");
      uds->Draw("same,p");
      c->Draw("same,p");
      b->Draw("same,p");
      leg->Draw(); 
      canvas->Write();
      delete leg;
      delete canvas;
    }//end of eta loop
  }//end of alpha loop

  //do the extrapolations
  for(int eta=0; eta<6; eta++)
    {
      undefined=new TH1F(TString::Format("undefined%d-extrapol",eta),"undefined flavour",Readscheme->GetNbinsX(),edges);
      undefined->GetXaxis()->SetTitle("#bar{p}_{T}");
      undefined->GetYaxis()->SetTitle("event fraction");
      gluon=new TH1F(TString::Format("gluon%d-extrapol",eta),"gluon flavour",Readscheme->GetNbinsX(),edges);
      gluon->GetXaxis()->SetTitle("#bar{p}_{T}");
      gluon->GetYaxis()->SetTitle("event fraction");
      uds=new TH1F(TString::Format("uds%d-extrapol",eta),"light flavour",Readscheme->GetNbinsX(),edges);
      uds->GetXaxis()->SetTitle("#bar{p}_{T}");
      uds->GetYaxis()->SetTitle("event fraction"); 
      c=new TH1F(TString::Format("c%d-extrapol",eta),"charmed flavour",Readscheme->GetNbinsX(),edges);
      c->GetXaxis()->SetTitle("#bar{p}_{T}");
      c->GetYaxis()->SetTitle("event fraction");
      b=new TH1F(TString::Format("b%d-extrapol",eta),"bottom flavour",Readscheme->GetNbinsX(),edges);
      b->GetXaxis()->SetTitle("#bar{p}_{T}");
      b->GetYaxis()->SetTitle("event fraction");
      undefined->Sumw2();
      gluon->Sumw2();
      uds->Sumw2();
      c->Sumw2();
      b->Sumw2();	    
      for(int bins=1; bins<Gluonalpha->GetNbinsY()+1; bins++)
        {
  
	  //make projections
	  TH1F *workUndefined=(TH1F*)Undefinedalpha->ProjectionX(TString::Format("workUndefined%d%d",eta,bins),bins,bins,eta+1,eta+1);
	  TH1F *workGluon=(TH1F*)Gluonalpha->ProjectionX(TString::Format("workGluon%d%d",eta,bins),bins,bins,eta+1,eta+1);
	  TH1F *workLight=(TH1F*)Lightalpha->ProjectionX(TString::Format("workLight%d%d",eta,bins),bins,bins,eta+1,eta+1);
	  TH1F *workCharm=(TH1F*)Charmalpha->ProjectionX(TString::Format("workCharm%d%d",eta,bins),bins,bins,eta+1,eta+1);
	  TH1F *workBottom=(TH1F*)Bottomalpha->ProjectionX(TString::Format("workBottom%d%d",eta,bins),bins,bins,eta+1,eta+1);
	  TF1 *lineU=new TF1("lineU", "[0]*x+[1]", 0, 0.3);
	  TF1 *lineG=new TF1("lineG", "[0]*x+[1]", 0, 0.3);
	  TF1 *lineUDS=new TF1("lineUDS", "[0]*x+[1]", 0, 0.3);
	  TF1 *lineC=new TF1("lineC", "[0]*x+[1]", 0, 0.3);
	  TF1 *lineB=new TF1("lineB", "[0]*x+[1]", 0, 0.3);
	  workUndefined->Fit(lineU);
	  TCanvas *fit=new TCanvas("fit","fit");
	  workGluon->Fit(lineG);
          if(verbose)fit->SaveAs(TString::Format("Output/Fit_eta%d_pt%.0f-%.0f.eps",eta,edges[bins-1],edges[bins]));
	  workLight->Fit(lineUDS);
	  workCharm->Fit(lineC);
	  workBottom->Fit(lineB);
	  double scale=lineU->GetParameter(1)+lineG->GetParameter(1)+lineUDS->GetParameter(1)+lineC->GetParameter(1)+lineB->GetParameter(1);
	  if(verbose)std::cout<<scale<<std::endl;
	  if(scale>0)
	    {
	      undefined->SetBinContent(bins, lineU->GetParameter(1)/scale);
	      gluon->SetBinContent(bins, lineG->GetParameter(1)/scale);
	      uds->SetBinContent(bins, lineUDS->GetParameter(1)/scale);
	      c->SetBinContent(bins, lineC->GetParameter(1)/scale);
	      b->SetBinContent(bins, lineB->GetParameter(1)/scale);
	      undefined->SetBinError(bins, lineU->GetParError(1)/scale);
	      gluon->SetBinError(bins, lineG->GetParError(1)/scale);
	      uds->SetBinError(bins, lineUDS->GetParError(1)/scale);
	      c->SetBinError(bins, lineC->GetParError(1)/scale);
	      b->SetBinError(bins, lineB->GetParError(1)/scale);	
	    }  
	}
      //all the canvas and legends settings
      TCanvas *canvas=new TCanvas(TString::Format("canvas%d-extrapol",eta),TString::Format("%.3f#leq|#eta|<%.3f",etabins[eta],etabins[eta+1]));      
      TLegend *leg=new TLegend(0.7,0.6,0.9,0.9,"flavor");
      undefined->SetLineColor(1);
      undefined->SetMarkerColor(1);
      undefined->SetMarkerSize(2);
      leg->AddEntry(undefined,"undefined", "lp");
      gluon->SetLineColor(2);
      gluon->SetMarkerColor(2);
      gluon->SetMarkerSize(2);
      leg->AddEntry(gluon,"gluon", "lp");     
      uds->SetLineColor(3);
      uds->SetMarkerColor(3);
      uds->SetMarkerSize(3);
      leg->AddEntry(uds,"uds", "lp");
      c->SetLineColor(4);
      c->SetMarkerColor(4);
      c->SetMarkerSize(4);
      leg->AddEntry(c,"c", "lp");
      b->SetLineColor(5);
      b->SetMarkerColor(5);
      b->SetMarkerSize(5);
      leg->AddEntry(b,"b", "lp");      
      leg->SetFillColor(0);
      undefined->SetMinimum(0);
      gluon->SetMinimum(0);  
      uds->SetMinimum(0);  
      c->SetMinimum(0);  
      b->SetMinimum(0);  
      undefined->SetMaximum(1);
      gluon->SetMaximum(1);  
      uds->SetMaximum(1);  
      c->SetMaximum(1);  
      b->SetMaximum(1);  
      //writing everything to file
      undefined->Write();
      gluon->Write();
      uds->Write();
      c->Write();
      b->Write();
      undefined->SetTitle(TString::Format("%.3f#leq|#eta|<%.3f flavor fraction summary",etabins[eta],etabins[eta+1])); 
      undefined->Draw("p");
      gluon->Draw("same,p");
      uds->Draw("same,p");
      c->Draw("same,p");
      b->Draw("same,p");
      leg->Draw(); 
      canvas->Write();
      delete leg;
      delete canvas;	
    }
  //get rid of unneeded histograms
  delete undefined;
  delete gluon;
  delete uds;
  delete c;
  delete b;
  save->Close();
  f->Close();
}
