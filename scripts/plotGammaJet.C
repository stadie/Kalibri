void plotGammaJet() {
  gStyle->SetPalette(1);
  TFile* f = TFile::Open("/scratch/current/cms/user/stadie/toy_photonjet2.root");
  
  TTree* tree = (TTree*)f->Get("GammaJetTree");

  double bins[] = {0,5,10,25,50,400};
  new TCanvas();
  TH1F* hntower = new TH1F("hntower",";E^{Had} [GeV];N_{tower}/Jet",5,bins);
  tree->Draw("TowHad >> hntower");
  hntower->SetNormFactor(hntower->Integral()/tree->GetEntries());
  hntower->DrawCopy(); 
  new TCanvas();
  TProfile* htowprof = new TProfile("htowprof",";E^{Had} [GeV];#frac{E_{T}^{Jet}}{E_{T}^{#gamma}} -1",
				    5,bins,-1,1);
  tree->Draw("JetCalPt/PhotonPt - 1:TowHad >> htowprof");
  htowprof->DrawCopy();
  new TCanvas();
  TH2F* htowresp = new TH2F("htowresp",";E^{Had} [GeV];#frac{E_{T}^{Jet}}{E_{T}^{#gamma}} -1",
			    5,bins,100,-0.1,0.1);
  tree->Draw("JetCalPt/PhotonPt -1:TowHad >> htowresp");
  gPad->SetLogz();
  htowresp->DrawCopy("COLZ");
  htowprof->DrawCopy("SAME");
  new TCanvas();
  TProfile2D* htowjetprof = new TProfile2D("htowjetprof",";E_{T}^{Jet};E^{Had} [GeV];#frac{E_{T}^{Jet}}{E_{T}^{#gamma}}",40,0,400,5,bins);
  tree->Draw("JetCalPt/PhotonPt:TowHad:JetCalE >> htowjetprof","JetCalE < 400");
  htowjetprof->SetMinimum(0.8);
  htowjetprof->Draw("COLZ");

  new TCanvas();
  TProfile* hjetprof = new TProfile("hjetprof",";E_T^{Jet} [GeV];#frac{E_{T}^{Jet}}{E_{T}^{#gamma}} -1",
				    40,0,400,-1,1);
  tree->Draw("JetCalPt/PhotonPt-1:JetCalPt >> hjetprof","JetCalE < 400");
  hjetprof->Draw();

  new TCanvas();
  TProfile2D* htowjetprof2 = new TProfile2D("htowjetprof2",";E_{T}^{Jet};E^{Had} [GeV];#frac{E_{T}^{Jet}}{E_{T}^{#gamma}}",40,0,400,5,bins);
  tree->Draw("JetCalPt/PhotonPt:TowHad:JetCalE >> htowjetprof2","JetCalE < 400 && JetCalPt > 40 && JetCalPt < 350");
  htowjetprof2->SetMinimum(0.8);
  htowjetprof2->Draw("COLZ");
  new TCanvas();
  TProfile* htowprof2 = new TProfile("htowprof2",";E^{Had} [GeV];#frac{E_{T}^{Jet}}{E_{T}^{#gamma}} -1",
				    5,bins,-1,1);
  tree->Draw("JetCalPt/PhotonPt - 1:TowHad >> htowprof2","JetCalE < 400 && JetCalPt > 40 && JetCalPt < 350");
  htowprof2->DrawCopy();
  new TCanvas();
  TH1F* hntower2 = new TH1F("hntower2","#frac{E_{T}^{Jet}}{E_{T}^{#gamma} < 0.9;E^{Had} [GeV];N_{tower}/Jet",5,bins);
  tree->Draw("TowHad >> hntower2","JetCalE > 200 && JetCalPt/PhotonPt < 0.9");
  hntower2->SetNormFactor(hntower2->Integral()/tree->GetEntries("JetCalE > 200 && JetCalPt/PhotonPt < 0.9"));
  hntower2->DrawCopy(); 
 
  new TCanvas();
  TH1F* hntower3 = new TH1F("hntower3","#frac{E_{T}^{Jet}}{E_{T}^{#gamma} > 1.1;E^{Had} [GeV];N_{tower}/Jet",5,bins);
  tree->Draw("TowHad >> hntower3","JetCalE > 200 && JetCalPt/PhotonPt > 1.1");
  hntower3->SetNormFactor(hntower3->Integral()/tree->GetEntries("JetCalE > 200 && JetCalPt/PhotonPt > 1.1"));
  hntower3->DrawCopy(); 

  new TCanvas();
  TProfile* hjetprof2 = new TProfile("hjetprof2",";E_{T}^{Jet};#frac{E_{T}^{Jet}}{E_{T}^{#gamma}} -1",
				    40,0,400,-1,1);
  tree->Draw("JetCalPt/PhotonPt - 1:JetCalPt >> hjetprof2","JetCalE < 400 && TowHad < 5");
  hjetprof2->DrawCopy();
  

  new TCanvas();
  TProfile* hntowjet = new TProfile("hntowjet",";E_{T}^{Jet};n_{tower}",40,0,400);
  
  tree->Draw("Length$(TowHad):JetCalPt>> hntowjet","JetCalPt/PhotonPt < 0.9");
  hntowjet->Draw();


  //only for MC 
  TCanvas* c = new TCanvas();
  c->Divide(3,2);
  TH1F* htowert = new TH1F("htowert",";E^{Had}_{true}/E^{Had} all;",100,0,2);
  TH1F* htower0t = new TH1F("htower0t",";E^{Had}_{true}/E^{Had} true bin 1;",100,0,2);
  TH1F* htower1t = new TH1F("htower1t",";E^{Had}_{true}/E^{Had} true bin 2;",100,0,2);
  TH1F* htower2t = new TH1F("htower2t",";E^{Had}_{true}/E^{Had} true bin 3;",100,0,2);
  TH1F* htower3t = new TH1F("htower3t",";E^{Had}_{true}/E^{Had} true bin 4;",100,0,2);
  TH1F* htower4t = new TH1F("htower4t",";E^{Had}_{true}/E^{Had} true bin 5;",100,0,2);
  c->cd(1);
  tree->Draw("TowHadTrue/TowHad >> htowert");
  htowert->SetMinimum(0);
  htowert->Fit("gaus");
  c->cd(2);
  tree->Draw("TowHadTrue/TowHad >> htower0t","TowHadTrue < 5");
  htower0t->SetMinimum(1);
  htower0t->Fit("gaus");
  c->cd(3);
  tree->Draw("TowHadTrue/TowHad >> htower1t","TowHadTrue > 5 && TowHadTrue < 10");
  htower1t->SetMinimum(1);
  htower1t->Fit("gaus");
  c->cd(4);
  tree->Draw("TowHadTrue/TowHad >> htower2t","TowHadTrue > 10 && TowHadTrue < 25");
  htower2t->SetMinimum(1);
  htower2t->Fit("gaus");
  c->cd(5);
  tree->Draw("TowHadTrue/TowHad >> htower3t","TowHadTrue > 25 && TowHadTrue < 50");
  htower3t->SetMinimum(1);
  htower3t->Fit("gaus");
  c->cd(6);  
  tree->Draw("TowHadTrue/TowHad >> htower4t","TowHadTrue > 50");
  htower4t->SetMinimum(1);
  htower4t->Fit("gaus");
  c = new TCanvas();
  c->Divide(3,2);
  TH1F* htower = new TH1F("htower",";E^{Had}_{true}/E^{Had} all;",100,0,2);
  TH1F* htower0 = new TH1F("htower0",";E^{Had}_{true}/E^{Had} bin 1;",100,0,2);
  TH1F* htower1 = new TH1F("htower1",";E^{Had}_{true}/E^{Had} bin 2;",100,0,2);
  TH1F* htower2 = new TH1F("htower2",";E^{Had}_{true}/E^{Had} bin 3;",100,0,2);
  TH1F* htower3 = new TH1F("htower3",";E^{Had}_{true}/E^{Had} bin 4;",100,0,2);
  TH1F* htower4 = new TH1F("htower4",";E^{Had}_{true}/E^{Had} bin 5;",100,0,2);
  c->cd(1);
  tree->Draw("TowHadTrue/TowHad >> htower");
  htower->SetMinimum(0);
  htower->Fit("gaus");
  c->cd(2);
  tree->Draw("TowHadTrue/TowHad >> htower0","TowHad < 5");
  htower0->SetMinimum(1);
  htower0->Fit("gaus");
  c->cd(3);
  tree->Draw("TowHadTrue/TowHad >> htower1","TowHad > 5 && TowHad < 10");
  htower1->SetMinimum(1);
  htower1->Fit("gaus");
  c->cd(4);
  tree->Draw("TowHadTrue/TowHad >> htower2","TowHad > 10 && TowHad < 25");
  htower2->SetMinimum(1);
  htower2->Fit("gaus");
  c->cd(5);
  tree->Draw("TowHadTrue/TowHad >> htower3","TowHad > 25 && TowHad < 50");
  htower3->SetMinimum(1);
  htower3->Fit("gaus");
  c->cd(6);  
  tree->Draw("TowHadTrue/TowHad >> htower4","TowHad > 50");
  htower4->SetMinimum(1);
  htower4->Fit("gaus");
}
