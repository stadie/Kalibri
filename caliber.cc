//
// Original Author:  Christian Autermann
//         Created:  Wed Jul 18 13:54:50 CEST 2007
// $Id: caliber.cc,v 1.67 2008/12/11 17:20:24 stadie Exp $
//
//
// for profiling:
//  - to prevent gprof from missing the threads: 
//      wget http://sam.zoy.org/writings/programming/gprof-helper.c
//      gcc -shared -fPIC gprof-helper.c -o gprof-helper.so -lpthread -ldl 
//      LD_PRELOAD=./gprof-helper.so ./junk
//
#include "caliber.h"

//C++ libs
#include <cmath>
#include <iomanip>
#include <set>
#include <time.h>

//Root libs
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TStyle.h"
//Boost
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
boost::mutex io_mutex;
// User
#include "ConfigFile.h"
#include "Parameters.h"
#include "ControlPlots.h"
#include "CalibMath.h"
#include "external.h"
#include "ToyMC.h"
#include "PhotonJetReader.h"
#include "DiJetReader.h"
#include "TriJetReader.h"

#include "TH1F.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TText.h"
#include "TCanvas.h"
#include "TPostScript.h"

using namespace std;

typedef std::vector<TData*>::iterator DataIter;
typedef std::vector<TData*>::const_iterator DataConstIter;

//Outlier Rejection
struct OutlierRejection {
  OutlierRejection(double cut):_cut(cut){};
  bool operator()(TData *d){
    if(d->GetType()==typeTowerConstraint) return true;
    return (d->chi2()/d->GetWeight())<_cut;
  }
  double _cut;
};

/* 
//ControlCut Selection
struct ControlCutSelection {
ControlCutSelection(double cut):_cut(cut){};
  bool operator()(TData *d){
    return d->GetTruth()>_cut && fabs(d->GetMess()->eta)<2.5;
  }
  double _cut;
};

*/
//"Not-Balanced" Rejection: Make average-fitting equal to peak-fitting
struct NotBalancedRejection {
  NotBalancedRejection(double *cut, double min, double max):
    _cut(cut),_min(min),_max(max){};
  bool operator()(TData *d){
    bool result = false;
    TAbstractData *ad = dynamic_cast<TAbstractData*>(d);
    if(! ad) return true;
    
    if(ad->GetType()!=GammaJet ||
       ad->GetMess()->pt<_min ||
       ad->GetMess()->pt>_max ||
       ad->GetMess()->pt==0.0 ) result = true;
    else
      result = (1.0-ad->GetTruth()/ad->GetMess()->pt) >
               _cut[(int)(ad->GetMess()->pt-_min)];
    return result;
  }
  double *_cut;
  double _min, _max;
};


class ComputeThread {
private:
  int npar;
  double chi2;
  double * td1;
  double * td2;
  double *parorig, *mypar;
  double *temp_derivative1;
  double *temp_derivative2;
  double epsilon;
  std::vector<TData*> data;
  bool data_changed;
  struct calc_chi2_on
  {
  private:
    ComputeThread *parent;
  public:
    calc_chi2_on(ComputeThread *parent) : parent(parent) {}
    void operator()()
    {
//       {
// 	boost::mutex::scoped_lock lock(io_mutex);
// 	std::cout << "start Thread for " << parent << std::endl; 
//       }   
      if(parent->data_changed) {
	for (DataIter it=parent->data.begin() ; it!= parent->data.end() ; ++it)
	  (*it)->ChangeParAddress(parent->parorig,parent->mypar); 
	parent->data_changed = false;
      }
      for (int param=0; param< parent->npar ; ++param) {
	parent->td1[param]= 0.0;
	parent->td2[param]= 0.0;
	parent->mypar[param] = parent->parorig[param];
      }
      parent->chi2 =0.0;   
      for (DataIter it=parent->data.begin() ; it!= parent->data.end() ; ++it) {
	parent->chi2 += (*it)->chi2_fast(parent->td1, parent->td2, parent->epsilon);
      } 
      boost::mutex::scoped_lock lock(io_mutex);
      for (int param=0; param< parent->npar ; ++param) {
	parent->temp_derivative1[param] += parent->td1[param];
	parent->temp_derivative2[param] += parent->td2[param];
      }
      //std::cout << "stop Thread with for " << parent << std::endl;
    }
  };
  boost::thread *thread;
  friend class calc_chi2_on;
public:
  ComputeThread(int npar,double *par, double *temp_derivative1, double *temp_derivative2, double epsilon) 
    : npar(npar), td1(new double[npar]), td2(new double[npar]), parorig(par),
      mypar(new double[npar]), temp_derivative1(temp_derivative1), 
      temp_derivative2(temp_derivative2), epsilon(epsilon), data_changed(false) {}
  ~ComputeThread() {
    ClearData();
    delete [] td1;
    delete [] td2;
    delete [] mypar;
  }
  void AddData(TData* d) { 
    //d->ChangeParAddress(parorig, mypar);
    data_changed = true;
    data.push_back(d);
  }
  void ClearData() {   
    for (DataIter it= data.begin() ; it!= data.end() ; ++it)  
      (*it)->ChangeParAddress(mypar,parorig);
    data.clear();
  }
  void Start() { thread = new boost::thread(calc_chi2_on(this)); }
  bool IsDone() { thread->join(); delete thread; return true;}
  void SyncParameters() {
    for (int param=0; param< npar ; ++param) mypar[param] = parorig[param];
  }
  double Chi2() const { return chi2;}
};


int TCaliber::GetSpectraBin(double m1, double m2=0., double m3=0.)
{
   //pt
  int bin1, bins1 = 7000; 
   double min1    = 0.;
   double max1    = 7000.;
   if      (m1<min1) bin1=0;
   else if (m1>max1) bin1=bins1+1;
   else              bin1=(int)(((m1-min1)/max1)*(double)bins1);
   //eta
   int bin2=0, bins2 = p->GetEtaGranularityJet()*p->GetPhiGranularityJet();
   double min2    = 0.;
   double max2    = 82.;
   if      (m2<min2) bin2=0;
   else if (m2>max2) bin2=bins2;
   else              bin2=(int)(((m2-min2)/max2)*(double)bins2);
   //EMF
   int bin3=0, bins3 = 100;
   double min3    = 0.;
   double max3    = 1.;
   if      (m3<min3) bin3=0;
   else if (m3>max3) bin3=bins3;
   else              bin3=(int)(((m3-min3)/max3)*(double)bins3);

   return bin1 *bins2*bins3 + bin2 *bins3 + bin3;
}

double gauss_step(double *x, double *par)
{
   return par[2]/(par[1]*2.5)*exp(-(x[0]-par[0])*(x[0]-par[0])/(2.0*par[1]*par[1]))
	   *
	  (1.0-1.0/(1.0+exp((par[3]-x[0])/par[4]) ) );  	
}

void FitWithoutBottom(TH1 * hist, TF1 * func, double bottom=0.33)
{
  TH1F * result=(TH1F*)hist->Clone();
  double maximum = hist->GetMaximum();
  int min=0, max=0;
  for (int i=0; i<hist->GetNbinsX(); ++i)
    if (hist->GetBinContent(i)>bottom*maximum){
      result->SetBinContent(i,hist->GetBinContent(i));
      max=i;
      if (min==0.) min=i;
    }      
  func->SetRange(hist->GetXaxis()->GetXmin()+(double)min/(double)hist->GetNbinsX()*(hist->GetXaxis()->GetXmax()-hist->GetXaxis()->GetXmin()),
                    hist->GetXaxis()->GetXmin()+(double)max/(double)hist->GetNbinsX()*(hist->GetXaxis()->GetXmax()-hist->GetXaxis()->GetXmin()));
  result->Fit("gauss_step","LLQNO","");
}

void TCaliber::FlattenSpectra()
{
 double allweights=0;
 //@@ Replace "7" with some more meaningful variable name!!!
 for(int i=0;i<7;++i)
   allweights += RelWeight[i];

 map<int,double> weights[7];
 double tot[7];
 double alltotal=0;
 for(int i=0;i<7;++i)
   tot[i]=0.;
 
   for (DataConstIter it = data.begin(); it!=data.end(); ++it) {
     TAbstractData* dt = dynamic_cast<TAbstractData*>(*it);
     if(! dt) continue;
     alltotal+=dt->GetWeight();
     for (int type=0; type<7; ++type){
       if (dt->GetType()!=type) continue;
       if (dt->GetType()==InvMass) continue;
       double em = 0.;
       double had = 0.;
       int index=0;
       double min_tower_dr = 10.;
       //double ptmax=0;
       
       TLorentzVector Ljet(0.,0.,0.,0.);
       Ljet.SetPtEtaPhiE(dt->GetMess()->pt,dt->GetMess()->eta,dt->GetMess()->phi,dt->GetMess()->E);
       for(std::vector<TAbstractData*>::const_iterator t=dt->GetRef().begin(); t!=dt->GetRef().end(); ++t) {
	 em  += (*t)->GetMess()->EMF;
	 had += (*t)->GetMess()->HadF;
	 had += (*t)->GetMess()->OutF;
	 TLorentzVector Ltower(0.,0.,0.,0.);
	 Ltower.SetPtEtaPhiE((*t)->GetMess()->pt,(*t)->GetMess()->eta,(*t)->GetMess()->phi,(*t)->GetMess()->E);
	 double dr = Ltower.DeltaR(Ljet);
	 if (dr<min_tower_dr) {
	   index = (*t)->GetIndex();
	   min_tower_dr = dr;
	 }
       }
       //int bin = GetSpectraBin( dt->GetScale(), index, em/(em+had)  );
       //int bin = GetSpectraBin( dt->GetScale(), index );
       int bin = GetSpectraBin( dt->GetScale() );
       double error = 1.;//dt->GetParametrizedErr( &dt->GetMess()->pt );
       weights[type][bin]+=dt->GetWeight()/error;
       tot[type]+=dt->GetWeight()/error;
     }
   }
   for (int type=0; type<7; ++type){
     if (tot[type]!=0.)
       for (DataIter it = data.begin(); it!=data.end(); ++it) {
	 TAbstractData* dt = dynamic_cast<TAbstractData*>(*it);
	 if(! dt) continue;
	 if (dt->GetType()!=type) continue;
         if (dt->GetType()==InvMass) continue;
	 
	 
	 double em = 0;
	 double had = 0;
	 int index=0;
	 double min_tower_dr = 10.;
	 TLorentzVector Ljet(0,0,0,0);
	 
	 Ljet.SetPtEtaPhiE(dt->GetMess()->pt,dt->GetMess()->eta,dt->GetMess()->phi,dt->GetMess()->E);
	 for(std::vector<TAbstractData*>::const_iterator t=dt->GetRef().begin(); t!=dt->GetRef().end(); ++t) {
	   em  += (*t)->GetMess()->EMF;
	   had += (*t)->GetMess()->HadF;
	   had += (*t)->GetMess()->OutF;
	   TLorentzVector Ltower(0.,0.,0.,0.);
	   Ltower.SetPtEtaPhiE((*t)->GetMess()->pt,(*t)->GetMess()->eta,(*t)->GetMess()->phi,(*t)->GetMess()->E);
	   double dr = Ltower.DeltaR(Ljet);
	   if (dr<min_tower_dr) {
	     index = (*t)->GetIndex();
	     min_tower_dr = dr;
	   }
	 }
	 //int bin = GetSpectraBin( dt->GetScale(), index, em/(em+had) );
	 //int bin = GetSpectraBin( dt->GetScale(), index );
	 int bin = GetSpectraBin( dt->GetScale() );
	 //dt->SetWeight(1);
	 
	 
	 //dt->SetWeight(dt->GetWeight()/weights[type][bin] * (double(tot[type]) / double(weights[type].size())));
	 
	 dt->SetWeight(dt->GetWeight()/weights[type][bin] * (alltotal / double(weights[type].size()) )  * RelWeight[type] / allweights );
	 
	 
	 //dt->SetWeight((1./weights[bin]) * (double(tot) / weights.size()));
       }
     
    /*
   // Old version using leading tower bin as jet bin
   for (DataConstIter it = data.begin(); it!=data.end(); ++it) {
   if (dt->GetType()!=type) continue;
   double em = 0;
   double had = 0;
   int index=0;
   double ptmax=0;
   for(std::vector<TData*>::const_iterator t=dt->GetRef().begin(); t!=dt->GetRef().end(); ++t) {
   em  += (*t)->GetMess()[1];
	had += (*t)->GetMess()[2];
	had += (*t)->GetMess()[3];
	if ((*t)->GetMess()[1]>ptmax) {
	ptmax=(*t)->GetMess()[1];
	index=(*t)->GetIndex();
	}
	}
	//int bin = GetSpectraBin( dt->GetScale(), index, em/(em+had)  );
	int bin = GetSpectraBin( dt->GetScale(), index );
	//int bin = GetSpectraBin( dt->GetScale() );
	weights[bin]+=dt->GetWeight();
	tot+=dt->GetWeight();
	}
	
	if (tot!=0.)
    for (DataIter it = data.begin(); it!=data.end(); ++it) {
    if (dt->GetType()!=type) continue;
    
      double em = 0;
      double had = 0;
      int index=0;
      double ptmax=0;
      for(std::vector<TData*>::const_iterator t=dt->GetRef().begin(); t!=dt->GetRef().end(); ++t) {
	em  += (*t)->GetMess()[1];
	had += (*t)->GetMess()[2];
	had += (*t)->GetMess()[3];
	if ((*t)->GetMess()[1]>ptmax) {
	  ptmax=(*t)->GetMess()[1];
	  index=(*t)->GetIndex();
	}
      }

      //int bin = GetSpectraBin( dt->GetScale(), index, em/(em+had) );
      int bin = GetSpectraBin( dt->GetScale(), index );
      //int bin = GetSpectraBin( dt->GetScale() );
      //dt->SetWeight(1);
      //dt->SetWeight(1./weights[bin]);
      dt->SetWeight((1./weights[bin]) * (double(tot) / weights.size()));
    }
    */

   } 
}
  
//further weighting.............................................................
void TCaliber::BalanceSpectra()
{
  cout<<"...further weighting"<<endl;
  double min = Et_cut_on_gamma;
  double max = 100.; //GeV
  int nbins = (int)(max-min);//one bin per GeV
  if (nbins<2) return;
  double EMF[nbins];
  double TOT[nbins];
  
  TCanvas * c1 = new TCanvas("controlplots","",600,600);
  TPostScript ps("balance_spectra.ps",111);
  TH1F * gauss_forpt[nbins];
  TH1F * gauss_forpt_truth[nbins];
  gauss_forpt[0] = new TH1F("hgauss","pT bin[20..21GeV];#frac{pT jet - pT truth}{pT jet}",600,-3,3);
  gauss_forpt_truth[0] = new TH1F("hgauss_truth","pT bin[20..21GeV];pT truth",400,0,200);
  char * name = new char[100];
  for(int i = 1 ; i < nbins; ++i) {
    gauss_forpt[i] = (TH1F*)gauss_forpt[0]->Clone();
    sprintf(name,"pT bin[%d..%dGeV]",(int)min+i,(int)min+i+1);
    gauss_forpt[i]->SetTitle(name);
    gauss_forpt_truth[i] = (TH1F*)gauss_forpt_truth[0]->Clone();
    sprintf(name,"pT bin[%d..%dGeV]",(int)min+i,(int)min+i+1);
    gauss_forpt_truth[i]->SetTitle(name);
  }

  cout<<"...fill truth histograms for each jet-pT bin"<<endl;
  //loop over all fit-events
  for ( std::vector<TData*>::iterator i = data.begin(); 
        i != data.end() ; ++i )  {
    TAbstractData* jg = dynamic_cast<TAbstractData*>(*i);
    if(! jg) continue;
    if (jg->GetType()!=GammaJet) continue;
    
    //double etjetcor = jg->GetParametrizedMess();
    if(jg->GetMess()->pt>min && jg->GetMess()->pt<max) {
      gauss_forpt[(int)(jg->GetMess()->pt-min)]->Fill( (jg->GetMess()->pt-jg->GetTruth())/jg->GetMess()->pt,jg->GetWeight() );
      gauss_forpt_truth[(int)(jg->GetMess()->pt-min)]->Fill( jg->GetTruth(),jg->GetWeight() );
      EMF[(int)(jg->GetMess()->pt-min)] += jg->GetWeight()*jg->GetMess()->EMF;
      TOT[(int)(jg->GetMess()->pt-min)] += jg->GetWeight()*jg->GetMess()->pt;
    }      
  }

cout<<"...fit the truth distributions"<<endl;
  double edge;
  TF1 * f = new TF1("gauss_step",gauss_step,-3,3,5);
  double * cuts = new double[nbins];
  TText * text = new TText();
  text->SetTextSize(0.03);
  text->SetTextColor(2);

  for(int i = 0; i < nbins; ++i) {
    if ( (i+1)%(nbins/10)==0) cout << (100*(i+1)/nbins)<<"% events weighted" << endl;  
    //TF1 *f=0;
    //gauss_forpt[i]->Fit("gaus","LLQNO","");
    //f = (TF1*)gROOT->GetFunction("gaus")->Clone();
    edge = 1.0-Et_cut_on_gamma/(((double)i)+min+0.5);
    f->SetParameters(-1.,2.0,3.0, edge, 0.0001);
    f->FixParameter(3, edge);
    f->FixParameter(4, 0.0001);
    FitWithoutBottom(gauss_forpt[i], f);
    //bla->Fit("gauss_step","LLQNO","");
    //delete bla;

    gauss_forpt[i]->Draw("h");
    f->SetLineColor(2);
    f->Draw("same");
    sprintf(name,"mean %f",f->GetParameter(0));
    text->DrawText(1.4,0.7*gauss_forpt[i]->GetMaximum(),name);

    sprintf(name,"average truth %f", (double)i+0.5+min-((double)i+0.5+min)*f->GetParameter(0));
    text->DrawText(0.50,0.65*gauss_forpt[i]->GetMaximum(),name);
    c1->Draw();

    if (TOT[i]!=0.0) EMF[i] = EMF[i]/TOT[i];
    sprintf(name,"average/truth %+f",(double)i+0.5+min-((double)i+0.5+min)*f->GetParameter(0) / 
                                     (1.3*((double)i+0.5+min)   )); //-EMF[i])  );
    text->DrawText(0.50,0.60*gauss_forpt[i]->GetMaximum(),name);
    c1->Draw();


    gauss_forpt_truth[i]->Draw("h");
    c1->Draw();
    
    cout<<"bin "<<i
      <<": mean="<<f->GetParameter(0)
      <<", sigma="<<f->GetParameter(1)
      <<", height="<<f->GetParameter(2)
      <<", edge("<<edge<<")="<<f->GetParameter(3)
      <<", width-edge="<<f->GetParameter(4)
      <<endl;
    cuts[i] = f->GetParameter(0)-fabs(f->GetParameter(0)-f->GetParameter(3))+0.2;
  }
  delete f;
  ps.Close();
 
  DataIter beg = partition(data.begin(), data.end(), 
                           NotBalancedRejection(cuts, min, max));
  cout<<"...remove " << int(data.end()-beg) << " events which are not 'balanced'";
  for(DataIter i = beg ; i != data.end() ; ++i) {
    cout<<".";
    delete *i;
  }
  cout<<endl;
  data.erase(beg,data.end());

cout<<"...cleaning up"<<endl;
  for(int i = 0; i < nbins; ++i){
    delete gauss_forpt[i];
  }  
  delete [] cuts;
  delete name;
}

void TCaliber::Run_ZJet()
//calculates from Z energy a truth value for one calo tower of the jet.
{
  //Run Z-Jet stuff  
  int nevent = zjet.fChain->GetEntries();
  for (int i=0;i<nevent;i++) {
    if(i%1000==0) cout<<"Z-Jet Event: "<<i<<endl;
    zjet.fChain->GetEvent(i); 
    if (zjet.NobjTowCal>200) {
      cerr<<"ERROR: Increase array sizes in ZJetSelector; NobjTowCal="
	  <<zjet.NobjTowCal<<"!"<<endl;
      exit(8);
    }
 
    //trivial cuts
    if (zjet.ZPt<Et_cut_on_Z || zjet.JetCalPt<Et_cut_on_jet) continue;
     
    //Find the jets eta & phi index using the nearest tower to jet axis:
    int jet_index=-1;
    double min_tower_dr = 10.0;
    double em = 0;
    double had = 0;
    TLorentzVector Ljet(0,0,0,0);
    Ljet.SetPtEtaPhiE(zjet.JetCalEt,zjet.JetCalEta,zjet.JetCalPhi,zjet.JetCalE);
    for (int n=0; n<zjet.NobjTowCal; ++n){
      em += zjet.TowEm[n];
      had +=  zjet.TowHad[n];
      TLorentzVector Ltower(0,0,0,0);
      Ltower.SetPtEtaPhiE(zjet.TowEt[n],zjet.TowEta[n],zjet.TowPhi[n],zjet.TowE[n]);
      double dr = Ltower.DeltaR(Ljet);
      if (dr<min_tower_dr) {
	jet_index = p->GetJetBin(p->GetJetEtaBin(zjet.TowId_eta[n]),
				 p->GetJetPhiBin(zjet.TowId_phi[n]));
	min_tower_dr = dr;
      }
    }

    if (jet_index<0){ cerr<<"WARNING: jet_index = " << jet_index << endl; continue; }
    if(em == 0) { continue;}
    //jet_index: p->eta_granularity*p->phi_granularity*p->GetNumberOfTowerParametersPerBin()
    //           has to be added for a correct reference to k[...].

    TMeasurement* jetp  = new TJet;
    jetp->pt  = zjet.JetCalEt;
    jetp->eta = zjet.JetCalEta;
    jetp->phi = zjet.JetCalPhi;
    jetp->E   = zjet.JetCalE;
    //Create an Z/Jet TData event
    TData_TruthMultMess * gj_data = new 
      TData_TruthMultMess(jet_index  * p->GetNumberOfJetParametersPerBin() + p->GetNumberOfTowerParameters(),
			  // zjet.ZEt,				    //truth//
			  zjet.JetGenPt,
			  sqrt(pow(0.5,2)+pow(0.10*zjet.ZEt,2)),    //error//
			  //zjet.EventWeight,                       //weight//
			  1.0,                                      //weight//
			  p->GetJetParRef( jet_index ),             //params
			  p->GetNumberOfJetParametersPerBin(),      //number of free jet param. p. bin
			  p->jet_parametrization,                   //function
			  jet_error_param,                          //error param. function
			  jetp
			  );
    
    //Add the jet's towers to "gj_data":
    for (int n=0; n<zjet.NobjTowCal; ++n){
      //if (zjet.TowEt[n]<0.01) continue;
      
      int index = p->GetBin(p->GetEtaBin(zjet.TowId_eta[n]),
			    p->GetPhiBin(zjet.TowId_phi[n]));
      if (index<0){ cerr<<"WARNING: towewer_index = " << index << endl; continue; }
      
      //double dR = deltaR(zjet.JetCalEta, zjet.JetCalPhi, zjet.TowEta[n], zjet.TowPhi[n]);
      
      double relativEt = zjet.TowEt[n]/zjet.JetCalEt;  
      //if (relativEt<=0) cerr << "relEt = " <<relativEt << endl; //continue;
      //This relativeE is used *only* for plotting! Therefore no cuts on this var!
      //create array with multidimensional measurement
      TMeasurement * mess = new TTower;
      mess->pt = double(zjet.TowEt[n]);
      double scale = zjet.TowEt[n]/zjet.TowE[n];
      mess->EMF = double(zjet.TowEm[n]*scale);
      mess->HadF = double(zjet.TowHad[n]*scale);
      mess->OutF = double(zjet.TowOE[n]*scale);
      mess->eta = double(zjet.TowEta[n]);
      mess->phi = double(zjet.TowPhi[n]);
      mess->E = double(zjet.TowE[n]);
      //mess[7] = double( cos( zjet.JetCalPhi-zjet.TowPhi[n] ) ); // Projection factor for summing tower Pt
      gj_data->AddMess(new TData_TruthMess(index,
					   mess,                                           //mess//
					   zjet.ZEt * relativEt,                           //truth//
					   sqrt(pow(0.5,2)+pow(0.1*zjet.ZEt*relativEt,2)), //error//
					   1.,                                             //weight//
					   p->GetTowerParRef( index ),                     //parameter//
					   p->GetNumberOfTowerParametersPerBin(),          //number of free tower param. p. bin//
					   p->tower_parametrization,                       //function//
					   tower_error_param                                        //error param.func.//
					   ));
    } 

    
    //Add the jet's tracks to "gj_data":
    for (int n=0; n<zjet.NobjTrack; ++n){
      
      int index = p->GetTrackBin(p->GetTrackEtaBin(zjet.TrackTowIdEta[n]),
                            p->GetTrackPhiBin(zjet.TrackTowIdPhi[n]));
      if (index<0){ cerr<<"WARNING: track_index = " << index << endl; continue; }
      //create array with multidimensional measurement
      //TMeasurement * Tmess = new TTrack;
      TTrack * Tmess = new TTrack;
      Tmess->TrackId = int(zjet.TrackId[n]);
      Tmess->TowerId = int(zjet.TrackTowId[n]);
      Tmess->pt = double(zjet.TrackPt[n]);
      double scale = zjet.TrackP[n]/zjet.TrackPt[n];
      Tmess->EM1 = double(zjet.TrackEMC1[n]*scale);
      Tmess->EMF = double(zjet.TrackEMC3[n]*scale);
      Tmess->EM5 = double(zjet.TrackEMC5[n]*scale);
      Tmess->Had1 = double(zjet.TrackHAC1[n]*scale);
      Tmess->HadF = double(zjet.TrackHAC3[n]*scale);
      Tmess->Had5 = double(zjet.TrackHAC5[n]*scale);
      Tmess->OutF = 0;
      Tmess->DR = double(zjet.TrackDR[n]);
      Tmess->DRout = double(zjet.TrackDROut[n]);
      Tmess->eta = double(zjet.TrackEta[n]);
      Tmess->etaOut = double(zjet.TrackEtaOut[n]);
      Tmess->phi = double(zjet.TrackPhi[n]);
      Tmess->phiOut = double(zjet.TrackPhiOut[n]);
      Tmess->E = double(zjet.TrackP[n]);
      Tmess->TrackChi2 = double(zjet.TrackChi2[n]);
      Tmess->NValidHits = int(zjet.TrackNHits[n]);
      Tmess->MuDR = double(zjet.MuDR[n]);
      Tmess->MuDE = double(zjet.MuDE[n]);
      //mess[7] = double( cos( zjet.JetCalPhi-zjet.TowPhi[n] ) ); // Projection factor for summing tower Pt
      //EM+=mess->EMF;
      //F+=mess->pt;
      gj_data->AddTrack(new TData_TruthMess(index,
					   Tmess,                                                    //mess//
					   0,                           //truth//
					   0.05 + 0.00015 * zjet.TrackPt[n], //error//
					   1.,                                                      //weight//
					   p->GetTowerParRef( index ),                              //parameter//
					   p->GetNumberOfTrackParametersPerBin(),                   //number of free tower param. p. bin//
					   p->track_parametrization,                                //function//
					   track_error_param                                        //error param.func.//
					   ));
    } 
    gj_data->UseTracks(useTracks);   //check if track information is sufficient to use Track Parametrization
    data.push_back( gj_data ); 
   
    if (n_zjet_events>=0 && i>=n_zjet_events-1)
      break;
  }
}

void TCaliber::Run_TrackTower()
{
  //Run Track-Tower stuff
  int nevent = tracktower.fChain->GetEntries();
  int evt=0;
  for (int i=0;i<nevent;i++) {
    tracktower.GetEntry(i); 
    if (tracktower.NobjTowCal>10000 || tracktower.NobjTrackCal>10000)
      cerr<<"ERROR: Increase array sizes in TrackTowerSelector; NobjTowCal="
	  <<tracktower.NobjTowCal<<", NobjTrackCal="<<tracktower.NobjTrackCal<<"!"<<endl;

    for (int n=0; n<tracktower.NobjTowCal && n<tracktower.NobjTrackCal; ++n){
      if (tracktower.TrackEt[n]<Et_cut_on_track ||
	  tracktower.TowEt[n]<Et_cut_on_tower)
	continue;

      int index=p->GetBin(p->GetEtaBin(tracktower.TowId_eta[n]),p->GetPhiBin(tracktower.TowId_phi[n]));
      if (index<0) {
	cerr << "INDEX = "<< index << endl;
	continue;
      }
      //create array with multidimensional measurement
      TMeasurement * mess = new TTower;
      mess->pt = double(tracktower.TowEt[n]);
      double scale = tracktower.TowEt[n]/tracktower.TowE[n];
      mess->EMF = double(tracktower.TowEm[n]*scale);
      mess->HadF = double(tracktower.TowHad[n]*scale);
      mess->OutF = double(tracktower.TowOE[n]*scale);
      mess->eta = double(tracktower.TowEta[n]);
      mess->phi = double(tracktower.TowPhi[n]);
      mess->E = double(tracktower.TowE[n]);
      //mess[7] = double( cos( tracktower.JetCalPhi-tracktower.TowPhi[n] ) ); // Projection factor for summing tower Pt
      data.push_back(new TData_TruthMess(index,
					 mess,                                                //mess//
					 tracktower.TrackEt[n],                               //truth//
					 //tracktower.TrackEterr[n],                          //error//
					 sqrt(pow(0.5,2)+ pow(0.1*tracktower.TrackEt[n] ,2)), //error//
					 //tracktower.EventWeight,                            //weight//
					 1.,                                                  //weight//
					 p->GetTowerParRef( index ),                          //parameter//
					 p->GetNumberOfTowerParametersPerBin(),               //number of free tower param. p. bin//
					 p->tower_parametrization,                            //function//
 				         tower_error_param                             //function//
					 ) );
      if((evt++)%1000==0) cout<<"Track-Tower Event: "<<evt<<endl;
      break;//use only one track-tower per event! ->bug in the producer
    }  
 
    if (n_tracktower_events>=0 && evt>=n_tracktower_events)
      break;
  }
}


void TCaliber::Run_TrackCluster()
{
  //Run Track-cluster stuff
  int nevent = trackcluster.fChain->GetEntries();
  int evt=0;
  for (int i=0;i<nevent;i++) {
    trackcluster.GetEntry(i); 
    if (trackcluster.NobjTowCal>200)
      cerr<<"ERROR: Increase array sizes in TrackClusterSelector; NobjTowCal="
	  <<trackcluster.NobjTowCal<<"!"<<endl;
     
    //Calculate cluster energy (needed for plotting)
    double cluster_energy = 0.0;	   
    for (int n=0; n<trackcluster.NobjTowCal; ++n)
      cluster_energy += trackcluster.TowEt[n];

    if (trackcluster.TrackEt < Et_cut_on_track || cluster_energy < Et_cut_on_cluster) continue;
    TMeasurement* clusterp  = new TJet;
    clusterp->pt  = cluster_energy;
    clusterp->eta = 0;
    clusterp->phi = 0;
    clusterp->E   = 0;
    //Define Track-Cluster event	
    TData_TruthMultMess * tc = new TData_TruthMultMess( 0,
							trackcluster.TrackEt,  			           //truth//
							sqrt(pow(0.5,2)+pow(0.10*trackcluster.TrackEt,2)), //error//
							//trackcluster.EventWeight,                        //weight//
							1.,                                                //weight//
							0,                                                 //params
							0,                                                 //number of free jet param. p. bin
							p->dummy_parametrization,                          //function
			                                jet_error_param,                                  //error param. function
							clusterp);
    tc->SetType( TrackCluster );
    //Add the towers to the event
    for (int n=0; n<trackcluster.NobjTowCal; ++n){
      //if (trackcluster.TrackEt[n]<Et_cut_on_track)
      //   continue;

      int index=p->GetBin(p->GetEtaBin(trackcluster.TowId_eta[n]),p->GetPhiBin(trackcluster.TowId_phi[n]));
      if (index<0) {
	cerr << "INDEX = "<< index << endl;
	continue;
      }
      //create array with multidimensional measurement
      TMeasurement * mess = new TTower;
      mess->pt = double(trackcluster.TowEt[n]);
      double scale = trackcluster.TowEt[n]/trackcluster.TowE[n];
      mess->EMF = double(trackcluster.TowEm[n]*scale);
      mess->HadF = double(trackcluster.TowHad[n]*scale);
      mess->OutF = double(trackcluster.TowOE[n]*scale);
      mess->eta = double(trackcluster.TowEta[n]);
      mess->phi = double(trackcluster.TowPhi[n]);
      mess->E = double(trackcluster.TowE[n]);
      //mess[7] = double( cos( trackcluster.JetCalPhi-trackcluster.TowPhi[n] ) ); // Projection factor for summing tower Pt

      TData_TruthMess * tower = new TData_TruthMess(index,
						    mess,                                                      //mess//
						    trackcluster.TrackEt*trackcluster.TowEt[n]/cluster_energy, //"truth" for plotting only!//
						    //trackcluster.TrackEterr[n],                              //error//
						    sqrt(pow(0.5,2)+ pow(0.1*trackcluster.TrackEt ,2)),        //error//
						    1.,                                                        //weight//
						    p->GetTowerParRef( index ),                                //parameter//
						    p->GetNumberOfTowerParametersPerBin(),                     //number of free cluster param. p. bin//
						    p->tower_parametrization,                                  //function//
					            tower_error_param                                        //error param.func.//
						    );
      tc->AddMess( tower );
    } 
     
    //Save event
    data.push_back( tc ); 

    if((evt++)%1000==0) cout<<"Track-Cluster Event: "<<evt<<endl;
    if (n_trackcluster_events>=0 && evt>=n_trackcluster_events)
      break;
  }
}
void TCaliber::AddTowerConstraint()
{

  for(std::vector<TowerConstraint>::const_iterator ic = tower_constraints.begin() ;
      ic != tower_constraints.end() ; ++ic) {
    std::cout << "adding constraint for towers " << ic->mineta << " to " << ic->maxeta
	      << " for em Et=" << ic->emEt << " and had Et=" << ic->hadEt 
	      << " with weight w=" << ic->weight << "\n";
    //constrain average tower response
    int ntowers= (ic->maxeta - ic->mineta + 1) * 72;
    if((ic->maxeta  > 0) && (ic->mineta < 0)) ntowers -= 72;
    double etsum = ic->hadEt + ic->emEt;
    TMeasurement* constraintp  = new TMeasurement;
    constraintp->pt  = etsum;
    constraintp->eta = 0;
    constraintp->phi = 0;
    constraintp->E   = etsum;
   
    TData_TruthMultMess * tc = new TData_TruthMultMess(0,
						       etsum * ntowers, //truth
						       sqrt(pow(0.5,2)+ pow(0.1*etsum * ntowers,2)), //error
						       ic->weight, //weight
						       0, //params
						       0, //number of free jet param. p. bin
						       p->dummy_parametrization, // function
						       p->const_error<10000>, // function
						       constraintp);
    tc->SetType( typeTowerConstraint );
    //Add the towers to the event
    for(int ideta = ic->mineta ; ideta <= ic->maxeta  ; ++ideta) {
      if(ideta == 0) ideta = 1;
      for(int idphi = 1 ; idphi <= 72  ; ++idphi) {
	int index=p->GetBin(p->GetEtaBin(ideta),p->GetPhiBin(idphi));
	if (index<0) {
	  cerr << "INDEX = "<< index << endl;
	  continue;
	}
	//create array with multidimensional measurement
	TMeasurement * mess = new TMeasurement;
	mess->pt = etsum;
	mess->EMF = ic->emEt;
	mess->HadF = ic->hadEt;
	mess->OutF = 0;
	mess->eta = 0;
	mess->phi = 0;
	mess->E = etsum;
	//mess[7] = double( cos( trackcluster.JetCalPhi-trackcluster.TowPhi[n] ) ); // Projection factor for summing tower Pt
	TData_TruthMess *tower = new TData_TruthMess(index,
						     mess, //mess
						     etsum, //"truth" for plotting only
						     sqrt(pow(0.5,2)+ pow(0.1*etsum,2)), //error
						     1.0, //weight ???
						     p->GetTowerParRef(index), //parameter
						     p->GetNumberOfTowerParametersPerBin(), //number of free cluster param. p. bin
						     p->tower_parametrization, //function
						     p->const_error<10> //error param.
						     );
	tc->AddMess(tower);
      } 
    } 
    data.push_back(tc);
  } 
}

void TCaliber::AddParameterLimits()
{

  for(std::vector<ParameterLimit>::const_iterator pl = par_limits.begin() ;
      pl != par_limits.end() ; ++pl) {
    std::cout << "adding limit for parameter " << pl->index+1 << " min:" 
	      << pl->min << " max:" << pl->max << " k:" << pl->k << '\n';
    TMeasurement* limitp  = new TMeasurement;
    limitp->pt  = pl->min;//@@Add new TMeasurement derivative
    limitp->eta = pl->max;
    TData_ParLimit * parlim = new TData_ParLimit(pl->index,limitp,pl->k,p->GetPars() + pl->index,p->parameter_limit);
    data.push_back(parlim);
  }
}


void TCaliber::Run_Top()
{
  //Run Top stuff  
  int nevent = top.fChain->GetEntries();
  int evt=0;
  
  for (int i=0;i<nevent;i++) {
    if((i+1)%1000==0) 
      cout<<"Top Event: "<<i+1<<endl;
    top.fChain->GetEvent(i); 
    if (top.NobjTow>1000 || top.NobjJet>8) {
      cerr << "ERROR: Increase array sizes in topSelector; NobjTow="
	   << top.NobjTow<<", NobjBJet="<<top.NobjJet<<"!"<<endl;
      exit(10);
    }
    //--------------

    if(useMassConstraintW) {
      bool goodevent = false;
      TData_InvMass2 * top_data[2]; //two W-jets
      top_data[0] = 0;
      int nstoredjets = 0;
      for (unsigned int ij = 0; ij<3; ++ij){
	if(top.JetPt[ij] < Et_cut_on_jet) continue;
	if((TJet::Flavor)top.JetFlavor[ij] != TJet::uds) continue;
	
	//Find the jets eta & phi index using the nearest tower to jet axis:
	int jet_index=-1;
	double min_tower_dr = 10.0;
	double em = 0;
	double had = 0;
	double out = 0;
	TLorentzVector Ljet(0,0,0,0);

	Ljet.SetPtEtaPhiE(top.JetPt[ij],top.JetEta[ij],top.JetPhi[ij],top.JetE[ij]);
	for (int n=0; n<top.NobjTow; ++n){
	  if (top.Tow_jetidx[n]!=(int)ij) continue;//look for ij-jet's towers
	  em += top.TowEm[n];
	  had += top.TowHad[n];
	  out += top.TowOE[n];
	  TLorentzVector Ltower(0,0,0,0);
	  Ltower.SetPtEtaPhiE(top.TowEt[n],top.TowEta[n],top.TowPhi[n],top.TowE[n]);
	  double dr = Ltower.DeltaR(Ljet);
	  if (dr<min_tower_dr) {
	    jet_index = p->GetJetBin(p->GetJetEtaBin(top.TowId_eta[n]),
				     p->GetJetPhiBin(top.TowId_phi[n]));
	    min_tower_dr = dr;
	  }
	}
	if (jet_index<0){
	  cerr<<"WARNING: Top jet_index = " << jet_index << endl;
	  continue;
	}
	
	double * direction = new double[3];
	direction[0] = sin(top.JetPhi[ij]);
	direction[1] = cos(top.JetPhi[ij]);
	direction[2] = sqrt(top.JetE[ij]*top.JetE[ij] - top.JetEt[ij]*top.JetEt[ij]);
	TJet* jetp  = new TJet;
	jetp->pt  = top.JetEt[ij];
	jetp->eta = top.JetEta[ij];
	jetp->phi = top.JetPhi[ij];
	jetp->E   = top.JetE[ij];
	jetp->flavor = (TJet::Flavor)top.JetFlavor[ij];
	//the following is not quite correct, as this factor is different for all towers.
	//These values should be in the n-tupel as well
	double factor =  top.JetEt[ij] /  top.JetE[ij];
	jetp->HadF = had * factor;
	jetp->EMF = em * factor;
	jetp->OutF = out * factor;
	//Create an Top TData event
	top_data[nstoredjets] = new TData_InvMass2(
          jet_index * p->GetNumberOfJetParametersPerBin() + p->GetNumberOfTowerParameters(),
	  direction,                                     //p_T direction of this jet
	  massConstraint_W,                                         //truth//
	  sqrt(pow(0.5,2)+pow(0.10*top.JetPt[ij],2)),    //error//
	  top.Weight,                                    //weight//
	  //1.,                                          //weight//
	  p->GetJetParRef( jet_index ),                  //params
	  p->GetNumberOfJetParametersPerBin(),           //number of free jet param. p. bin
	  p->jet_parametrization,                        //function
	  //p->dummy_parametrization,
          jet_error_param,                               //error param. function
	  jetp                                           //jet momentum for plotting and scale
	  );
	
	//Add the jet's towers to "top_data":
	for (int n=0; n<top.NobjTow; ++n){
	  if (top.Tow_jetidx[n]!=(int)ij) continue;//look for ij-jet's towers
	  //if (top.TowEt[n]<0.01) continue;

	  int index = p->GetBin(p->GetEtaBin(top.TowId_eta[n]),
				p->GetPhiBin(top.TowId_phi[n]));
	  //std::cout << "jet:" << ij << "bin index:" << index << "\n";
	  if (index<0){ cerr<<"WARNING: Top tower_index = " << index << endl; continue; }
	  
	  double relativEt = top.TowEt[n]/top.JetEt[ij];
	  //if (relativEt<=0) cerr << "relEt = " <<relativEt << endl; //continue;
	  //This relativeE is used *only* for plotting! Therefore no cuts on this var!
	  //create array with multidimensional measurement
	  TMeasurement * mess = new TTower;
	  mess->pt = double(top.TowEt[n]);
	  double scale = 0.;
	  if (top.TowE[n]!=0.) scale = top.TowEt[n]/top.TowE[n];
	  mess->EMF = double(top.TowEm[n]*scale);
	  mess->HadF = double(top.TowHad[n]*scale);
	  mess->OutF = double(top.TowOE[n]*scale);
	  mess->eta = double(top.TowEta[n]);
	  mess->phi = double(top.TowPhi[n]);
	  mess->E = double(top.TowE[n]);
	  //mess[7] = double( cos( top.JetCalPhi-top.TowPhi[n] ) ); // Projection factor for summing tower Pt
	  
	  top_data[nstoredjets]->AddMess(new TData_TruthMess(
	     index,
	     mess,                                                   //mess//
	     top.JetPt[ij] * relativEt,                             //truth//
	     sqrt(pow(0.5,2)+pow(0.1*top.JetPt[ij]*relativEt,2)),   //error//
	     //1., //weight//
	     top.Weight,                                            //weight//
	     p->GetTowerParRef( index ),                             //parameter//
	     p->GetNumberOfTowerParametersPerBin(),                  //number of free tower param. p. bin//
	     p->tower_parametrization,                               //function//
	     tower_error_param                                       //error param. function//
	     ));
	}
	if(nstoredjets> 0) {
          top_data[0]->AddNewMultMess( top_data[nstoredjets] );
	}
	++nstoredjets;
	goodevent = (nstoredjets==2 ? true : false);
      }//loop over all top-jets

      if (goodevent) {
	data.push_back( top_data[0] );
      }
    }

    if(useMassConstraintTop) {
      bool goodevent = false;
      TData_InvMass2 * top_data[3]; //two W-jets and one b-jet
      top_data[0] = 0;
      int nstoredjets = 0;
      for (unsigned int ij = 0; ij<3; ++ij){
	if(top.JetPt[ij] < Et_cut_on_jet) continue;
	
	//Find the jets eta & phi index using the nearest tower to jet axis:
	int jet_index=-1;
	double min_tower_dr = 10.0;
	double em = 0;
	double had = 0;
	double out = 0;
	TLorentzVector Ljet(0,0,0,0);

	Ljet.SetPtEtaPhiE(top.JetPt[ij],top.JetEta[ij],top.JetPhi[ij],top.JetE[ij]);
	for (int n=0; n<top.NobjTow; ++n){
	  if (top.Tow_jetidx[n]!=(int)ij) continue;//look for ij-jet's towers
	  em += top.TowEm[n];
	  had += top.TowHad[n];
	  out += top.TowOE[n];
	  TLorentzVector Ltower(0,0,0,0);
	  Ltower.SetPtEtaPhiE(top.TowEt[n],top.TowEta[n],top.TowPhi[n],top.TowE[n]);
	  double dr = Ltower.DeltaR(Ljet);
	  if (dr<min_tower_dr) {
	    jet_index = p->GetJetBin(p->GetJetEtaBin(top.TowId_eta[n]),
				     p->GetJetPhiBin(top.TowId_phi[n]));
	    min_tower_dr = dr;
	  }
	}
	if (jet_index<0){
	  cerr<<"WARNING: Top jet_index = " << jet_index << endl;
	  continue;
	}
	
	double * direction = new double[3];
	direction[0] = sin(top.JetPhi[ij]);
	direction[1] = cos(top.JetPhi[ij]);
	direction[2] = sqrt(top.JetE[ij]*top.JetE[ij] - top.JetEt[ij]*top.JetEt[ij]);
	TJet* jetp  = new TJet;
	jetp->pt  = top.JetEt[ij];
	jetp->eta = top.JetEta[ij];
	jetp->phi = top.JetPhi[ij];
	jetp->E   = top.JetE[ij];
	jetp->flavor = (TJet::Flavor)top.JetFlavor[ij];
	//the following is not quite correct, as this factor is different for all towers.
	//These values should be in the n-tupel as well
	double factor =  top.JetEt[ij] /  top.JetE[ij];
	jetp->HadF = had * factor;
	jetp->EMF = em * factor;
	jetp->OutF = out * factor;
	//Create an Top TData event
	top_data[nstoredjets] = new TData_InvMass2(
          jet_index * p->GetNumberOfJetParametersPerBin() + p->GetNumberOfTowerParameters(),
	  direction,                                     //p_T direction of this jet
	  massConstraint_Top,                                         //truth//
	  sqrt(pow(0.5,2)+pow(0.10*top.JetPt[ij],2)),    //error//
	  top.Weight,                                    //weight//
	  //1.,                                          //weight//
	  p->GetJetParRef( jet_index ),                  //params
	  p->GetNumberOfJetParametersPerBin(),           //number of free jet param. p. bin
	  p->jet_parametrization,                        //function
	  //p->dummy_parametrization,
          jet_error_param,                               //error param. function
	  jetp                                           //jet momentum for plotting and scale
	  );
	
	//Add the jet's towers to "top_data":
	for (int n=0; n<top.NobjTow; ++n){
	  if (top.Tow_jetidx[n]!=(int)ij) continue;//look for ij-jet's towers
	  //if (top.TowEt[n]<0.01) continue;

	  int index = p->GetBin(p->GetEtaBin(top.TowId_eta[n]),
				p->GetPhiBin(top.TowId_phi[n]));
	  //std::cout << "jet:" << ij << "bin index:" << index << "\n";
	  if (index<0){ cerr<<"WARNING: Top tower_index = " << index << endl; continue; }
	  
	  double relativEt = top.TowEt[n]/top.JetEt[ij];
	  //if (relativEt<=0) cerr << "relEt = " <<relativEt << endl; //continue;
	  //This relativeE is used *only* for plotting! Therefore no cuts on this var!
	  //create array with multidimensional measurement
	  TMeasurement * mess = new TTower;
	  mess->pt = double(top.TowEt[n]);
	  double scale = 0.;
	  if (top.TowE[n]!=0.) scale = top.TowEt[n]/top.TowE[n];
	  mess->EMF = double(top.TowEm[n]*scale);
	  mess->HadF = double(top.TowHad[n]*scale);
	  mess->OutF = double(top.TowOE[n]*scale);
	  mess->eta = double(top.TowEta[n]);
	  mess->phi = double(top.TowPhi[n]);
	  mess->E = double(top.TowE[n]);
	  //mess[7] = double( cos( top.JetCalPhi-top.TowPhi[n] ) ); // Projection factor for summing tower Pt
	  
	  top_data[nstoredjets]->AddMess(new TData_TruthMess(
	     index,
	     mess,                                                   //mess//
	     top.JetPt[ij] * relativEt,                             //truth//
	     sqrt(pow(0.5,2)+pow(0.1*top.JetPt[ij]*relativEt,2)),   //error//
	     //1., //weight//
	     top.Weight,                                            //weight//
	     p->GetTowerParRef( index ),                             //parameter//
	     p->GetNumberOfTowerParametersPerBin(),                  //number of free tower param. p. bin//
	     p->tower_parametrization,                               //function//
	     tower_error_param                                       //error param. function//
	     ));
	}
	if(nstoredjets> 0) {
          top_data[0]->AddNewMultMess( top_data[nstoredjets] );
	}
	++nstoredjets;
	goodevent = (nstoredjets==3 ? true : false);
      }//loop over all top-jets

      if (goodevent) {
	data.push_back( top_data[0] );
      }
    }

    ++evt;
    if (evt>=n_top_events && n_top_events>=0)
      break;
  }
}


//--------------------------------------------------------------------------------------------
void TCaliber::Run()
{
  if (fit_method!=3){

    time_t start = time(0);
    
    if (n_tracktower_events!=0)       Run_TrackTower();
    if (n_trackcluster_events!=0)     Run_TrackCluster();
    if (n_zjet_events!=0)             Run_ZJet();
    if (n_top_events!=0)              Run_Top();
    if (flatten_spectra){
      FlattenSpectra();
      //BalanceSpectra();
    }  
    if (! tower_constraints.empty())  AddTowerConstraint();
    if (! par_limits.empty())         AddParameterLimits();

    if (fit_method==1) Run_Lvmini();

    time_t end = time(0);
    cout << "Done, fitted " << p->GetNumberOfParameters() << " parameters in " << difftime(end,start) << " sec." << endl;
  } 
  //Dummy Configuration: Nothing to be done, start-values are written to file
}

void TCaliber::Run_Lvmini()
{ 
  //int naux = 1000000, niter=1000, iret=0;
  int naux = 3000000, niter=1000, iret=0;
  //int mvec = 29;
  int mvec = 6;
  //int mvec = 2;
  
  int npar = p->GetNumberOfParameters();

  naux = lvmdim_(npar,mvec);
  cout<<"array of size "<<naux<<" needed."<<endl;

  double aux[naux], fsum = 0;

  double *temp_derivative1 = new double[npar];
  double *temp_derivative2 = new double[npar];
  double epsilon = 1.E-3;

  cout << "\nFitting " << npar << " parameters; \n";
  p->Print();
  cout << " with LVMINI.\n" << "Using " << data.size() << " total events and ";
  cout << nthreads << " threads.\n";
  
  ComputeThread *t[nthreads];
  for (int ithreads=0; ithreads<nthreads; ++ithreads){
    t[ithreads] = new ComputeThread(npar, p->GetPars(),temp_derivative1,temp_derivative2,epsilon);
  }

  float eps =float(1.E-2);//-4
  float wlf1=1.E-4;
  float wlf2=0.9;
  lvmeps_(eps,wlf1,wlf2);

  //Set errors per default to 0 //@@ doesn't seem to work...
  int error_index=2;
  error_index = lvmind_(error_index);
  p->FillErrors(aux+error_index);

  for( unsigned int loop = 0; loop < _residualScalingScheme.size() ; ++loop ) {
    cout<<"Updating Di-Jet Errors"<<endl;
    for(DataIter it = data.begin()  ; it < data.end() ; ++it) {
      (*it)->UpdateError();
    }

    // Setting function to scale residuals in chi2 calculation
    cout << loop+1 << flush;
    if(  loop+1 == 1  ) cout << "st" << flush;
    else if(  loop+1 == 2  ) cout << "nd" << flush;
    else if(  loop+1 == 3  ) cout << "rd" << flush;
    else cout << "th" << flush;
    cout << " of " << _residualScalingScheme.size() <<" iteration(s): " << flush;
    if(  _residualScalingScheme.at(loop) == 0  ) {
	TData::ScaleResidual = &TData::ScaleNone;	
	cout << "no scaling of residuals." << endl;

	cout << "Rejecting outliers " << flush;
	DataIter beg = partition(data.begin(), data.end(), OutlierRejection(OutlierChi2Cut));
	for(DataIter i = beg ; i != data.end() ; ++i) {
	  delete *i;
	}
	data.erase(beg,data.end());
	cout << "and using " << data.size() << " events." << endl;
      }
    else if(  _residualScalingScheme.at(loop) == 1  ) {
	TData::ScaleResidual = &TData::ScaleCauchy;	
	cout << "scaling of residuals with Cauchy-Function." << endl;
      }
    else if(  _residualScalingScheme.at(loop) == 2  ) {
	TData::ScaleResidual = &TData::ScaleHuber;	
	cout << "scaling of residuals with Huber-Function." << endl;
      }
    else {
      cerr << "ERROR: " << _residualScalingScheme.at(loop) << " is not a valid scheme for resdiual scaling! Breaking iteration!" << endl;
      break;
    }
    if(lvmdim_(npar,mvec) > naux)
      cout<<"Aux field too small. "<<lvmdim_(npar,mvec)<<" enntires needed."<<endl;
    if (npar>0) npar*=-1; //Show output
    //initialization
    lvmini_( npar, mvec, niter, aux);
    npar=abs(npar);

    int n = 0;

    for(DataIter it = data.begin()  ; it < data.end() ; ++it) {
      t[n]->AddData(*it);
      n++;
      if(n == nthreads) n = 0;
    }

    do {
      //set storage for temporary derivative storage to zero
      for (int param=0; param< npar ; ++param) {
	temp_derivative1[param]=0.0;
	temp_derivative2[param]=0.0;
      }  
      fsum = 0;
      for (int ithreads=0; ithreads<nthreads; ++ithreads) t[ithreads]->Start();
      for (int ithreads=0; ithreads<nthreads; ++ithreads){
	if(t[ithreads]->IsDone()) fsum += t[ithreads]->Chi2();
      }
      //fast derivative calculation:
      for( int param = 0 ; param < std::abs(npar) ; ++param ) {
	aux[param]           = temp_derivative1[param]/(2.0*epsilon);
	aux[param+abs(npar)] = temp_derivative2[param]/(epsilon*epsilon);
      }
	
      lvmfun_(p->GetPars(),fsum,iret,aux);
      //p->SetParameters(aux + par_index); 
      lvmprt_(2,aux,2); //print out
    }
    while (iret<0);
    lvmprt_(2,aux,2); //print out
    for (int ithreads=0; ithreads<nthreads; ++ithreads){
      t[ithreads]->ClearData();
    }  
    int par_index = 1;
    par_index = lvmind_(par_index);
    p->SetParameters(aux + par_index);
  }
  //Copy Parameter errors from aux array to the TParameter::e array
  error_index=2;
  error_index = lvmind_(error_index);
  p->SetErrors(aux+error_index); 
  p->SetFitChi2(fsum);
  
  for (int ithreads=0; ithreads<nthreads; ++ithreads){
    delete t[ithreads];
  }
  delete [] temp_derivative1;
  delete [] temp_derivative2;
}
//--------------------------------------------------------------------------------------------


void TCaliber::Done()
{
  // write calibration to cfi output file if ending is cfi
  bool cfi=false;
  bool txt=false;
  std::string fileName(GetOutputFile());
  if( fileName.find(".cfi")!=std::string::npos ){
    if( fileName.substr(fileName.find(".cfi")).compare(".cfi")==0 ){
      p->Write_CalibrationCfi( fileName.c_str() );
      cfi=true; // file has a real .cfi ending
    }
  }
  // write calibration to cfi output file if ending is txt
  if( fileName.find(".txt")!=std::string::npos ){
    if( fileName.substr(fileName.find(".txt")).compare(".txt")==0 ){
      p->Write_CalibrationTxt( fileName.c_str() );
      txt=true; // file has a real .txt ending
    }
  }
  // write calibration to cfi & txt output file if w/o ending
  if( !cfi && !txt ){
    p->Write_CalibrationCfi( (fileName+".cfi").c_str() );
    p->Write_CalibrationTxt( (fileName+".txt").c_str() );
  }

  // Do Plots
  if(plots) {
    cout << endl << "Writing control plots in .ps " << flush;
    if( plots->OutputFormatRoot() ) cout << "and .root " << flush;
    cout << "format:" << endl;

      if( makeControlPlotsTowers )
	{
	  cout << "Creating tower control plots... " << flush;
	  plots->MakeControlPlotsTowers();
	  cout << "ok" << endl;
	}
      if(n_gammajet_events!=0)
	{
	  if( makeControlPlotsGammaJet )
	    {
	      cout << "Creating more gamma jet control plots... " << flush;
	      plots->MakeControlPlotsGammaJet(mPlottedQuant);
	      cout << "ok" << endl;
	    }
	  if( makeControlPlotsGammaJet2 )
	    {
	      cout << "Creating gamma jet (tower bin) control plots... " << flush;
	      plots->MakeControlPlotsGammaJetPerTowerBin();
	      cout << "ok" << endl;

	      cout << "Creating gamma jet (jet bin) control plots... " << flush;
	      plots->MakeControlPlotsGammaJetPerJetBin();
	      cout << "ok" << endl;
	      
	      cout << "Creating even more gamma jet control plots (sigmas)... " << flush;
	      plots->MakeControlPlotsGammaJetSigmas();
	      cout << "ok" << endl;
	    }
      }
    if (n_dijet_events!=0)   
      {
	if( makeControlPlotsDiJet )
	  {
	    cout << "Creating di-jet control plots... " << flush;
	    plots->MakeControlPlotsDiJet();
	    cout << "ok" << endl;
	  }
      }
    if( makeControlPlotsParScan ) 
      {
	cout << "Creating parameter scan control plots... " << flush;
	plots->MakeControlPlotsParameterScan();
	cout << "ok" << endl;
      }
  }
  // Clean-up
  cout << endl << "Cleaning up... " << flush;
  delete plots; 
  for(DataIter i = data.begin() ; i != data.end() ; ++i) {
    delete *i;
  }
  data.clear();
  cout << "Done" << endl;
}


void TCaliber::Init()
{
  ConfigFile config(configfile.c_str() );

  p = TParameters::CreateParameters(configfile);

  if(config.read<bool>("create plots",1))
    {
      plots = new TControlPlots(&data, p, config.read<bool>("plot output format",0));
      makeControlPlotsGammaJet = config.read<bool>("create gamma jet plots",false);
      if( makeControlPlotsGammaJet )
	{ 
	  vector<std::string> tmpPlottedQuant = bag_of_string(config.read<std::string>( "gamma jet plotted quantities", ""));
	  for(std::vector<std::string>::const_iterator it = tmpPlottedQuant.begin(); it < tmpPlottedQuant.end(); it++)
	    {
	      mPlottedQuant.insert(*it);
	    }
	}  				 
      makeControlPlotsGammaJet2 = config.read<bool>("create more gamma jet plots",false);
      makeControlPlotsDiJet = config.read<bool>("create dijet plots",false);
      makeControlPlotsTowers = config.read<bool>("create tower plots",false);
      makeControlPlotsParScan = config.read<bool>("create parameter scan plots",false);
    }

  //initialize temp arrays for fast derivative calculation
  TAbstractData::total_n_pars     = p->GetNumberOfParameters();
  //--------------------------------------------------------------------------
  //read config file
  fit_method = config.read<int>("Fit method",1);
  nthreads = config.read<int>("Number of Threads",1);
  flatten_spectra = config.read<int>("Flatten Spectra",1); 
  vector<double> limits = bag_of<double>(config.read<string>( "Jet Parameter Limits",""));

  if(limits.size() % 4 == 0) {
    for(unsigned int i = 0 ; i < limits.size() ; i += 4) {
      int index = (int)limits[i];
      for(int j = p->GetNumberOfTowerParameters() + index; 
	  j <  p->GetNumberOfParameters() ; 
	  j += p->GetNumberOfJetParametersPerBin()) {
	par_limits.push_back(ParameterLimit(j,limits[i+1],limits[i+2],
					    limits[i+3]));
      }
    }
  } else if(limits.size() > 1) {
    std::cout << "wrong number of arguments for Jet Parameter Limits:" 
	      << limits.size() << '\n';
  }
  limits.clear();
  limits = bag_of<double>(config.read<string>( "Tower Parameter Limits",""));

  if(limits.size() % 4 == 0) {
    for(unsigned int i = 0 ; i < limits.size() ; i += 4) {
      int index = (int)limits[i];
      for(int j = index; 
	  j <  p->GetNumberOfTowerParameters() ; 
	  j += p->GetNumberOfTowerParametersPerBin()) {
	par_limits.push_back(ParameterLimit(j,limits[i+1],limits[i+2],
					    limits[i+3]));
      }
    }
  } else if(limits.size() > 1) {
    std::cout << "wrong number of arguments for Tower Parameter Limits:" 
	      << limits.size() << '\n';
  }
  
  //Error Parametrization...
  //...for tracks:
  track_error_param = p->track_error_parametrization;
  //...for tower:
  string te = config.read<string>("tower error parametrization","standard"); 
  if (te=="standard")
    tower_error_param = p->tower_error_parametrization;
  else if (te=="fast")
    tower_error_param = p->fast_error_parametrization;
  else if (te=="Jans E parametrization")
    tower_error_param = p->jans_E_tower_error_parametrization;
  else if(te=="const")
    tower_error_param = p->const_error_parametrization;
  else if(te=="toy")
    tower_error_param = p->toy_tower_error_parametrization;
  else if(te=="jet")
    tower_error_param = p->jet_only_tower_error_parametrization;
  else  
    tower_error_param = p->tower_error_parametrization;
  //...for jets:
  string je = config.read<string>("jet error parametrization","standard");
  if (je=="standard")
    jet_error_param   = p->jet_error_parametrization;
  else if (je=="fast")
    jet_error_param   = p->fast_error_parametrization;
  else if (je=="dummy")
    jet_error_param   = p->dummy_error_parametrization;
  else if(je=="const")
    jet_error_param = p->const_error_parametrization;
  else if(je=="toy")
    jet_error_param   = p->toy_jet_error_parametrization;
  else if(je=="jet et")
    jet_error_param   = p->jet_only_jet_error_parametrization_et;
  else if(je=="jet energy")
    jet_error_param   = p->jet_only_jet_error_parametrization_energy;

  else  
    jet_error_param   = p->jet_error_parametrization;

  //last minute kinematic cuts
  //Et_cut_on_jet   = config.read<double>("Et cut on jet",0.0); 
  //Et_cut_on_gamma = config.read<double>("Et cut on gamma",0.0); 
  Et_cut_on_Z     = config.read<double>("Et cut on Z",0.0); 
  Et_cut_on_track = config.read<double>("Et cut on track",0.0); 
  Et_cut_on_tower = config.read<double>("Et cut on tower",0.0);
  Et_cut_on_cluster = config.read<double>("Et cut on cluster",0.0);
  //Et_cut_nplus1Jet = config.read<double>("Et cut on n+1 Jet",10.0);
  //Rel_cut_on_nJet  =  config.read<double>("Relative n+1 Jet Et Cut",0.2);
  //Rel_cut_on_gamma =  config.read<double>("Relative Rest Jet Cut",0.2);
  //relative sample weight 
  RelWeight[2]        = config.read<double>("Gamma-Jet weight",1.0);
  RelWeight[1]        = config.read<double>("Track-Tower weight",1.0);
  RelWeight[3]        = config.read<double>("Track-Cluster weight",1.0);
  RelWeight[4]        = config.read<double>("Di-Jet weight",1.0);
  RelWeight[5]        = config.read<double>("Multi-Jet weight",1.0);
  //RelWeight[2]        = config.read<double>("Z-Jet weight",1.0);

  //top quark specific settings
  useMassConstraintW   = config.read<bool>("use W mass constraint",true);
  useMassConstraintTop = config.read<bool>("use Top mass constraint",false);
  massConstraint_W   = config.read<double>("W mass",80.4);
  massConstraint_Top = config.read<double>("Top mass",172.4);
  if(!useMassConstraintW && !useMassConstraintTop) {
    std::cout << "W mass constraint and Top mass constraint were both turned off," << std::endl
	      << "you should enable at least one of these constraints if you want to run with Top...!" << std::endl;
  }

  //specify constraints
  vector<double> tower_constraint = bag_of<double>(config.read<string>( "Tower Constraint",""));
  if(tower_constraint.size() % 5 == 0) {
    for(unsigned int i = 0 ; i < tower_constraint.size() ; i += 5) {
      tower_constraints.push_back(TowerConstraint((int)tower_constraint[i],(int)tower_constraint[i+1],
						  tower_constraint[i+2],tower_constraint[i+3],
						  tower_constraint[i+4]));
    } 
  } else if(tower_constraint.size() > 1) {
    std::cout << "wrong number of arguments for tower constraint:" << tower_constraint.size() << '\n';
  }

  // Residual scaling
  const char* resScheme = ( config.read<string>("Residual Scaling Scheme","221").c_str() );
  while(  *resScheme != 0  )
    {
      int scheme = static_cast<int>(*resScheme - '0');
      if(  scheme < 0  ||  scheme > 2  )
	{
	  cerr << "ERROR: " << scheme << " is not a valid scheme for resdiual scaling! Using default scheme 111." << endl << endl;
	  _residualScalingScheme.clear();
	  _residualScalingScheme.push_back(1);
	  _residualScalingScheme.push_back(1);
	  _residualScalingScheme.push_back(1);
	  break;
	}

      _residualScalingScheme.push_back( static_cast<int>(*resScheme - '0') );
      resScheme++;
    }
  OutlierChi2Cut        = config.read<double>("Outlier Cut on Chi2",100.0);

  //input/output
  //n_gammajet_events     = config.read<int>("use Gamma-Jet events",-1);
  n_zjet_events         = config.read<int>("use Z-Jet events",-1);
  n_tracktower_events   = config.read<int>("use Track-Tower events",-1);
  n_trackcluster_events = config.read<int>("use Track-Cluster events",-1);
  //n_dijet_events        = config.read<int>("use Di-Jet events",-1);
  //n_trijet_events       = config.read<int>("use Tri-Jet events",-1);
  n_top_events          = config.read<int>("use Top events",-1);
  string default_tree_name = config.read<string>( "Default Tree Name","CalibTree");
  output_file = config.read<string>( "Output file", "calibration_k.cfi" );

  useTracks = config.read<bool>("use Tracks",true);
  if(p->GetNumberOfTrackParameters() < 1) useTracks = false;
  if(useTracks) cout<<"Tracks are used to calibrate jets"<<endl;
  else cout<<"Only Calorimeter information is used"<<endl;


  PhotonJetReader pjr(configfile,p);
  n_gammajet_events = pjr.readEvents(data);

  //Read Track-Tower Tree:
  string treename_tracktower = config.read<string>( "Track-Tower tree", default_tree_name );
  TChain * tchain_tracktower = new TChain( treename_tracktower.c_str() );
  vector<string> input_tracktower = bag_of_string( 
						  config.read<string>( "Track-Tower input file", "input/tracktower.root" ) );
  for (bag_of_string::const_iterator it = input_tracktower.begin(); it!=input_tracktower.end(); ++it){
    cout << "...opening root-file " << (*it) << " for Track-Tower analysis." << endl;
    tchain_tracktower->Add( it->c_str() );
  }  
  tracktower.Init( tchain_tracktower );
  
  //Read Track-Cluster Tree:
  string treename_trackcluster    = config.read<string>( "Track-Cluster tree", default_tree_name );
  TChain * tchain_trackcluster = new TChain( treename_trackcluster.c_str() );
  vector<string> input_trackcluster = bag_of_string( 
						    config.read<string>( "Track-Cluster input file", "input/trackcluster.root" ) );
  for (bag_of_string::const_iterator it = input_trackcluster.begin(); it!=input_trackcluster.end(); ++it){
    cout << "...opening root-file " << (*it) << " for Track-Cluster analysis." << endl;
    tchain_trackcluster->Add( it->c_str() );
  }  
  trackcluster.Init( tchain_trackcluster );
  
  DiJetReader djr(configfile,p);
  n_dijet_events = djr.readEvents(data);

  TriJetReader tjr(configfile,p);
  tjr.readEvents(data);

  //Read Z-Jet Tree:
  string treename_zjet      = config.read<string>( "Z-Jet tree", default_tree_name );
  TChain * tchain_zjet      = new TChain( treename_zjet.c_str() );
  vector<string> input_zjet = bag_of_string( config.read<string>( "Z-Jet input file", "input/zjet.root" ) );
  for (bag_of_string::const_iterator it = input_zjet.begin(); it!=input_zjet.end(); ++it){
    cout << "...opening root-file " << (*it) << " for Z-Jet analysis." << endl;
    tchain_zjet->Add( it->c_str() );
  }  
  zjet.Init( tchain_zjet );

  //Read Top Tree:
  string treename_top    = config.read<string>( "Top tree", default_tree_name );
  TChain * tchain_top = new TChain( treename_top.c_str() );
  vector<string> input_top = bag_of_string( config.read<string>( "Top input file", "input/top.root" ) );
  for (bag_of_string::const_iterator it = input_top.begin(); it!=input_top.end(); ++it){
    cout << "...opening root-file " << (*it) << " for Top analysis." << endl;
    tchain_top->Add( it->c_str() );
  }  
  top.Init( tchain_top );

}

//--^-TCaliber class-^------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------

int caliber(int argc, char *argv[])
{
  std::cout << "The University Hamburg Calorimeter Calibration Tool, 2007/08/15." << std::endl;
  
  TCaliber * Calibration;
  if (argc>1)
    Calibration = new TCaliber( argv[1] );
  else  
    Calibration = new TCaliber("config/calibration.cfg"); //Read input defined in config file
  
  Calibration->Init();
  Calibration->Run();  //Run Fit
  Calibration->Done(); //Do Plots & Write Calibration to file
  
  delete Calibration;    

  return 0;
}

void PrintUsage()
{
  std::cerr << "ERROR: You did something wrong! Better fix it." << std::endl;
}

int main(int argc, char *argv[])
{
  if (argc>2) {
    PrintUsage();
    exit(EXIT_FAILURE);
  }
  return caliber(argc, argv);
}

