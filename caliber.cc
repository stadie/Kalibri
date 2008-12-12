//
// Original Author:  Christian Autermann
//         Created:  Wed Jul 18 13:54:50 CEST 2007
// $Id: caliber.cc,v 1.69 2008/12/12 17:06:00 stadie Exp $
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
#include <ctime>

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
#include "ZJetReader.h"
#include "TopReader.h"
#include "TrackClusterReader.h"
#include "ParameterLimitsReader.h"
#include "TowerConstraintsReader.h"

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

//--------------------------------------------------------------------------------------------
void TCaliber::Run()
{
  if (fit_method!=3){

    time_t start = time(0);
    
    if (flatten_spectra){
      FlattenSpectra();
      //BalanceSpectra();
    }  
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
    plots->MakePlots();
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
      plots = new TControlPlots(configfile,&data, p);
    }

  //initialize temp arrays for fast derivative calculation
  TAbstractData::total_n_pars     = p->GetNumberOfParameters();
  //--------------------------------------------------------------------------
  //read config file
  fit_method = config.read<int>("Fit method",1);
  nthreads = config.read<int>("Number of Threads",1);
  flatten_spectra = config.read<int>("Flatten Spectra",1); 

 
  //last minute kinematic cuts
  Et_cut_on_jet   = config.read<double>("Et cut on jet",0.0); 
  Et_cut_on_gamma = config.read<double>("Et cut on gamma",0.0); 
  //relative sample weight 
  RelWeight[2]        = config.read<double>("Gamma-Jet weight",1.0);
  RelWeight[1]        = config.read<double>("Track-Tower weight",1.0);
  RelWeight[3]        = config.read<double>("Track-Cluster weight",1.0);
  RelWeight[4]        = config.read<double>("Di-Jet weight",1.0);
  RelWeight[5]        = config.read<double>("Multi-Jet weight",1.0);
  //RelWeight[2]        = config.read<double>("Z-Jet weight",1.0);

  // Residual scaling
  const char* resScheme = ( config.read<string>("Residual Scaling Scheme","221").c_str() );
  while(  *resScheme != 0  )
    {
      int scheme = static_cast<int>(*resScheme - '0');
      if(  scheme < 0  ||  scheme > 2  )
	{
	  cerr << "ERROR: " << scheme << " is not a valid scheme for resdiual scaling! Using default scheme 221." << endl << endl;
	  _residualScalingScheme.clear();
	  _residualScalingScheme.push_back(2);
	  _residualScalingScheme.push_back(2);
	  _residualScalingScheme.push_back(1);
	  break;
	}

      _residualScalingScheme.push_back( static_cast<int>(*resScheme - '0') );
      resScheme++;
    }
  OutlierChi2Cut        = config.read<double>("Outlier Cut on Chi2",100.0);

 
  output_file = config.read<string>( "Output file", "calibration_k.cfi" );

  //fill data vector
  PhotonJetReader pjr(configfile,p);
  n_gammajet_events = pjr.readEvents(data);
  
  DiJetReader djr(configfile,p);
  n_dijet_events = djr.readEvents(data);

  TriJetReader tjr(configfile,p);
  n_trijet_events = tjr.readEvents(data);

  ZJetReader zjr(configfile,p);
  n_zjet_events = zjr.readEvents(data);

  TopReader tr(configfile,p);
  n_top_events = tr.readEvents(data);
  
  TrackClusterReader tcr(configfile,p);
  n_trackcluster_events = tcr.readEvents(data);

  ParameterLimitsReader plr(configfile,p);
  plr.readEvents(data);

  TowerConstraintsReader cr(configfile,p);
  cr.readEvents(data);
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

