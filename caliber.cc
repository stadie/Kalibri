//
// Original Author:  Christian Autermann
//         Created:  Wed Jul 18 13:54:50 CEST 2007
// $Id: caliber.C,v 1.17 2008/04/18 10:26:18 auterman Exp $
//
#include "caliber.h"

//C++ libs
#include <cmath>
#include <iomanip>
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
#include "CalibData.h"
#include "external.h"

using namespace std;

typedef std::vector<TData*>::iterator DataIter;
//Outlier Rejection
struct OutlierRejection {
  OutlierRejection(double cut):_cut(cut){};
  bool operator()(TData *d){return (d->chi2()/d->GetWeight())<_cut;}
  double _cut;
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
      for (int param=0; param< parent->npar ; ++param) {
	parent->td1[param]= 0.0;
	parent->td2[param]= 0.0;
	parent->mypar[param] = parent->parorig[param];
      }
      parent->chi2 =0.0;   
      for (DataIter it=parent->data.begin() ; it!= parent->data.end() ; ++it) {
	parent->chi2 += (*it)->chi2_fast(parent->td1,parent->td2,parent->epsilon); 
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
  ComputeThread(int npar,double *par, double *temp_derivative1, double *temp_derivative2,
                            double epsilon) : npar(npar), td1(new double[npar]),
	            td2(new double[npar]), parorig(par),mypar(new double[npar]),
	            temp_derivative1(temp_derivative1), temp_derivative2(temp_derivative2), 
                            epsilon(epsilon) {}
  ~ComputeThread() {
    ClearData();
    delete [] td1;
    delete [] td2;
    delete [] mypar;
  }
  void AddData(TData* d) { 
    d->ChangeParAddress(parorig, mypar);
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

void TCaliber::global_fit(int &npar, double *gin, double &f, double *allpar, int iflag) 
//usage inadvisable -> derivative will be cached if chi2_fast() is used.
//This function can be used for fitting with Minuit which does derivatives 
//first or for plotting.
{
  double chisq = 0.0;
  std::vector<TData*>::const_iterator data_it, it;
  for (data_it=data.begin(); data_it!=data.end(); ++data_it)
    chisq += (*data_it)->GetWeight()*(*data_it)->chi2();     //standard
  f = chisq;
}

double TCaliber::numeric_derivate( void (*func)(int&,double*,double&,double*,int), 
                                   double * pars, int npar, int index)
//usage inadvisable -> derivative will be cached if chi2_fast() is used
{
  //double const epsilon = 1.E-5;//defined in CalibData.h as global static
  double * gin=0;
  int i=0;
  double x, x0, x1, result=0.0;
  double epsilon = 1e-03;
  if (index<npar){//first derivative
    pars[index]+=epsilon;
    func(npar,gin,x1,pars,i);
    pars[index]-=2*epsilon;
    func(npar,gin,x0,pars,i);
    pars[index]+=epsilon;
    result = (x1-x0)/(2.0*epsilon);
  } else if (npar<=index && index<2*npar){//second derivative
    pars[index-npar]+=epsilon;
    func(npar,gin,x1,pars,i);
    pars[index-npar]-=2*epsilon;
    func(npar,gin,x0,pars,i);
    pars[index-npar]+=epsilon;
    func(npar,gin,x,pars,i);
    result = (x0+x1-2*x)/(epsilon*epsilon);
  } 
  return result;
}


void TCaliber::Run_GammaJet()
//calculates from photon energy a truth value for one calo tower of the jet.
{
  //Run Gamma-Jet stuff  
  int nevent = gammajet.fChain->GetEntries();
  for (int i=0;i<nevent;i++) {
    if(i%1000==0) cout<<"Gamma-Jet Event: "<<i<<endl;
    gammajet.fChain->GetEvent(i); 
    if (gammajet.NobjTowCal>200) {
      cerr<<"ERROR: Increase array sizes in GammaJetSelector; NobjTowCal="
	  <<gammajet.NobjTowCal<<"!"<<endl;
      exit(8);
    }
 
    //trivial cuts
    if (gammajet.PhotonPt<Et_cut_on_gamma || gammajet.JetCalPt<Et_cut_on_jet) continue;
     
    //Find the jets eta & phi index using the leading (ET) tower:
    int jet_index=0;
    double max_tower_et = 0.0;
    double em = 0;
    for (int n=0; n<gammajet.NobjTowCal; ++n){
      em += gammajet.TowEm[n];
      if (gammajet.TowEt[n]>max_tower_et) {
	jet_index = p->GetJetBin(p->GetJetEtaBin(gammajet.TowId_eta[n]),
				 p->GetJetPhiBin(gammajet.TowId_phi[n]));
	max_tower_et = gammajet.TowEt[n];
      }
    }
    if (jet_index<0){ cerr<<"WARNING: jet_index = " << jet_index << endl; continue; }
    if(em == 0) { continue;}
    //jet_index: p->eta_granularity*p->phi_granularity*p->GetNumberOfTowerParametersPerBin()
    //           has to be added for a correct reference to k[...].
    double* jetp  = new double[3];
    jetp[0] = gammajet.JetCalEt;
    jetp[1] = gammajet.JetCalEta;
    jetp[2] = gammajet.JetCalPhi;
    //Create an Gamma/Jet TData event
    TData_TruthMultMess * gj_data = new 
      TData_TruthMultMess(jet_index + p->GetNumberOfTowerParameters(),
			  // gammajet.PhotonEt,				    //truth//
			  gammajet.JetGenPt,
			  sqrt(pow(0.5,2)+pow(0.10*gammajet.PhotonEt,2)),   //error//
			  //gammajet.EventWeight,                           //weight//
			  1.0,                                              //weight//
			  p->GetJetParRef( jet_index ),                     //params
			  p->GetNumberOfJetParametersPerBin(),              //number of free jet param. p. bin
			  p->jet_parametrization,                           //function
			  p->jet_error_parametrization,                     //function
			  jetp
			  );

    //Add the jet's towers to "gj_data":
    for (int n=0; n<gammajet.NobjTowCal; ++n){
      //if (gammajet.TowEt[n]<0.01) continue;
   
      int index = p->GetBin(p->GetEtaBin(gammajet.TowId_eta[n]),
			    p->GetPhiBin(gammajet.TowId_phi[n]));
      if (index<0){ cerr<<"WARNING: towewer_index = " << index << endl; continue; }

      //double dR = deltaR(gammajet.JetCalEta, gammajet.JetCalPhi, gammajet.TowEta[n], gammajet.TowPhi[n]);
	      
      double relativEt = gammajet.TowEt[n]/gammajet.JetCalEt;  
      //if (relativEt<=0) cerr << "relEt = " <<relativEt << endl; //continue;
      //This relativeE is used *only* for plotting! Therefore no cuts on this var!
      //create array with multidimensional measurement
      double * mess = new double[__DimensionMeasurement]; //__DimensionMeasurement difined in CalibData.h
      mess[0] = double(gammajet.TowEt[n]);
      double scale = gammajet.TowEt[n]/gammajet.TowE[n];
      mess[1] = double(gammajet.TowEm[n]*scale);
      mess[2] = double(gammajet.TowHad[n]*scale);
      mess[3] = double(gammajet.TowOE[n]*scale);
      gj_data->AddMess(new TData_TruthMess(index,
					   mess,                                                    //mess//
					   gammajet.PhotonEt * relativEt,                           //truth//
					   sqrt(pow(0.5,2)+pow(0.1*gammajet.PhotonEt*relativEt,2)), //error//
					   1.,                                                      //weight//
					   p->GetTowerParRef( index ),                              //parameter//
					   p->GetNumberOfTowerParametersPerBin(),                   //number of free tower param. p. bin//
					   p->tower_parametrization,                                //function//
					   p->tower_error_parametrization                           //function//
					   ));
    } 
 
    data.push_back( gj_data ); 
   
    if (n_gammajet_events>=0 && i==n_gammajet_events-1)
      break;
  }
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
     
    //Find the jets eta & phi index using the leading (ET) tower:
    int jet_index=0;
    double max_tower_et = 0.0;
    double em = 0;
    for (int n=0; n<zjet.NobjTowCal; ++n){
      em += zjet.TowEm[n];
      if (zjet.TowEt[n]>max_tower_et) {
	jet_index = p->GetJetBin(p->GetJetEtaBin(zjet.TowId_eta[n]),
				 p->GetJetPhiBin(zjet.TowId_phi[n]));
	max_tower_et = zjet.TowEt[n];
      }
    }
    if (jet_index<0){ cerr<<"WARNING: jet_index = " << jet_index << endl; continue; }
    if(em == 0) { continue;}
    //jet_index: p->eta_granularity*p->phi_granularity*p->GetNumberOfTowerParametersPerBin()
    //           has to be added for a correct reference to k[...].
    double* jetp  = new double[3];
    jetp[0] = zjet.JetCalEt;
    jetp[1] = zjet.JetCalEta;
    jetp[2] = zjet.JetCalPhi;
    //Create an Z/Jet TData event
    TData_TruthMultMess * gj_data = new 
      TData_TruthMultMess(jet_index + p->GetNumberOfTowerParameters(),
			  // zjet.ZEt,				    //truth//
			  zjet.JetGenPt,
			  sqrt(pow(0.5,2)+pow(0.10*zjet.ZEt,2)),   //error//
			  //zjet.EventWeight,                           //weight//
			  1.0,                                              //weight//
			  p->GetJetParRef( jet_index ),                     //params
			  p->GetNumberOfJetParametersPerBin(),              //number of free jet param. p. bin
			  p->jet_parametrization,                           //function
			  p->jet_error_parametrization,                     //function
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
      double * mess = new double[__DimensionMeasurement]; //__DimensionMeasurement difined in CalibData.h
      mess[0] = double(zjet.TowEt[n]);
      double scale = zjet.TowEt[n]/zjet.TowE[n];
      mess[1] = double(zjet.TowEm[n]*scale);
      mess[2] = double(zjet.TowHad[n]*scale);
      mess[3] = double(zjet.TowOE[n]*scale);
      gj_data->AddMess(new TData_TruthMess(index,
					   mess,                                                    //mess//
					   zjet.ZEt * relativEt,                           //truth//
					   sqrt(pow(0.5,2)+pow(0.1*zjet.ZEt*relativEt,2)), //error//
					   1.,                                                      //weight//
					   p->GetTowerParRef( index ),                              //parameter//
					   p->GetNumberOfTowerParametersPerBin(),                   //number of free tower param. p. bin//
					   p->tower_parametrization,                                //function//
					   p->tower_error_parametrization                           //function//
					   ));
    } 
 
    data.push_back( gj_data ); 
   
    if (n_zjet_events>=0 && i==n_zjet_events-1)
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
      double * mess = new double[__DimensionMeasurement]; //__DimensionMeasurement difined in CalibData.h
      mess[0] = double(tracktower.TowEt[n]);
      double scale = tracktower.TowEt[n]/tracktower.TowE[n];
      mess[1] = double(tracktower.TowEm[n])*scale;
      mess[2] = double(tracktower.TowHad[n]*scale);
      mess[3] = double(tracktower.TowOE[n])*scale;
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
					 p->tower_error_parametrization                       //function//
					 ) );
      if((evt++)%1000==0) cout<<"Track-Tower Event: "<<evt<<endl;
      break;//use only one track-tower per event! ->bug in the producer
    }  
 
    if (n_tracktower_events>=0 && evt==n_tracktower_events)
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
 
    //Define Track-Cluster event	
    TData_TruthMultMess * tc = new TData_TruthMultMess( 0,
							trackcluster.TrackEt,  			           //truth//
							sqrt(pow(0.5,2)+pow(0.10*trackcluster.TrackEt,2)), //error//
							//trackcluster.EventWeight,                        //weight//
							1.,                                                //weight//
							0,                                                 //params
							0,                                                 //number of free jet param. p. bin
							p->dummy_parametrization,                          //function
							p->jet_error_parametrization                       //function
							);
    tc->SetType( TypeTrackCluster );
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
      double * mess = new double[__DimensionMeasurement]; //__DimensionMeasurement difined in CalibData.h
      mess[0] = double(trackcluster.TowEt[n]);
      double scale = trackcluster.TowEt[n]/trackcluster.TowE[n];
      mess[1] = double(trackcluster.TowEm[n]*scale);
      mess[2] = double(trackcluster.TowHad[n]*scale);
      mess[3] = double(trackcluster.TowOE[n]*scale);
      TData_TruthMess * tower = new TData_TruthMess(index,
						    mess,                                                      //mess//
						    trackcluster.TrackEt*trackcluster.TowEt[n]/cluster_energy, //"truth" for plotting only!//
						    //trackcluster.TrackEterr[n],                              //error//
						    sqrt(pow(0.5,2)+ pow(0.1*trackcluster.TrackEt ,2)),        //error//
						    1.,                                                        //weight//
						    p->GetTowerParRef( index ),                                //parameter//
						    p->GetNumberOfTowerParametersPerBin(),                     //number of free cluster param. p. bin//
						    p->tower_parametrization,                                  //function//
						    p->tower_error_parametrization                             //function//
						    );
      tc->AddMess( tower );
    } 
     
    //Save event
    data.push_back( tc ); 

    if((evt++)%1000==0) cout<<"Track-Cluster Event: "<<evt<<endl;
    if (n_trackcluster_events>=0 && evt==n_trackcluster_events)
      break;
  }
}

void TCaliber::Run_NJet(NJetSel & njet)
{
  //@@ TODO: CHECK ERROR CALCULATION!!!
  
  
  //Run jet-Jet stuff  
  int nevent = njet.fChain->GetEntries();
  for (int i=0;i<nevent;i++) {
    if(i%1==0) cout<<"Jet-Jet Event: "<<i<<endl;
    njet.fChain->GetEvent(i); 
    if (njet.NobjTow>10000 || njet.NobjJet>100) {
      cerr << "ERROR: Increase array sizes in NJetSelector; NobjTow="
	   << njet.NobjTow<<", NobjJet="<<njet.NobjJet<<"!"<<endl;
      exit(9);
    }
    //--------------
    //  n - Jet
    //--------------
    TData_PtBalance * jj_data[ njet.NobjJet ];
    for (unsigned int ij = 0; (int)ij<njet.NobjJet; ++ij){
      //Find the jets eta & phi index using the leading (ET) tower:
      int jet_index = -1;
      double max_tower_et = 0.0;
      for (int n=0; n<njet.NobjTow; ++n){
        if (njet.Tow_jetidx[n]!=(int)ij) continue;//look for ij-jet's towers
//cout << ij<<". jet, "<<n<<". tow with ieta="<<njet.TowId_eta[n]
//     <<", iphi="<<njet.TowId_phi[n]
//     <<",  jetidx = " << p->GetJetBin(p->GetJetEtaBin(njet.TowId_eta[n]),
//				 p->GetJetPhiBin(njet.TowId_phi[n]))
//     <<endl;
	if (njet.TowEt[n]>max_tower_et) {
	  jet_index = p->GetJetBin(p->GetJetEtaBin(njet.TowId_eta[n]),
				   p->GetJetPhiBin(njet.TowId_phi[n]));
	  max_tower_et = njet.TowEt[n];
	}
      }
      if (jet_index<0){ 
	 cerr<<"WARNING: JJ jet_index = " << jet_index << endl; 
	 continue; 
      }

      double * direction = new double[2];
      direction[0] = sin(njet.JetPhi[ij]);
      direction[1] = cos(njet.JetPhi[ij]);

      //Create an jet/Jet TData event
      jj_data[ij] = new TData_PtBalance( 
          jet_index + p->GetNumberOfTowerParameters(),
	  direction,                                    //p_T direction of this jet
	  0.0,                                           //truth//
	  sqrt(pow(0.5,2)+pow(0.10*njet.JetPt[ij],2)), //error//
	  1.,                                            //weight//
	  p->GetJetParRef( jet_index ),                  //params
	  p->GetNumberOfJetParametersPerBin(),           //number of free jet param. p. bin
	  p->jet_parametrization,                        //function
	  p->jet_error_parametrization                   //function
        );

      //Add the jet's towers to "jj_data":
      for (int n=0; n<njet.NobjTow; ++n){
        if (njet.Tow_jetidx[n]!=(int)ij) continue;//look for ij-jet's towers
	//if (njet.TowEt[n]<0.01) continue;

	int index = p->GetBin(p->GetEtaBin(njet.TowId_eta[n]),
			      p->GetPhiBin(njet.TowId_phi[n]));
	if (index<0){ cerr<<"WARNING: JJ tower_index = " << index << endl; continue; }

	double relativEt = njet.TowEt[n]/njet.JetEt[ij];  
	//if (relativEt<=0) cerr << "relEt = " <<relativEt << endl; //continue;
	//This relativeE is used *only* for plotting! Therefore no cuts on this var!
	//create array with multidimensional measurement
	double * mess = new double[__DimensionMeasurement]; //__DimensionMeasurement difined in CalibData.h
	mess[0] = double(njet.TowEt[n]);
	double scale = njet.TowEt[n]/njet.TowE[n];
	mess[1] = double(njet.TowEm[n]*scale);
	mess[2] = double(njet.TowHad[n]*scale);
	mess[3] = double(njet.TowOE[n]*scale);
	jj_data[ij]->AddMess(new TData_TruthMess(
	    index,
	    mess,                                                   //mess//
	    njet.JetPt[ij] * relativEt,                           //truth//
	    sqrt(pow(0.5,2)+pow(0.1*njet.JetPt[ij]*relativEt,2)), //error//
            1.,                                                     //weight//
	    p->GetTowerParRef( index ),                             //parameter//
	    p->GetNumberOfTowerParametersPerBin(),                  //number of free tower param. p. bin//
	    p->tower_parametrization,                               //function//
	    p->tower_error_parametrization                          //function//
	  ));
	if (ij>0) 
  	  jj_data[0]->AddNewMultMess( jj_data[ij] );
      }  
    }//loop over all n-jets
    data.push_back( jj_data[0] ); 
  }
}



/*
void TCaliber::Run_NJet()
{
  //Run jet-Jet stuff  
  int nevent = njet.fChain->GetEntries();
  for (int i=0;i<nevent;i++) {
    if(i%1000==0) cout<<"Jet-Jet Event: "<<i<<endl;
    njet.fChain->GetEvent(i); 
    if (njet.NobjTowJ1Cal>200 || njet.NobjTowJ2Cal>200) {
      cerr << "ERROR: Increase array sizes in NJetSelector; NobjTowJ1Cal="
	   << njet.NobjTowJ1Cal<<", NobjTowJ2Cal="<<njet.NobjTowJ2Cal<<"!"<<endl;
      exit(9);
    }
    //--------------
    //   1. Jet
    //--------------
    //Find the jets eta & phi index using the leading (ET) tower:
    int jet_index = -1;
    double max_tower_et = 0.0;
    for (int n=0; n<njet.NobjTowJ1Cal; ++n){
      if (njet.TowJ1Et[n]>max_tower_et) {
	jet_index = p->GetJetBin(p->GetJetEtaBin(njet.TowJ1Id_eta[n]),
				 p->GetJetPhiBin(njet.TowJ1Id_phi[n]));
	max_tower_et = njet.TowJ1Et[n];
      }
    }
    if (jet_index<0){ 
       cerr<<"WARNING: JJ jet1_index = " << jet_index << endl; 
       continue; 
    }
     
    double * direction1 = new double[2];
    direction1[0] = sin(njet.FirstJetPhi);
    direction1[1] = cos(njet.FirstJetPhi);
     
    //Create an jet/Jet TData event
    TData_PtBalance * jj_data = new TData_PtBalance( jet_index + p->GetNumberOfTowerParameters(),
						     direction1,                                    //p_T direction of this jet
						     0.0,                                           //truth//
						     sqrt(pow(0.5,2)+pow(0.10*njet.ScndJetPt,2)), //error//
						     //njet.EventWeight,                          //weight//
						     1.,                                            //weight//
						     p->GetJetParRef( jet_index ),                  //params
						     p->GetNumberOfJetParametersPerBin(),           //number of free jet param. p. bin
						     p->jet_parametrization,                        //function
						     p->jet_error_parametrization                   //function
						     );

    //Add the jet's towers to "jj_data":
    double control_sum=0.0;
    for (int n=0; n<njet.NobjTowJ1Cal; ++n){
      //if (njet.TowJ1Et[n]<0.01) continue;

      int index = p->GetBin(p->GetEtaBin(njet.TowJ1Id_eta[n]),
			    p->GetPhiBin(njet.TowJ1Id_phi[n]));
      if (index<0){ cerr<<"WARNING: JJ towerj1_index = " << index << endl; continue; }

      double relativEt = njet.TowJ1Et[n]/njet.FirstJetEt;  
      //if (relativEt<=0) cerr << "relEt = " <<relativEt << endl; //continue;
      //This relativeE is used *only* for plotting! Therefore no cuts on this var!
      //create array with multidimensional measurement
      double * mess = new double[__DimensionMeasurement]; //__DimensionMeasurement difined in CalibData.h
      mess[0] = double(njet.TowJ1Et[n]);
      double scale = njet.TowJ1Et[n]/njet.TowJ1E[n];
      mess[1] = double(njet.TowJ1Em[n]*scale);
      mess[2] = double(njet.TowJ1Had[n]*scale);
      mess[3] = double(njet.TowJ1OE[n]*scale);
      jj_data->AddMess(new TData_TruthMess(index,
					   mess,                                                   //mess//
					   njet.ScndJetPt * relativEt,                           //truth//
					   sqrt(pow(0.5,2)+pow(0.1*njet.ScndJetPt*relativEt,2)), //error//
                                           1.,                                                     //weight//
					   p->GetTowerParRef( index ),                             //parameter//
					   p->GetNumberOfTowerParametersPerBin(),                  //number of free tower param. p. bin//
					   p->tower_parametrization,                               //function//
					   p->tower_error_parametrization                          //function//
					   ));
    }  

    //--------------
    //   2. Jet
    //--------------
    //Find the jets eta & phi index using the leading (ET) tower:
    max_tower_et = 0.0;
    for (int n=0; n<njet.NobjTowJ2Cal; ++n){
      if (njet.TowJ2Et[n]>max_tower_et) {
	jet_index = p->GetJetBin(p->GetJetEtaBin(njet.TowJ2Id_eta[n]),
				 p->GetJetPhiBin(njet.TowJ2Id_phi[n]));
	max_tower_et = njet.TowJ2Et[n];
      }
    }
    if (jet_index<0){ cerr<<"WARNING: JJ jet2_index = " << jet_index << endl; continue; }
     
    //jet_index: p->eta_granularity*p->phi_granularity*p->GetNumberOfTowerParametersPerBin()
    //           has to be added for a correct reference to k[...].

    double * direction2 = new double[2];
    direction2[0] = sin(njet.ScndJetPhi);
    direction2[1] = cos(njet.ScndJetPhi);
     
    //Create an jet/Jet TData event
    TData_PtBalance * jj2_data = new TData_PtBalance( jet_index + p->GetNumberOfTowerParameters(),
						      direction2,
						      0.0,				             //truth//
						      sqrt(pow(0.5,2)+pow(0.10*njet.ScndJetPt,2)), //error//
                                                      //njet.EventWeight,                          //weight//
                                                      1.,                                            //weight//
						      p->GetJetParRef( jet_index ),                  //params
						      p->GetNumberOfJetParametersPerBin(),           //number of free jet param. p. bin
						      p->jet_parametrization,                        //function
						      p->jet_error_parametrization                   //function
						      );

    //Add the jet's towers to "jj_data":
    control_sum=0.0;
    for (int n=0; n<njet.NobjTowJ2Cal; ++n){
      //if (njet.TowJ2Et[n]<0.01) continue;

      int index = p->GetBin(p->GetEtaBin(njet.TowJ2Id_eta[n]),
			    p->GetPhiBin(njet.TowJ2Id_phi[n]));
      if (index<0){ cerr<<"WARNING: JJ towerj1_index = " << index << endl; continue; }

      double relativEt = njet.TowJ2Et[n]/njet.ScndJetEt;  
      //if (relativEt<=0) cerr << "relEt = " <<relativEt << endl; //continue;
      //This relativeE is used *only* for plotting! Therefore no cuts on this var!
      //create array with multidimensional measurement
      double * mess = new double[4];
      mess[0] = double(njet.TowJ2Et[n]);
      double scale = njet.TowJ2Et[n]/njet.TowJ2E[n];
      mess[1] = double(njet.TowJ2Em[n]*scale);
      mess[2] = double(njet.TowJ2Had[n]*scale);
      mess[3] = double(njet.TowJ2OE[n]*scale);
      jj2_data->AddMess(new TData_TruthMess(index,
					    mess, //mess
					    njet.FirstJetPt * relativEt,                           //truth//
					    sqrt(pow(0.5,2)+pow(0.1*njet.FirstJetPt*relativEt,2)), //error//
					    1.,                                                      //weight//
					    p->GetTowerParRef( index ),                              //parameter//
					    p->GetNumberOfTowerParametersPerBin(),                   //number of free tower param. p. bin//
					    p->tower_parametrization,                                //function//
					    p->tower_error_parametrization                           //function//
					    ));
    }
    jj_data->AddNewMultMess( jj2_data );
    data.push_back( jj_data ); 

    if (n_njet_events>=0 && i==n_njet_events-1)
      break;
  }
}
*/

//--------------------------------------------------------------------------------------------
void TCaliber::Run()
{
  if (fit_method!=3){
    if (n_gammajet_events!=0)     Run_GammaJet();
    if (n_tracktower_events!=0)   Run_TrackTower();
    if (n_trackcluster_events!=0) Run_TrackCluster();
    if (n_dijet_events!=0)        Run_NJet( dijet);
    if (n_trijet_events!=0)       Run_NJet( trijet);
    if (n_zjet_events!=0)         Run_ZJet();

    if (fit_method==1) Run_Lvmini();
  } 
  //Dummy Configuration: Nothing to be done, start-values are written to file
}

void TCaliber::Run_Lvmini()
{
  int naux = 1000000, niter=1000, iret=0;
  //int mvec = 29;
  int mvec = 6;
  //int mvec = 2;
  double aux[naux], fsum = 0;

  int npar = p->GetNumberOfParameters();
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

  float eps =float(1.E-3*data.size());
  float wlf1=1.E-4;
  float wlf2=0.9;
  lvmeps_(eps,wlf1,wlf2);

  //Set errors per default to 0 //@@ doesn't seem to work...
  int error_index=2;
  error_index = lvmind_(error_index);
  p->FillErrors(aux+error_index);
  //for (int n=0; n<naux; ++n) aux[n]=0.0; 

  //outlier rejection before first iteration
  {
    DataIter beg =  partition(data.begin(), data.end(), OutlierRejection(OutlierChi2CutPresel));
    for(DataIter i = beg ; i != data.end() ; ++i) {
      delete *i;
    }
    data.erase(beg,data.end());  
  }
  for (int i=0; i<OutlierIterationSteps; ++i) {
    cout << i+1 << "th of "<<OutlierIterationSteps<<" iteration using " << data.size() << " events (chi2 terms)." << endl;
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
      for (int  ithreads=0; ithreads<nthreads; ++ithreads) t[ithreads]->Start();
      for (int ithreads=0; ithreads<nthreads; ++ithreads){
	if(t[ithreads]->IsDone()) fsum += t[ithreads]->Chi2();
      }
      
      //fast derivative calculation:
      for (unsigned param=0; param<abs(npar); ++param) {
	aux[param]           = (temp_derivative2[param]-temp_derivative1[param])/(2.0*epsilon);
	aux[param+abs(npar)] = (temp_derivative2[param]+temp_derivative1[param]-2*fsum)/(epsilon*epsilon);
      }
	
      lvmfun_(p->GetPars(),fsum,iret,aux);
      //p->SetParameters(aux + par_index); 
      lvmprt_(2,aux,2); //Has any effect?
    }
    while (iret<0);
    lvmprt_(2,aux,2); //Has any effect?
    //outlier rejection
    for (int ithreads=0; ithreads<nthreads; ++ithreads){
      t[ithreads]->ClearData();
    }  
    int par_index = 1;
    par_index = lvmind_(par_index);
    p->SetParameters(aux + par_index);
    if (i+1!=OutlierIterationSteps)  {
      std::vector<TData*>::iterator beg =  partition(data.begin(), data.end(), OutlierRejection(OutlierChi2Cut));
      for(std::vector<TData*>::iterator i = beg ; i != data.end() ; ++i) {
	delete *i;
      }
      data.erase(beg,data.end());
    }
  }
  //Copy Parameter errors from aux array to the TParameter::e array
  error_index=2;
  error_index = lvmind_(error_index);
  p->SetErrors(aux+error_index); 
  p->SetFitChi2(fsum);
  for (int ithreads=0; ithreads<nthreads; ++ithreads){
    delete t[ithreads];
  }
  delete []  temp_derivative1;
  delete []  temp_derivative2;
}
//--------------------------------------------------------------------------------------------


void TCaliber::Done()
{
  //Write calibration to file
  cout << "Writing calibration to file '"<<GetOutputFile()<<"',"<<endl;
  ofstream outfile (this->GetOutputFile(),ofstream::binary);
  outfile << (*p);
  outfile.close();
  
  //Do Plots
  if(plots) {
    cout << "Creating tower control plots,"<<endl;
    plots->FitControlPlots();
    cout << "Creating gamma jet (tower bin) control plots,"<<endl;
    plots->GammaJetControlPlots();
    cout << "Creating gamma jet (jet bin) control plots,"<<endl;
    plots->GammaJetControlPlotsJetBin();
    cout << "Creating more gamma jet control plots,"<<endl;
    plots->GammaJetControlPlotsJetJEC();
    cout << "Creating track tower control plots,"<<endl;
    plots->TrackTowerControlPlots();
    cout << "Creating track cluster control plots,"<<endl;
    plots->TrackClusterControlPlots();
  }
  //Clean-up
  delete plots; 
  for(DataIter i = data.begin() ; i != data.end() ; ++i) {
    delete *i;
  }
  data.clear();
  cout << "Done, cleaning up."<<endl;
}


void TCaliber::Init(string file)
{
  //set root style
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetTitleFillColor(0);  
  gStyle->SetTitleBorderSize(1);
  gStyle->SetStatColor(0);          
  gStyle->SetStatBorderSize(0);     
  gStyle->SetStatX(0.89);              
  gStyle->SetStatY(0.89);              
  gStyle->SetStatW(0.2);              
  gStyle->SetStatH(0.2);              


  
  ConfigFile config( file.c_str() );

  p = TParameters::CreateParameters(file);

  if(config.read<bool>("create plots",1)) {
    plots = new TControlPlots(file, &data, p);
  }
  //initialize temp arrays for fast derivative calculation
  TData::total_n_pars     = p->GetNumberOfParameters();
  //--------------------------------------------------------------------------
  //read config file
  fit_method = config.read<int>("Fit method",1);
  nthreads = config.read<int>("Number of Threads",1);
  //last minute kinematic cuts
  Et_cut_on_jet   = config.read<double>("Et cut on jet",0.0); 
  Et_cut_on_gamma = config.read<double>("Et cut on gamma",0.0); 
  Et_cut_on_Z     = config.read<double>("Et cut on Z",0.0); 
  Et_cut_on_track = config.read<double>("Et cut on track",0.0); 
  Et_cut_on_tower = config.read<double>("Et cut on tower",0.0);
  Et_cut_on_cluster = config.read<double>("Et cut on cluster",0.0);
  //outlier rejection
  OutlierIterationSteps = config.read<int>("Outlier Iteration Steps",3);
  OutlierChi2Cut        = config.read<double>("Outlier Cut on Chi2",100.0);
  OutlierChi2CutPresel  = config.read<double>("Outlier Cut on Chi2 presel",400.0);
  //input/output
  n_gammajet_events     = config.read<int>("use Gamma-Jet events",-1);
  n_zjet_events         = config.read<int>("use Z-Jet events",-1);
  n_tracktower_events   = config.read<int>("use Track-Tower events",-1);
  n_trackcluster_events = config.read<int>("use Track-Cluster events",-1);
  n_dijet_events        = config.read<int>("use Di-Jet events",-1);
  n_trijet_events       = config.read<int>("use Tri-Jet events",-1);
  string default_tree_name = config.read<string>( "Default Tree Name","CalibTree");
  output_file = config.read<string>( "Output file", "calibration_k.cfi" );

  //--------------------------------------------------------------------------
  //Read Gamma-Jet Tree:
  string treename_gammajet = config.read<string>( "Gamma-Jet tree", default_tree_name );
  TChain * tchain_gammajet = new TChain( treename_gammajet.c_str() );
  vector<string> input_gammajet = bag_of_string( 
						config.read<string>( "Gamma-Jet input file", "input/gammajet.root" ) );
  for (bag_of_string::const_iterator it = input_gammajet.begin(); it!=input_gammajet.end(); ++it){
    cout << "...opening root-file " << (*it) << " for Gamma-Jet analysis." << endl;
    tchain_gammajet->Add( it->c_str() );
  }  
  gammajet.Init( tchain_gammajet );

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
  
  //Read Di-Jet Tree:
  string treename_dijet    = config.read<string>( "Di-Jet tree", default_tree_name );
  TChain * tchain_dijet = new TChain( treename_dijet.c_str() );
  vector<string> input_dijet = bag_of_string( 
					      config.read<string>( "Di-Jet input file", "input/dijet.root" ) );
  for (bag_of_string::const_iterator it = input_dijet.begin(); it!=input_dijet.end(); ++it){
    cout << "...opening root-file " << (*it) << " for Di-Jet analysis." << endl;
    tchain_dijet->Add( it->c_str() );
  }  
  dijet.Init( tchain_dijet );

  //Read Tri-Jet Tree:
  string treename_trijet    = config.read<string>( "Tri-Jet tree", default_tree_name );
  TChain * tchain_trijet = new TChain( treename_trijet.c_str() );
  vector<string> input_trijet = bag_of_string( 
					      config.read<string>( "Tri-Jet input file", "input/trijet.root" ) );
  for (bag_of_string::const_iterator it = input_trijet.begin(); it!=input_trijet.end(); ++it){
    cout << "...opening root-file " << (*it) << " for Tri-Jet analysis." << endl;
    tchain_trijet->Add( it->c_str() );
  }  
  trijet.Init( tchain_trijet );

  //Read Z-Jet Tree:
  string treename_zjet      = config.read<string>( "Z-Jet tree", default_tree_name );
  TChain * tchain_zjet      = new TChain( treename_zjet.c_str() );
  vector<string> input_zjet = bag_of_string( 
					      config.read<string>( "Z-Jet input file", "input/zjet.root" ) );
  for (bag_of_string::const_iterator it = input_zjet.begin(); it!=input_zjet.end(); ++it){
    cout << "...opening root-file " << (*it) << " for Z-Jet analysis." << endl;
    tchain_zjet->Add( it->c_str() );
  }  
  zjet.Init( tchain_zjet );
}

//--^-TCaliber class-^------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------

int caliber(int argc, char *argv[])
{
  std::cout << "The University Hamburg Calorimeter Calibration Tool, 2007/08/15." << std::endl;
  
  TCaliber * Calibration = new TCaliber();
  if (argc>1)
    Calibration->Init( argv[1] );
  else  
    Calibration->Init("config/calibration.cfg"); //Read input defined in config file
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

