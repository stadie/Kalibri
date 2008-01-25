#include "caliber.h"

//C++ libs
#include <cmath>
#include <iomanip>
//Root libs
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TStyle.h"
// User
#include "CalibMath.h"
#include "CalibData.h"
#include "external.h"
using namespace std;

//Global defined function; necessary for TMinuit
void (*fitfunction)(int &npar, double *gin, double &f, double *allpar, int iflag);

//The data (needs to be global, since it's used in the static member function
//          global fit. Thats needed for TMinuit.)
std::vector<TData*> data;

//Outlier Rejection
struct OutlierRejection {
   OutlierRejection(double cut):_cut(cut){};
   bool operator()(TData *d){return d->chi2()<_cut;}
   double _cut;
};


//--------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------
//--v-TCaliber class-v------------------------------------------------------------------------

void TCaliber::global_fit_fast(int &npar, double *gin, double &f, double *allpar, int iflag) 
{
   //set storage for temporary derivative storage to zero
   for (unsigned param=0; param<abs(npar); ++param) {
     TData::temp_derivative1[param]=0.0;
     TData::temp_derivative2[param]=0.0;
   }

   double chisq = 0.0;
   std::vector<TData*>::const_iterator data_it, it;
   for (data_it=data.begin(); data_it!=data.end(); ++data_it){
     chisq += (*data_it)->chi2_fast();  //caches derivatives ->fast
   }
   f = chisq;
}

void TCaliber::global_fit(int &npar, double *gin, double &f, double *allpar, int iflag) 
//usage inadvisable -> derivative will be cached if chi2_fast() is used.
//This function can be used for fitting with Minuit which does derivatives 
//first or for plotting.
{
   double chisq = 0.0;
   std::vector<TData*>::const_iterator data_it, it;
   for (data_it=data.begin(); data_it!=data.end(); ++data_it)
     chisq += (*data_it)->chi2();     //standard
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

  if (index<npar){//first derivative
    pars[index]+=TData::epsilon;
    func(npar,gin,x1,pars,i);
    pars[index]-=2*TData::epsilon;
    func(npar,gin,x0,pars,i);
    pars[index]+=TData::epsilon;
    result = (x1-x0)/(2.0*TData::epsilon);
  } else if (npar<=index && index<2*npar){//second derivative
    pars[index-npar]+=TData::epsilon;
    func(npar,gin,x1,pars,i);
    pars[index-npar]-=2*TData::epsilon;
    func(npar,gin,x0,pars,i);
    pars[index-npar]+=TData::epsilon;
    func(npar,gin,x,pars,i);
    result = (x0+x1-2*x)/(TData::epsilon*TData::epsilon);
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
     int jet_index;
     double max_tower_et = 0.0;
     for (int n=0; n<gammajet.NobjTowCal; ++n){
       if (gammajet.TowEt[n]>max_tower_et) {
         jet_index = p->GetJetBin(p->GetJetEtaBin(gammajet.TowId_eta[n]),
	                          p->GetJetPhiBin(gammajet.TowId_phi[n]));
         max_tower_et = gammajet.TowEt[n];
       }
     }
     if (jet_index<0){ cerr<<"WARNING: jet_index = " << jet_index << endl; continue; }
     
     //jet_index: p->eta_granularity*p->phi_granularity*p->free_pars_per_bin
     //           has to be added for a correct reference to k[...].
     
     //Create an Gamma/Jet TData event
     TData_TruthMultMess * gj_data = new TData_TruthMultMess( jet_index + p->GetNumberOfTowerParameters(),
			    gammajet.PhotonEt,				    //truth//
			    sqrt(pow(0.5,2)+pow(0.10*gammajet.PhotonEt,2)), //error//
        		    p->GetJetParRef( jet_index ),                   //params
			    p->free_pars_per_bin_jet,                       //number of free jet param. p. bin
							      p->jet_parametrization,                            //function
							      p->jet_error_parametrization                       //function
     );

     //Add the jet's towers to "gj_data":
     double control_sum=0.0;
     for (int n=0; n<gammajet.NobjTowCal; ++n){
       //if (gammajet.TowEt[n]<0.01) continue;

       int index = p->GetBin(p->GetEtaBin(gammajet.TowId_eta[n]),
	                       p->GetPhiBin(gammajet.TowId_phi[n]));
       if (index<0){ cerr<<"WARNING: towewer_index = " << index << endl; continue; }

       double dR = deltaR(gammajet.JetCalEta, gammajet.JetCalPhi, 
                          gammajet.TowEta[n], gammajet.TowPhi[n]);
	      
       double relativEt = gammajet.TowEt[n]/gammajet.JetCalEt;  
       //if (relativEt<=0) cerr << "relEt = " <<relativEt << endl; //continue;
       //This relativeE is used *only* for plotting! Therefore no cuts on this var!
       //create array with multidimensional measurement
       double * mess = new double[4];
       mess[0] = double(gammajet.TowEt[n]);
       double scale = gammajet.TowEt[n]/gammajet.TowE[n];
       mess[1] = double(gammajet.TowEm[n]*scale);
       mess[2] = double(gammajet.TowHad[n]*scale);
       mess[3] = double(gammajet.TowOE[n]*scale);
       gj_data->AddMess(new TData_TruthMess(index,
	 mess,                                                    //mess//
 	 gammajet.PhotonEt * relativEt,                           //truth//
  	 sqrt(pow(0.5,2)+pow(0.1*gammajet.PhotonEt*relativEt,2)), //error//
	 p->GetTowerParRef( index ),                              //parameter//
         p->free_pars_per_bin,                                    //number of free tower param. p. bin//
         p->tower_parametrization,                                   //function//
	 p->tower_error_parametrization                              //function//
       ));
     }  
     data.push_back( gj_data ); 
   
     if (n_gammajet_events>=0 && i==n_gammajet_events-1)
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
       double * mess = new double[4];
       mess[0] = double(tracktower.TowEt[n]);
       double scale = tracktower.TowEt[n]/tracktower.TowE[n];
       mess[1] = double(tracktower.TowEm[n])*scale;
       mess[2] = double(tracktower.TowHad[n]*scale);
       mess[3] = double(tracktower.TowOE[n])*scale;
       data.push_back(new TData_TruthMess(index,
			 mess,                                                //mess//
			 tracktower.TrackEt[n],                               //truth//
		       //tracktower.TrackEterr[n],                            //error//
  		         sqrt(pow(0.5,2)+ pow(0.1*tracktower.TrackEt[n] ,2)), //error//
			 p->GetTowerParRef( index ),                          //parameter//
	                 p->free_pars_per_bin,                                //number of free tower param. p. bin//
					  p->tower_parametrization,                               //function//
					  p->tower_error_parametrization                          //function//
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
       trackcluster.TrackEt,  			          //truth//
       sqrt(pow(0.5,2)+pow(0.10*trackcluster.TrackEt,2)), //error//
       0,                                                 //params
       0,                                                 //number of free jet param. p. bin
							 p->dummy_parametrization,                             //function
							 p->jet_error_parametrization                          //function
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
       double * mess = new double[4];
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
			 p->GetTowerParRef( index ),                                //parameter//
	                 p->free_pars_per_bin,                                      //number of free cluster param. p. bin//
						     p->tower_parametrization,                                     //function//
						     p->tower_error_parametrization                                //function//
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

//--------------------------------------------------------------------------------------------
void TCaliber::Run()
{
  if (fit_method!=3){
    if (n_gammajet_events!=0)     Run_GammaJet();
    if (n_tracktower_events!=0)   Run_TrackTower();
    if (n_trackcluster_events!=0) Run_TrackCluster();

    if (fit_method==1) Run_Lvmini();
  } 
  //Dummy Configuration: Nothing to be done, start-values are written to file
}

void TCaliber::Run_Lvmini()
{
  int naux = 1000000, niter=1000, iflag, iret;
  //int mvec = 6;
  int mvec = 29;
  double aux[naux], fsum, fopt, fedm, dummy;

  int npar = p->GetNumberOfParameters();

  cout << "\nFitting " << npar << " parameters; \n" 
       << p->eta_granularity << " x " << p->phi_granularity << " tower bins with " 
       << p->free_pars_per_bin << " free parameters each, or " 
       << p->GetNumberOfTowerParameters() << " in total, and\n"
       << p->eta_granularity_jet << " x " << p->phi_granularity_jet << " JES bins with " 
       << p->free_pars_per_bin_jet << " free parameters each, or " 
       << p->GetNumberOfJetParameters() << " in total with LVMINI.\n" << endl;

  fitfunction    = this->global_fit_fast;

  float eps =float(1.E-3*data.size());
  float wlf1=1.E-4;
  float wlf2=0.9;
  lvmeps_(eps,wlf1,wlf2);

  //Set errors per default to 0 //@@ doesn't seem to work...
  int error_index=2;
  error_index = lvmind_(error_index);
  std::memcpy(aux+error_index,p->e,npar*sizeof(double));
  //for (int n=0; n<naux; ++n) aux[n]=0.0; 

  //outlier rejection before first iteration
  {
    std::vector<TData*>::iterator beg =  partition(data.begin(), data.end(), OutlierRejection(OutlierChi2CutPresel));
    for(std::vector<TData*>::iterator i = beg ; i != data.end() ; ++i) {
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
    do {
      fitfunction(npar, aux, fsum, p->k, iflag);  

      ////classic derivative calculation: (inadvisable !)
      //for (unsigned param=0; param<abs(npar); ++param) {
      ////first numeric derivative of fitfunction w.r.t. par[p]     
      //aux[param] = numeric_derivate(fitfunction, p->k, npar, param); 
      //
      ////first analytic derivative of fitfunction w.r.t. par[p]     
      ////aux[p] = analytic_derivate(p->k, npar, p); 
      //
      ////second numeric derivative of fitfunction w.r.t. par[p]     
      //	aux[param+npar] = numeric_derivate(fitfunction, p->k, npar, param+npar); 
      //}

      //fast derivative calculation:
      for (unsigned param=0; param<abs(npar); ++param) {
         aux[param]           = (TData::temp_derivative2[param]-TData::temp_derivative1[param])/(2.0*TData::epsilon);
         aux[param+abs(npar)] = (TData::temp_derivative2[param]+TData::temp_derivative1[param]-2*fsum)/(TData::epsilon*TData::epsilon);
      	}
      lvmfun_(p->k,fsum,iret,aux);
      
    }
    while (iret<0);
    //lvmprt_(2,aux,2); //Has any effect?

    //outlier rejection
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
  std::memcpy(p->e,aux+error_index,npar*sizeof(double));
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
  cout << "Creating control plots,"<<endl;
  plots->FitControlPlots();
  plots->GammaJetControlPlots();
  plots->GammaJetControlPlotsJetBin();
  plots->TrackTowerControlPlots();
  plots->TrackClusterControlPlots();
  
  //Clean-up
  cout << "Done, cleaning up."<<endl;
  delete p;
  delete plots;
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

  //init parameter and plot classes
  string param_class = config.read<string>("Parametrization Class","TParameters");

  //This is hard coded, no change of parametrization by config file possible
  if      (param_class=="TStepEfracParameters") p = new TStepEfracParameters( file );

//  if      (param_class=="TStepParameters") p = new TStepParameters( file );
//  else if (param_class=="TMyParameters")   p = new TMyParameters( file );
//  else                                     p = new TParameters( file );

  plots = new TControlPlots(file, &data, p);

  //initialize temp arrays for fast derivative calculation
  TData::total_n_pars     = p->GetNumberOfParameters();
  TData::temp_derivative1 = new double[p->GetNumberOfParameters()];
  TData::temp_derivative2 = new double[p->GetNumberOfParameters()];

  //--------------------------------------------------------------------------
  //read config file
  fit_method = config.read<int>("Fit method",1); 
  //last minute kinematic cuts
  Et_cut_on_jet   = config.read<double>("Et cut on jet",5.0); 
  Et_cut_on_gamma = config.read<double>("Et cut on gamma",20.0); 
  Et_cut_on_track = config.read<double>("Et cut on track",4.0); 
  Et_cut_on_tower = config.read<double>("Et cut on tower",1.0);
  Et_cut_on_cluster = config.read<double>("Et cut on cluster",1.0);
  //outlier rejection
  OutlierIterationSteps = config.read<int>("Outlier Iteration Steps",3);
  OutlierChi2Cut        = config.read<double>("Outlier Cut on Chi2",100.0);
  OutlierChi2CutPresel  = config.read<double>("Outlier Cut on Chi2 presel",400.0);
  //input/output
  n_gammajet_events     = config.read<int>("use Gamma-Jet events",-1);
  n_tracktower_events   = config.read<int>("use Track-Tower events",-1);
  n_trackcluster_events = config.read<int>("use Track-Cluster events",-1);
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
  string treename_tracktower    = config.read<string>( "Track-Tower tree", default_tree_name );
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
  for(std::vector<TData*>::iterator i = data.begin() ; i < data.end() ; ++i) {
    delete *i;
  }
  data.clear();
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
