//
// Original Author:  Christian Autermann
//         Created:  Wed Jul 18 13:54:50 CEST 2007
// $Id: caliber.cc,v 1.72 2008/12/27 16:35:13 stadie Exp $
//
//
// for profiling:
//  - to prevent gprof from missing the threads: 
//      wget http://sam.zoy.org/writings/programming/gprof-helper.c
//      gcc -shared -fPIC gprof-helper.c -o gprof-helper.so -lpthread -ldl 
//      LD_PRELOAD=./gprof-helper.so ./junk
//
#include "caliber.h"


#include <algorithm>

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
#include "EventProcessor.h"
#include "Jet.h"

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
	boost::mutex::scoped_lock lock(io_mutex);
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
      temp_derivative2(temp_derivative2), epsilon(epsilon), data_changed(false) {
    //std::cout << "threads par array:" << mypar << '\n';
  }
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



//--------------------------------------------------------------------------------------------
void TCaliber::Run()
{
  if (fit_method!=3){

    time_t start = time(0);
    
    //calls FlattenSpectra and BalanceSpectra if enabled in config
    EventProcessor ep(configfile,p);
    ep.process(data);

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

  double* aux = new double[naux], fsum = 0;

  double *temp_derivative1 = new double[npar];
  double *temp_derivative2 = new double[npar];

  cout << "\nFitting " << npar << " parameters; \n";
  p->Print();
  cout << " with LVMINI.\n" << "Using " << data.size() << " total events and ";
  cout << nthreads << " threads.\n";
  
  ComputeThread *t[nthreads];
  for (int ithreads=0; ithreads<nthreads; ++ithreads){
    t[ithreads] = new ComputeThread(npar, p->GetPars(),temp_derivative1,temp_derivative2,deriv_step);
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
    npar=std::abs(npar);

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
      for( int param = 0 ; param < npar ; ++param ) {
	aux[param]      = temp_derivative1[param]/(2.0*deriv_step);
	aux[param+npar] = temp_derivative2[param]/(deriv_step*deriv_step);
      }
      
      //print derivatives:
      if(print_parnderiv) {
	for( int param = 0 ; param < npar ; ++param ) {
	  std::cout << "par: " << param << ": p = " <<  p->GetPars()[param] << " dp/dx = " 
		    <<  aux[param] << " dp/dx^2 = " << aux[param+npar] << std::endl;
	}
      }
      lvmfun_(p->GetPars(),fsum,iret,aux);
      //p->SetParameters(aux + par_index); 
      lvmprt_(2,aux,2); //print out
    } while (iret<0);
   
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
  delete [] aux;
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

  //BFGS fit parameters
  deriv_step = config.read<double>("BFGS derivative step",1e-03);
  eps        = config.read<double>("BFGS eps",1e-02);
  wlf1       = config.read<double>("BFGS 1st wolfe parameter",1e-04);
  wlf2       = config.read<double>("BFGS 2nd wolfe parameter",0.9);
  print_parnderiv = config.read<bool>("BFGS print derivatives",false);

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
  Jet::printInversionStats();
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

