//  $Id: Kalibri.cc,v 1.17 2011/05/26 07:42:53 mschrode Exp $

#include "Kalibri.h"

#include <algorithm>
#include <iomanip>

//Boost
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
boost::mutex io_mutex;
// User
#include "ConfigFile.h"
#include "Parameters.h"
#include "ControlPlots.h"
#include "ControlPlotsResolution.h"
#include "CalibMath.h"
#include "external.h"
#include "ToyMC.h"
#include "PhotonJetReader.h"
#include "ThreadedDiJetReader.h"
#include "TriJetReader.h"
#include "ZJetReader.h"
#include "TopReader.h"
#include "ParameterLimitsReader.h"
#include "EventProcessor.h"
#include "EventWeightProcessor.h"
#include "EventBinning.h"
#include "DiJetEventWeighting.h"

#include <dlfcn.h>

using namespace std;

typedef std::vector<Event*>::iterator DataIter;
typedef std::vector<Event*>::const_iterator DataConstIter;



//!  \brief Outlier Rejection
// -----------------------------------------------------------------
struct OutlierRejection {
  OutlierRejection(double cut):_cut(cut){};
  bool operator()(Event *d){
    if(d->type()==typeTowerConstraint) return true;
    return (d->chi2()/d->weight())<_cut;
  }
  double _cut;
};



// -----------------------------------------------------------------
class ComputeThread {
private:
  int npar_;
  double chi2_;
  double * td1_;
  double * td2_;
  Parameters *parorig_, *mypar_;
  const double *epsilon_;
  std::vector<Event*> data_;
  bool data_changed_;
  struct calc_chi2_on
  {
  private:
    ComputeThread *parent_;
  public:
    calc_chi2_on(ComputeThread *parent) : parent_(parent) {}
    void operator()()
    {
      {
       	//boost::mutex::scoped_lock lock(io_mutex);
      	//std::cout << "start Thread for " << parent_ << "  " 
	//	  << parent_->data_changed_ << " " << parent_->mypar_->parameters() << std::endl;
      }   
      if(parent_->data_changed_) {
	for (DataIter it=parent_->data_.begin() ; it!= parent_->data_.end() ; ++it) {
	  (*it)->setParameters(parent_->mypar_);
	} 
	parent_->data_changed_ = false;
      }
      for (int param = 0 ; param < parent_->npar_ ; ++param) {
	parent_->td1_[param]= 0.0;
	parent_->td2_[param]= 0.0;
	parent_->mypar_->parameters()[param] = parent_->parorig_->parameters()[param];
      }
      parent_->chi2_ =0.0;   
      for (DataIter it=parent_->data_.begin() ; it != parent_->data_.end() ; ++it) { 
	//boost::mutex::scoped_lock lock(io_mutex);
	parent_->chi2_ += (*it)->chi2_fast(parent_->td1_, parent_->td2_, parent_->epsilon_);
      }
    }
  };
  boost::thread *thread_;
  friend class calc_chi2_on;
public:
  ComputeThread(int npar,Parameters *par, const double *epsilon) 
    : npar_(npar), td1_(new double[npar]), td2_(new double[npar]), parorig_(par),
      mypar_(par->clone()), epsilon_(epsilon), data_changed_(false) {
    //std::cout << "threads par array:" << mypar << '\n';
  }
  ~ComputeThread() {
    clearData();
    delete [] td1_;
    delete [] td2_;
    Parameters::removeClone(mypar_);
  }
  void addData(Event* d) { 
    //d->ChangeParAddress(parorig, mypar);
    data_changed_ = true;
    data_.push_back(d);
  }
  void clearData() {   
    for (DataIter it= data_.begin() ; it!= data_.end() ; ++it) {
      (*it)->setParameters(parorig_);
    }
    data_.clear();
  }
  void start() { thread_ = new boost::thread(calc_chi2_on(this)); }
  bool isDone() { thread_->join(); delete thread_; return true;}
  void syncParameters() {
    for (int param = 0 ; param < npar_ ; ++param) mypar_->parameters()[param] = parorig_->parameters()[param];
  }
  double chi2() const { return chi2_;}
  double tempDeriv1(int i) const { return td1_[i];}
  double tempDeriv2(int i) const { return td2_[i];}
};


//--------------------------------------------------------------------------------------------
void Kalibri::run()
{
  if (fitMethod_!=3){
    time_t start = time(0);
    
    std::vector<EventProcessor*> processors;
    processors.push_back(new EventWeightProcessor(configFile_,par_));
    processors.push_back(new EventBinning(configFile_,par_));
    processors.push_back(new DiJetEventWeighting(configFile_,par_));

    for(std::vector<EventProcessor*>::iterator i = processors.begin() ; i != processors.end() ; ++i) {
      (*i)->preprocess(data_,control_[0],control_[1]);
    } 

    if(! data_.size()) {
      std::cout << "Warning: No events to perform the fit!\n";
      return;
    } 
    if((fitMethod_==1) || (fitMethod_==0)) {
      int npar = par_->numberOfParameters();
      epsilon_ = new double[npar];
      temp_derivative1_ = new double[npar];
      temp_derivative2_ = new double[npar];
      
      cout << "\nFitting " << npar << " parameters; \n";
      par_->print();
      if(fitMethod_==1) 
	cout << " with LVMINI.\n"; 
      else 
	cout << " with lbfgs.\n"; 
      cout << "Using " << data_.size() << " total events and ";
      cout << nThreads_ << " threads.\n";
      
      // Fixed pars
      if( fixedJetPars_.size() > 0 ) cout << "Fixed jet parameters:\n";
      for(unsigned int i = 0; i < fixedJetPars_.size(); i++) {
	int idx = fixedJetPars_.at(i);
	cout << "  " << idx+1 << ": " << par_->parameters()[idx] << endl;
      }
      if( fixedGlobalJetPars_.size() > 0 ) cout << "Fixed global jet parameters:\n";
      for(unsigned int i = 0; i < fixedGlobalJetPars_.size(); i++) {
	int idx = fixedGlobalJetPars_.at(i);
	cout << "  " << idx+1 << ": " << par_->parameters()[idx] << endl;
      }
      threads_ = new ComputeThread*[nThreads_];
      for (int ithreads=0; ithreads<nThreads_; ++ithreads){
	threads_[ithreads] = new ComputeThread(npar, par_,epsilon_);
      }
      
      if(fitMethod_==1) {
	run_Lvmini();
	time_t end = time(0);
	cout << "Done, fitted " << par_->numberOfParameters() << " parameters in " << difftime(end,start) << " sec." << endl;
      } else if(fitMethod_==0) {
	run_lbfgs();
	time_t end = time(0);
	cout << "Done, fitted " << par_->numberOfParameters() << " parameters in " << difftime(end,start) << " sec." << endl;
      }
      for (int ithreads=0; ithreads<nThreads_; ++ithreads){
	delete threads_[ithreads];
      }
      delete [] threads_;
      delete [] epsilon_;
      delete [] temp_derivative1_;
      delete [] temp_derivative2_;
    }  
    else {
      if( par_->needsUpdate() ) par_->update();
    } 
    for(std::vector<EventProcessor*>::iterator i = processors.begin() ; i != processors.end() ; ++i) {
      (*i)->postprocess(data_,control_[0],control_[1]);
      delete *i;
    } 
  } 
}

// code for lbfgs

lbfgsfloatval_t Kalibri::lbfgs_evaluate(void *instance,
					const lbfgsfloatval_t *x,
					lbfgsfloatval_t *g,
					const int npar,
					const lbfgsfloatval_t step)
{
  Kalibri *k = static_cast<Kalibri*>(instance);

  lbfgsfloatval_t fsum = 0.0;
   
  for (int param=0; param< npar ; ++param) {
    k->temp_derivative1_[param]=0.0;
    k->temp_derivative2_[param]=0.0;
  } 
  //computed step sizes for derivative calculation
  if(k->printParNDeriv_) std::cout << "new par:\n";
  for(int param = 0 ; param < npar ; ++param) {
    k->par_->parameters()[param] = x[param]; 
    if(k->printParNDeriv_)  {
      std::cout << std::setw(5) << param;
      std::cout << std::setw(15) << k->par_->parameters()[param];
    }
    //k->epsilon_[param] = step ? step : k->derivStep_;
    k->epsilon_[param] =  k->derivStep_ * std::abs( k->par_->parameters()[param]);
    if( k->epsilon_[param] < 1e-06)  k->epsilon_[param] = 1e-06;
  }
  if(k->printParNDeriv_) std::cout << std::endl;
  //use zero step for fixed pars
  for( std::vector<int>::const_iterator iter = k->fixedJetPars_.begin();
       iter != k->fixedJetPars_.end() ; ++ iter) {
    k->epsilon_[*iter] = 0;
  }

  for (int ithreads=0; ithreads < k->nThreads_; ++ithreads) k->threads_[ithreads]->start();
  for (int ithreads=0; ithreads < k->nThreads_; ++ithreads){
    if(k->threads_[ithreads]->isDone()) {
      fsum += k->threads_[ithreads]->chi2();
      for (int param=0 ; param < npar ; ++param) {
	assert(k->threads_[ithreads]->tempDeriv1(param) == k->threads_[ithreads]->tempDeriv1(param));
	k->temp_derivative1_[param] += k->threads_[ithreads]->tempDeriv1(param);
	k->temp_derivative2_[param] += k->threads_[ithreads]->tempDeriv2(param);
      }
    }
  }
  for( std::vector<int>::const_iterator iter = k->fixedGlobalJetPars_.begin();
       iter != k->fixedGlobalJetPars_.end() ; ++ iter) {
    k->temp_derivative1_[*iter] = 0;
    k->temp_derivative2_[*iter] = 0;
  }
  fsum *= 0.5;//lvmini uses log likelihood not chi2
  //fast derivative calculation:
  for( int param = 0 ; param < npar ; ++param ) {
    //std::cout << "hier:" << step << " " << k->epsilon_[param] << " " 
    //	      << k->temp_derivative1_[param] << '\n';
    if(k->epsilon_[param] > 0) {
      g[param]      = 0.5 * k->temp_derivative1_[param]/(2.0*k->epsilon_[param]);
    } else {
      g[param] = 0;
    }
    if(g[param] != g[param]) {
      std::cout << "bad derivative: par = " << param << " td = " << k->temp_derivative1_[param] << " epsilon = " << k->epsilon_[param] << '\n';
    }
    assert(g[param] == g[param]);
  }
  //print derivatives:
  if(k->printParNDeriv_) {
    std::cout << std::setw(5) << "\npar";
    std::cout << std::setw(15) << "p";
    std::cout << std::setw(15) << "dp/dx";
    std::cout << std::setw(15) << "d^2p/dx^2\n";
    for( int param = 0 ; param < npar ; ++param ) {
      std::cout << std::setw(5) << param;
      std::cout << std::setw(15) << k->par_->parameters()[param];
      std::cout << std::setw(15) << g[param] 
		<< std::setw(15) << 0 << std::endl;
    }
    std::cout << "fsum:" << fsum << std::endl;
  }
  assert( fsum > 0 );
  return fsum;
}

int Kalibri::lbfgs_progress(void *instance,
			    const lbfgsfloatval_t *x,
			    const lbfgsfloatval_t *g,
			    const lbfgsfloatval_t fx,
			    const lbfgsfloatval_t xnorm,
			    const lbfgsfloatval_t gnorm,
			    const lbfgsfloatval_t step,
			    int n, int k, int ls)
{
  printf("Iteration %d:\n", k);
  printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
  printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
  printf("\n");
  return 0;
}

void Kalibri::run_lbfgs()
{
  int npar = par_->numberOfParameters();
  lbfgsfloatval_t fx;
  lbfgsfloatval_t *x = lbfgs_malloc(npar);
  lbfgs_parameter_t param;
  
  if (x == NULL) {
    printf("ERROR: Failed to allocate a memory block for variables.\n");
    return;
  }

  /*
    Start the L-BFGS optimization; this will invoke the callback functions
    evaluate() and progress() when necessary.
  */  
  for( unsigned int loop = 0; loop < residualScalingScheme_.size() ; ++loop ) {
    cout<<"Updating Di-Jet Errors"<<endl;
    for(DataIter it = data_.begin()  ; it < data_.end() ; ++it) {
      (*it)->updateError();
    }
    if( par_->needsUpdate() ) par_->update();
    // Setting function to scale residuals in chi2 calculation
    cout << loop+1 << flush;
    if(  loop+1 == 1  ) cout << "st" << flush;
    else if(  loop+1 == 2  ) cout << "nd" << flush;
    else if(  loop+1 == 3  ) cout << "rd" << flush;
    else cout << "th" << flush;
    cout << " of " << residualScalingScheme_.size() <<" iteration(s)" << flush;
    if( mode_ == 0 ) {
      if(  residualScalingScheme_.at(loop) == 0  ) {
	Event::scaleResidual = &Event::scaleNone;	
	cout << ": no scaling of residuals." << endl;
	
	cout << "Rejecting outliers " << flush;
	DataIter beg = partition(data_.begin(), data_.end(), OutlierRejection(outlierChi2Cut_));
	for(DataIter i = beg ; i != data_.end() ; ++i) {
	  delete *i;
	}
	data_.erase(beg,data_.end());
	cout << "and using " << data_.size() << " events." << endl;
      }
      else if(  residualScalingScheme_.at(loop) == 1  ) {
	Event::scaleResidual = &Event::scaleCauchy;	
	cout << ": scaling of residuals with Cauchy-Function." << endl;
      }
      else if(  residualScalingScheme_.at(loop) == 2  ) {
	Event::scaleResidual = &Event::scaleHuber;	
	cout << ": scaling of residuals with Huber-Function." << endl;
      }
      else if(  residualScalingScheme_.at(loop) == 3  ) {
	Event::scaleResidual = &Event::scaleTukey;	
	cout << ": scaling of residuals a la Tukey." << endl;
      }
      else {
	cerr << "ERROR: " << residualScalingScheme_.at(loop) << " is not a valid scheme for resdiual scaling! Breaking iteration!" << endl;
	break;
      }
    } else {
      std::cout << std::endl;
    }
    /*param.linesearch = LBFGS_LINESEARCH_BACKTRACKING;*/
    /* Initialize the parameters for the L-BFGS optimization. */ 
    for(int i = 0 ; i < npar ; ++i) {
      x[i] = par_->parameters()[i];
    } 
    lbfgs_parameter_init(&param);
    int n = 0;
    for(DataIter it = data_.begin()  ; it < data_.end() ; ++it) {
      threads_[n]->addData(*it);
      n++;
      if(n == nThreads_) n = 0;
    }
    
    int ret = lbfgs(npar, x, &fx, lbfgs_evaluate, lbfgs_progress, this, &param);
    /* Report the result. */
    printf("L-BFGS optimization terminated with status code = %d\n", ret);
    printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
    
    for (int ithreads=0; ithreads<nThreads_; ++ithreads){
      threads_[ithreads]->clearData();
    }  
  }
  lbfgs_free(x); 
}



// -----------------------------------------------------------------
void Kalibri::run_Lvmini()
{ 
  int naux = 3000000, iret=0;
  
  int npar = par_->numberOfParameters();
  int mvec = mvec_;
  if( calcCov_ ) mvec = -mvec_;

  naux = lvmdim_(npar,mvec);
  cout<<"array of size "<<naux<<" needed."<<endl;

  double* aux = new double[naux], fsum = 0;
 

 
  //lvmeps_(data_.size()*eps_,wlf1_,wlf2_);
  lvmeps_(eps_,wlf1_,wlf2_);

  //Set errors per default to 0 //@@ doesn't seem to work...
  int error_index=2;
  error_index = lvmind_(error_index);
  par_->fillErrors(aux+error_index);

  for( unsigned int loop = 0; loop < residualScalingScheme_.size() ; ++loop ) {
    cout<<"Updating Di-Jet Errors"<<endl;
    for(DataIter it = data_.begin()  ; it < data_.end() ; ++it) {
      (*it)->updateError();
    }
    if( par_->needsUpdate() ) par_->update();
    // Setting function to scale residuals in chi2 calculation
    cout << loop+1 << flush;
    if(  loop+1 == 1  ) cout << "st" << flush;
    else if(  loop+1 == 2  ) cout << "nd" << flush;
    else if(  loop+1 == 3  ) cout << "rd" << flush;
    else cout << "th" << flush;
    cout << " of " << residualScalingScheme_.size() <<" iteration(s)" << flush;
    if( mode_ == 0 ) {
      if(  residualScalingScheme_.at(loop) == 0  ) {
	Event::scaleResidual = &Event::scaleNone;	
	cout << ": no scaling of residuals." << endl;
	
	cout << "Rejecting outliers " << flush;
	DataIter beg = partition(data_.begin(), data_.end(), OutlierRejection(outlierChi2Cut_));
	for(DataIter i = beg ; i != data_.end() ; ++i) {
	  delete *i;
	}
	data_.erase(beg,data_.end());
	cout << "and using " << data_.size() << " events." << endl;
      }
      else if(  residualScalingScheme_.at(loop) == 1  ) {
	Event::scaleResidual = &Event::scaleCauchy;	
	cout << ": scaling of residuals with Cauchy-Function." << endl;
      }
      else if(  residualScalingScheme_.at(loop) == 2  ) {
	Event::scaleResidual = &Event::scaleHuber;	
	cout << ": scaling of residuals with Huber-Function." << endl;
      }
      else if(  residualScalingScheme_.at(loop) == 3  ) {
	Event::scaleResidual = &Event::scaleTukey;	
	cout << ": scaling of residuals a la Tukey." << endl;
      }
      else {
	cerr << "ERROR: " << residualScalingScheme_.at(loop) << " is not a valid scheme for resdiual scaling! Breaking iteration!" << endl;
	break;
      }
    } else {
      std::cout << std::endl;
    }
    if(lvmdim_(npar,mvec) > naux)
      cout<<"Aux field too small. "<<lvmdim_(npar,mvec)<<" enntires needed."<<endl;
  
    //initialization
    int nparm = -npar; //Show output
    lvmini_(nparm, mvec, nIter_, aux);
    
    int n = 0;
    
    for(DataIter it = data_.begin()  ; it < data_.end() ; ++it) {
      threads_[n]->addData(*it);
      n++;
      if(n == nThreads_) n = 0;
    }
    do {
      //set storage for temporary derivative storage to zero
      for (int param=0; param< npar ; ++param) {
	temp_derivative1_[param]=0.0;
	temp_derivative2_[param]=0.0;
      } 
      //computed step sizes for derivative calculation
      if(printParNDeriv_) std::cout << "new par:\n";
      for(int param = 0 ; param < npar ; ++param) {
	if(printParNDeriv_)  {
	  std::cout << std::setw(5) << param;
	  std::cout << std::setw(15) << par_->parameters()[param];
	}
	epsilon_[param] = derivStep_ * std::abs(par_->parameters()[param]);
	if(epsilon_[param] < 1e-06) epsilon_[param] = 1e-06;
      }
      if(printParNDeriv_) std::cout << std::endl;
      //use zero step for fixed pars
      for( std::vector<int>::const_iterator iter = fixedJetPars_.begin();
	   iter != fixedJetPars_.end() ; ++ iter) {
	epsilon_[*iter] = 0;
      }
      fsum = 0;
      for (int ithreads=0; ithreads<nThreads_; ++ithreads) threads_[ithreads]->start();
      for (int ithreads=0; ithreads<nThreads_; ++ithreads){
	if(threads_[ithreads]->isDone()) {
	  fsum += threads_[ithreads]->chi2();
	  for (int param=0 ; param < npar ; ++param) {
	    assert(threads_[ithreads]->tempDeriv1(param) == threads_[ithreads]->tempDeriv1(param));
	    temp_derivative1_[param] += threads_[ithreads]->tempDeriv1(param);
	    temp_derivative2_[param] += threads_[ithreads]->tempDeriv2(param);
	  }
	}
      }
      for( std::vector<int>::const_iterator iter = fixedGlobalJetPars_.begin();
	   iter != fixedGlobalJetPars_.end() ; ++ iter) {
	temp_derivative1_[*iter] = 0;
	temp_derivative2_[*iter] = 0;
      }
      fsum *= 0.5;//lvmini uses log likelihood not chi2
      //fast derivative calculation:
      for( int param = 0 ; param < npar ; ++param ) {
	if(epsilon_[param] > 0) {
	aux[param]      = 0.5 * temp_derivative1_[param]/(2.0*epsilon_[param]);
	aux[param+npar] = 0.5 * temp_derivative2_[param]/(epsilon_[param]*epsilon_[param]);
	} else {
	  aux[param] = 0;
	  aux[param+npar] = 0;
	}
	if(aux[param] != aux[param]) {
	  std::cout << "bad derivative: par = " << param << " td = " << temp_derivative1_[param] << " epsilon_ = " << epsilon_[param] << '\n';
	}
 	assert(aux[param] == aux[param]);
 	assert(aux[param+npar] == aux[param+npar]);
      }
      //lvmout_(npar,mvec_,aux);
      //print derivatives:
      if(printParNDeriv_) {
	std::cout << std::setw(5) << "\npar";
	std::cout << std::setw(15) << "p";
	std::cout << std::setw(15) << "dp/dx";
	std::cout << std::setw(15) << "d^2p/dx^2\n";
	for( int param = 0 ; param < npar ; ++param ) {
	  std::cout << std::setw(5) << param;
	  std::cout << std::setw(15) << par_->parameters()[param];
	  std::cout << std::setw(15) << aux[param];
	  std::cout << std::setw(15) << aux[param+npar] << std::endl;
	}
	std::cout << "fsum:" << fsum << std::endl;
      }
      assert( fsum > 0 );
      lvmfun_(par_->parameters(),fsum,iret,aux);
      //lvmout_(npar,mvec_,aux);
    } while (iret<0); 
    for (int ithreads=0; ithreads<nThreads_; ++ithreads){
      threads_[ithreads]->clearData();
    }  
    int par_index = 1;
    par_index = lvmind_(par_index);
    par_->setParameters(aux + par_index);
  }
  //Copy Parameter errors from aux array to the TParameter::e array
  if( !calcCov_ ) {
    error_index=2;
    error_index = lvmind_(error_index);
    par_->setErrors(aux+error_index);
  } else {
    // Retrieve parameter errors
    error_index = 3;
    error_index = lvmind_(error_index);
    par_->setErrors(aux+error_index);
    // Retrieve global parameter correlation coefficients
    error_index = 4;
    error_index = lvmind_(error_index);
    bool nanOccured = false;
    for(int i = 0; i < par_->numberOfParameters(); i++) {
      if( aux[error_index+i] != aux[error_index+i] ) { // Check for NAN
	if( !nanOccured ) {
	  nanOccured = true;
	  std::cout << "The following global correlation coefficients are NAN and set to 0:\n";
	}
	std::cout << i << std::endl;
	aux[error_index+i] = 0.;
      }
    }
    par_->setGlobalCorrCoeff(aux+error_index);
    // Retrieve parameter covariances
    error_index = 5;
    error_index = lvmind_(error_index);
    // Set cov = 0 for fixed parameters
    nanOccured = false;
    for(int i = 0; i < par_->numberOfCovCoeffs(); i++) {
      if( aux[error_index+i] != aux[error_index+i] ) { // Check for NAN
	if( !nanOccured ) {
	  nanOccured = true;
	  std::cout << "The following covariance elements are NAN and set to 0:\n";
	}
	std::cout << i << std::endl;

	aux[error_index+i] = 0.;
      }
    }
    par_->setCovCoeff(aux+error_index);
  }
  delete [] aux;  
}



//--------------------------------------------------------------------------------------------
void Kalibri::done()
{
  std::cout << "Kalibri::done()\n";
  if( par_->needsUpdate() ) par_->update();
  
  ConfigFile config( configFile_.c_str() );

  bool txt=false;
  bool tex=false;
  std::string fileName(getOutputFile());
  // write calibration to txt output file if ending is txt
  if( fileName.find(".txt")!=std::string::npos ){
    if( fileName.substr(fileName.find(".txt")).compare(".txt")==0 ){
      par_->writeCalibrationTxt( fileName.c_str() );
      txt=true; // file has a real .txt ending
    }
  }
  // write calibration to tex output file if ending is tex
  if( fileName.find(".tex")!=std::string::npos ){
    if( fileName.substr(fileName.find(".tex")).compare(".tex")==0 ){
      par_->writeCalibrationTex( fileName.c_str(), config );
      tex=true; // file has a real .txt ending
    }
  }

  // write calibration to txt and tex output file if w/o ending
  if( !txt && !tex ){
    par_->writeCalibrationTxt( (fileName+".txt").c_str() );
    par_->writeCalibrationTex( (fileName+".tex").c_str(), config );
  }


  // Make control plots
  if( config.read<bool>("create plots",0) ) {
    if( mode_ == 0 ) {  // Control plots for calibration
      std::vector<std::vector<Event*>* > samples;
      samples.push_back(&data_);
      samples.push_back(&control_[0]);
      samples.push_back(&control_[1]);      
      ControlPlots * plots = new ControlPlots(&config,samples);
      plots->makePlots();
      delete plots;
    } else if( mode_ == 1 ) {  // Control plots for jetsmearing
      ControlPlotsResolution * plotsjs = new ControlPlotsResolution(configFile_,&data_,par_);
      plotsjs->makePlots();
      delete plotsjs;
    }
  }
  
  // Clean-up
  cout << endl << "Cleaning up... " << flush;
  for(DataIter i = data_.begin() ; i != data_.end() ; ++i) {
    delete *i;
  }
  data_.clear(); 
  for(DataIter i = control_[0].begin() ; i != control_[0].end() ; ++i) {
    delete *i;
  }
  control_[0].clear(); 
  for(DataIter i = control_[1].begin() ; i != control_[1].end() ; ++i) {
    delete *i;
  }
  control_[1].clear();
  cout << "Done" << endl;
}



//--------------------------------------------------------------------------------------------
void Kalibri::init()
{
  ConfigFile config(configFile_.c_str() );

  par_ = Parameters::createParameters(configFile_);


  //--------------------------------------------------------------------------
  //read config file
  mode_ = config.read<int>("Mode",0);
  fitMethod_ = config.read<int>("Fit method",1);
  nThreads_ = config.read<int>("Number of Threads",1);
  
  // Residual scaling
  const char* resScheme = ( config.read<string>("Residual Scaling Scheme","221").c_str() );
  while(  *resScheme != 0  )
    {
      int scheme = static_cast<int>(*resScheme - '0');
      if(  scheme < 0  ||  scheme > 3  )
	{
	  cerr << "ERROR: " << scheme << " is not a valid scheme for resdiual scaling! Using default scheme 221." << endl << endl;
	  residualScalingScheme_.clear();
	  residualScalingScheme_.push_back(2);
	  residualScalingScheme_.push_back(2);
	  residualScalingScheme_.push_back(1);
	  break;
	}

      residualScalingScheme_.push_back( static_cast<int>(*resScheme - '0') );
      resScheme++;
    }


  outlierChi2Cut_        = config.read<double>("Outlier Cut on Chi2",100.0);

  //BFGS fit parameters
  derivStep_ = config.read<double>("BFGS derivative step",1e-03);
  mvec_       = config.read<int>("BFGS mvec",6);
  nIter_      = config.read<int>("BFGS niter",100);
  eps_        = config.read<double>("BFGS eps",1e-02);
  wlf1_       = config.read<double>("BFGS 1st wolfe parameter",1e-04);
  wlf2_       = config.read<double>("BFGS 2nd wolfe parameter",0.9);
  calcCov_    = config.read<bool>("BFGS calculate covariance",false);
  printParNDeriv_ = config.read<bool>("BFGS print derivatives",false);

  //fixed jet parameters
  std::vector<int> fixJetPars = bag_of<int>(config.read<string>("fixed jet parameters",""));
  if(fixJetPars.size() % 3 == 0) {
    // Fix the specified parameters
    for(unsigned int i = 0 ; i < fixJetPars.size() ; i += 3) {
      int etaid = fixJetPars[i];
      int phiid = fixJetPars[i+1];
      int parid = fixJetPars[i+2];
      if(parid >= par_->numberOfJetParametersPerBin()) continue;
      int jetbin = par_->jetBin(par_->jetEtaBin(etaid),par_->jetPhiBin(phiid));
      if(jetbin < 0) {
	std::cerr<<"WARNING: fixed jet parameter bin index = " << jetbin << endl; 
	exit(-2);  
      }
      //std::cout << "jetbin:" << jetbin << '\n';
      fixedJetPars_.push_back(jetbin * par_->numberOfJetParametersPerBin() + par_->numberOfTowerParameters() + parid);
    }
  }
  else if( fixJetPars.size() % 2 == 0 ) {
    // Fix all parameters in the specified jet bins
    for(unsigned int i = 0 ; i < fixJetPars.size() ; i += 2) {
      int etaid = fixJetPars[i];
      int phiid = fixJetPars[i+1];
      int jetbin = par_->jetBin(par_->jetEtaBin(etaid),par_->jetPhiBin(phiid));
      if(jetbin < 0) {
	std::cerr<<"WARNING: fixed jet parameter bin index = " << jetbin << endl; 
	exit(-2);  
      }

      for(int parIdx = 0; parIdx < par_->numberOfJetParametersPerBin(); parIdx++) {
	fixedJetPars_.push_back( par_->numberOfTowerParameters() + 
				 jetbin * par_->numberOfJetParametersPerBin() +
				 parIdx );
      }
    }
  } else {
//     cerr << "ERROR: Syntax error for fixed jet parameters. Syntax is:\n";
//     cerr << "       'fixed jet parameter = { <eta_id> <phi_id> <par_id> }' or\n"; 
//     cerr << "       'fixed jet parameter = { <eta_id> <phi_id> }'\n"; 
    // Fix specified parameters in all jet bins
    for(int jetBin = 0; jetBin < par_->etaGranularityJet(); jetBin++) {
      for(unsigned int i = 0 ; i < fixJetPars.size() ; i++) {
	int parIdx = par_->numberOfTowerParameters()
	  + jetBin * par_->numberOfJetParametersPerBin()
	  + fixJetPars[i];
	fixedJetPars_.push_back(parIdx);
      }
    }
  }
  for( std::vector<int>::const_iterator iter = fixedJetPars_.begin();
       iter != fixedJetPars_.end() ; ++ iter) {
    par_->fixPar(*iter);
  }

  //fixed global parameters
  std::vector<int> fixGlobalJetPars = bag_of<int>(config.read<string>("fixed global jet parameters",""));
  for(size_t globalJetBin = 0; globalJetBin < fixGlobalJetPars.size(); globalJetBin++) {
      if( globalJetBin < 0 ) {
	std::cerr << "ERROR: fixed global jet parameter bin index = " << globalJetBin << std::endl; 
	exit(-2);  
      } else if( static_cast<int>(globalJetBin) > par_->numberOfGlobalJetParameters() ) {
	std::cerr << "ERROR: fixed global jet parameter bin index = " << globalJetBin;
	std::cerr << " which is larger than the max number ";
	std::cerr << par_->numberOfGlobalJetParameters() << " of global parameters." << std::endl;
	exit(-2);  
      } else {
	fixedGlobalJetPars_.push_back( par_->numberOfTowerParameters() +
				       par_->numberOfJetParameters()   +
				       par_->numberOfTrackParameters() +
				       globalJetBin );
      }
  }
  for( std::vector<int>::const_iterator iter = fixedGlobalJetPars_.begin();
       iter != fixedGlobalJetPars_.end() ; ++ iter) {
    par_->fixPar(*iter);
  }
  
  //load plugin
  std::string jcs = config.read<string>("jet correction source","");
  if(jcs != "") {
    std::string libname = "lib/lib"+jcs+".so";
    //std::cout << "loading lib " << libname << '\n';
    void *hndl = dlopen(libname.c_str(), RTLD_NOW);
    if(hndl == NULL){
      std::cerr << "failed to load plugin: " << dlerror() << std::endl;
      exit(-1);
    }
  }
  outputFile_ = config.read<string>( "Output file", "calibration_k.cfi" );

  int niothreads = config.read<int>("Number of IO Threads",1); 

  //int djdc = config.read<int>("Di-Jet data class", 0);
  //fill data vector  
  std::vector<EventReader*> readers;
  readers.push_back(new PhotonJetReader(configFile_,par_));
  if((niothreads < 1)) {
    readers.push_back(new DiJetReader(configFile_,par_));
  } else {
    std::cout << "reading DiJet data in " << niothreads << " threads.\n";
    readers.push_back(new ThreadedDiJetReader(configFile_,par_,niothreads));
  }
  readers.push_back(new TriJetReader(configFile_,par_));
  readers.push_back(new ZJetReader(configFile_,par_));
  readers.push_back(new TopReader(configFile_,par_));
  readers.push_back(new ParameterLimitsReader(configFile_,par_));
  for(std::vector<EventReader*>::iterator i = readers.begin() ; 
      i != readers.end() ; ++i) {
    (*i)->readEvents(data_);
    (*i)->readControlEvents(control_[0],1);
    (*i)->readControlEvents(control_[1],2);
  }
  for(std::vector<EventReader*>::iterator i = readers.begin() ; 
     i != readers.end() ; ++i) {
    delete *i;
  }
  EventReader::addConstraints(data_);
}
//--^-Kalibri class-^------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------

