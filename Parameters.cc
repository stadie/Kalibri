#include <fstream>
#include <pwd.h>
#include <unistd.h>
#include <time.h>

#include "Parameters.h"
#include "ConfigFile.h"


using namespace std;

TParameters* TParameters::instance = 0;

Parametrization* TParameters::CreateParametrization(const std::string& name) {
  if(name == "StepParametrization") {
    return new StepParametrization();
  } else if(name == "MyParametrization") {
    return new MyParametrization();
  } else if(name == "StepEfracParametrization") {
    return new StepEfracParametrization();
  }  else if(name == "JetMETParametrization") {
    return new JetMETParametrization();
  }
  return 0;
}

TParameters* TParameters::CreateParameters(const std::string& configfile) 
{
  static Cleaner cleanup;
  assert(instance == 0);
 
  ConfigFile config( configfile.c_str() );
  
  string parclass = config.read<string>("Parametrization Class","");
  //create Parameters
  if(parclass == "TStepParameters") {
    parclass = "StepParametrization";
  } else if(parclass == "TMyParameters") {
    parclass = "MyParametrization";
  } else if(parclass == "TStepEfracParameters") {
    parclass = "StepEfracParametrization";
  }  else if(parclass == "TJetMETParameters") {
    parclass = "JetMETParametrization";
  }
  Parametrization *param = CreateParametrization(parclass);
  if(! param) {
    cerr << "TParameters::CreateParameters: could not instantiate class " << parclass << '\n';
    exit(1);
  }
  instance = new TParameters(param);
  instance->Init(config);
  return instance;
}

void TParameters::Init(const ConfigFile& config)
{
  eta_ntwr_used   = config.read<unsigned>("maximum eta twr used",82); 
  eta_granularity = config.read<unsigned>("granularity in eta",1); 
  phi_granularity = config.read<unsigned>("granularity in phi",1); 
  eta_symmetry    = config.read<bool>("symmetry in eta",false);
  eta_granularity_jet = config.read<unsigned>("jet granularity in eta",1); 
  phi_granularity_jet = config.read<unsigned>("jet granularity in phi",1); 
  input_calibration   = config.read<string>("input calibration","");

  if (eta_ntwr_used%2 !=0){
    cerr << "WARNING: Use even number of eta towers! Forced exit."<< endl;    
    exit(1);
  }

  if (phi_ntwr%phi_granularity!=0) {
    cerr << "WARNING: Check phi granularity! Forced exit."<< endl;
    exit(1);
  }
      
  if (eta_symmetry && (eta_granularity!=1 && eta_granularity!=3 &&eta_granularity!=5&&eta_granularity!=11&&
      eta_granularity!=21 && eta_granularity!=41 )){
    cerr << "WARNING: Check eta granularity! Should be 1, 3, 5, 11, 21, or 41: Forced exit."<< endl;
    exit(1);
  }

  if (!eta_symmetry && (eta_granularity!=2 && eta_granularity!=6 &&eta_granularity!=10&&eta_granularity!=22&&
      eta_granularity!=42 && eta_granularity!=82 )){
    cerr << "WARNING: Check eta granularity! Should be 2, 6, 10, 22, 42 or 82: Forced exit."<< endl;
    exit(1);
  }

  start_values = bag_of<double>(config.read<string>("start values","3.8 -0.7  0.04")); 
  if ( start_values.size()< p->nTowerPars()){
    cerr<< "ERROR: Number of start values and free parameters does not match!"<<endl
        << "       There must be at least " << p->nTowerPars() << " parameters!" << endl;
    exit(2);    
  }
  jet_start_values = bag_of<double>(config.read<string>("jet start values","3.8 -0.7  0.04")); 
  if ( jet_start_values.size()< p->nJetPars()){
    cerr<< "ERROR: Number of jet start values and free jet parameters does not match!"<<endl
        << "       There must be at least " << p->nJetPars() << " parameters!" << endl;
    exit(3);
  }
  
  k = new double[GetNumberOfParameters()];
  e = new double[GetNumberOfParameters()];

  for (unsigned int bin=0; bin<eta_granularity*phi_granularity; ++bin){
    for (unsigned int tp=0; tp < p->nTowerPars(); ++tp){
      //step[ bin*free_pars_per_bin + tp ]   = step_sizes[ tp ];
      k[ bin*p->nTowerPars() + tp ] = start_values[ tp ];
      e[ bin*p->nTowerPars() + tp ] = 0.0;
    }
  }
  for (unsigned int bin=0; bin<eta_granularity_jet*phi_granularity_jet; ++bin){
    for (unsigned int jp=0; jp < p->nJetPars(); ++jp){
      int i = GetNumberOfTowerParameters() + bin*p->nTowerPars() + jp;   
      k[i] = jet_start_values[jp];
      e[i] = 0.0;
    }
  }
  Read_Calibration(input_calibration);
}

std::string TParameters::trim(std::string const& source, char const* delims) {
  std::string result(source);
  std::string::size_type index = result.find_last_not_of(delims);
  if(index != std::string::npos)
    result.erase(++index);

  index = result.find_first_not_of(delims);
  if(index != std::string::npos)
    result.erase(0, index);
  else
    result.erase();
  
  //replace all "," by " "  :
  std::string::size_type  pos = result.find(",");
  while(pos != string::npos) {
    result.replace(pos,1," ");
    pos = result.find(",",pos);
  }
    
  return result;
}

void TParameters::Read_Calibration(std::string const& configFile) {
  std::ifstream file(configFile.c_str());

  std::string line, name;
  char * dummy = new char[28];
  std::vector<int> eta, phi;
  std::vector<double> param[p->nTowerPars()], error[p->nTowerPars()];
  std::vector<int> eta_jet, phi_jet;
  std::vector<double> param_jet[p->nJetPars()], error_jet[p->nJetPars()];
  int posEqual;
  while (std::getline(file,line)) {
    if (! line.length()) continue;
    if( line.find("#") != string::npos) continue;
    
    //Read Tower Calibration: ---------------------------------------------------
    //if ( line.find("module ccctm = CalibratedCaloTowerMaker") != string::npos ) {
    if ( line.find("block TowerCalibConstants = {") != string::npos ) {
      while( std::getline(file,line) ) {
	if( line.find("block") != string::npos) break;
	posEqual=line.find('=');
	name  = line.substr(0,posEqual);
	if( name.find("TowMapEta") != string::npos) 
  	  eta = bag_of<int>(trim(line.substr(posEqual+1)));
	if( name.find("TowMapPhi") != string::npos) 
  	  phi = bag_of<int>(trim(line.substr(posEqual+1)));
	for (unsigned i=0; i < p->nTowerPars() ; ++i) {
	  sprintf(dummy,"TowerParam%d ",i);
	  if( name.find(dummy) != string::npos) 
  	    param[i] = bag_of<double>(trim(line.substr(posEqual+1)));
	  sprintf(dummy,"TowerError%d ",i);
	  if( name.find(dummy) != string::npos) 
  	    error[i] = bag_of<double>(trim(line.substr(posEqual+1)));
	}
      }
    }
    //Read Jet Calibration --------------------------------------------------------
    if ( line.find("block JetCalibConstants = {") != string::npos ) {
      while (std::getline(file,line)) {
	if( line.find("block") != string::npos) break;
	posEqual=line.find('=');
	name  = line.substr(0,posEqual);
	std::cout << name << ".\n";
	if( name.find("JetMapEta") != string::npos) 
	  eta_jet = bag_of<int>(trim(line.substr(posEqual+1)));
	if( name.find("JetMapPhi") != string::npos) 
  	  phi_jet = bag_of<int>(trim(line.substr(posEqual+1)));
	for (unsigned i=0; i<p->nJetPars(); ++i) {
	  sprintf(dummy,"JetParam%d ",i);
	  if( name.find(dummy) != string::npos) 
  	    param_jet[i] = bag_of<double>(trim(line.substr(posEqual+1)));
	  sprintf(dummy,"JetError%d ",i);
	  if( name.find(dummy) != string::npos) 
  	    error_jet[i] = bag_of<double>(trim(line.substr(posEqual+1)));
	}
      }
    }
  }
  //check if the read calibration is ok:
  bool ok=eta.size()==phi.size();
  for (unsigned i=0; i < p->nTowerPars(); ++i){
    ok &= eta.size()==param[i].size();
    ok &= eta.size()==error[i].size();
  }
  //fill tower parameters and errors:  
  if (ok) {
    for (unsigned i=0; i<eta.size(); ++i){
      int index = GetBin(GetEtaBin(eta[i]),GetPhiBin(phi[i]));
      if (index<0) continue;
      for (unsigned n=0; n< p->nTowerPars(); ++n) {
        k[index*p->nTowerPars()+n] = param[n][i];
        e[index*p->nTowerPars()+n] = error[n][i];
      }
    }
  }

  //check if the read calibration is ok:
  ok=eta_jet.size()==phi_jet.size();
  for (unsigned i=0; i<p->nJetPars(); ++i){
    ok &= eta_jet.size()==param_jet[i].size();
    ok &= eta_jet.size()==error_jet[i].size();
  } 
  //fill Jet parameters and errors:  
  if (ok) {
    for (unsigned i=0; i<eta_jet.size(); ++i){
      int index = GetJetBin(GetJetEtaBin(eta_jet[i]),GetJetPhiBin(phi_jet[i]));
      if (index<0) continue;
      for (unsigned n=0; n<p->nJetPars(); ++n) {
        k[GetNumberOfTowerParameters() + index*p->nJetPars()+n] = param_jet[n][i];
        e[GetNumberOfTowerParameters() + index*p->nJetPars()+n] = error_jet[n][i];
      }
    }
  }

  delete[] dummy;
}

int TParameters::GetEtaBin(int const eta_id) const
{
//This function knows the number of wanted eta-bins and returns 
//in which eta-bin the tower with eta-ID "eta_id" is located.
  //Case 1 bin:
//cout << "eta="<<eta_id<<", eta_granularity:"<< eta_granularity<< ", eta_ntwr_used:"<< eta_ntwr_used<<endl;
  if (eta_granularity<=1) return 0;
  if (eta_granularity==2) return (eta_id < 0) ? 0 : 1;

  //check if tower is within wanted etarange:
  if ( eta_symmetry && abs(eta_id)*2>(int)eta_ntwr_used)   return -2; 
  
  //calculate an index:
  unsigned index=(unsigned)(41+eta_id);
  if (eta_id>0) --index;
  if (index<0 || index>81) return -3;
   
  //Lookup tables for some binning scenarios:                      eta 1   ->|<-  eta 2                                                |<- eta 3                         
  unsigned ta_42[82]={ 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,20,20,21,21,22,23,23,24,24,25,25,26,26,27,27,28,28,29,29,30,30,31,31,32,32,33,33,34,34,35,35,36,36,37,37,38,38,39,39,40,40,41,41};
  unsigned ta_22[82]={ 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9,10,10,10,10,11,11,11,11,12,12,12,12,13,13,13,13,14,14,14,14,15,15,15,15,16,16,16,16,17,17,17,17,18,18,18,18,19,19,19,19,20,20,20,21,21};
  unsigned ta_10[82]={ 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9};
  unsigned ta_6[82]= { 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5};

  unsigned ts_21[41]={ 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,20,20};
  unsigned ts_11[41]={ 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9,10,10};
  unsigned ts_5[41]= { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4};
  unsigned ts_3[41]= { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2};
  if (!eta_symmetry){
    if (eta_granularity==82) return index;
    else if (eta_granularity==42) return ta_42[index];
    else if (eta_granularity==22) return ta_22[index];
    else if (eta_granularity==10) return ta_10[index];
    else if (eta_granularity==6) return ta_6[index];
  } else {
    if (eta_granularity==41) return abs(eta_id)-1;
    else if (eta_granularity==21) return ts_21[abs(eta_id)-1];
    else if (eta_granularity==11) return ts_11[abs(eta_id)-1];
    else if (eta_granularity== 5) return ts_5[ abs(eta_id)-1];
    else if (eta_granularity== 3) return ts_3[ abs(eta_id)-1];
  }
  
  //Default value, should never be returned!
  return -4;
}

int TParameters::GetPhiBin(int const phi_id) const
//This function knows the number of wanted phi-bins and returns 
//in which phi-bin the tower with eta-ID "phi_id" is located.
{
  return (phi_id-1)*phi_granularity/phi_ntwr;
}

void TParameters::Print() const
{
  std::cout  << p->name() << " resulting in:\n "
    << eta_granularity << " x " << phi_granularity << " tower bins with " 
    << GetNumberOfTowerParametersPerBin() << " free parameters each, or " 
    << GetNumberOfTowerParameters() << " in total, and\n"
    << eta_granularity_jet << " x " << phi_granularity_jet << " JES bins with " 
    << GetNumberOfJetParametersPerBin() << " free parameters each, or " 
    << GetNumberOfJetParameters() << " in total \n";
}

std::ostream& operator<<( std::ostream& os, const TParameters& cal )
{
  time_t rawtime = time(0);
  struct tm * timeinfo;
  char buffer [80];
  timeinfo = localtime ( &rawtime );
  strftime (buffer,80," # Hamburg Calorimeter Calibration Tool, created %c",timeinfo);
  struct passwd* pw = getpwuid(getuid());	
  
  //double aux[10000], fsum;
  //int npar = cal.GetNumberOfParameters(), iflag=0;
  //cal.fitfunction(npar, aux, fsum, cal.k, iflag);  

  os << buffer << " by " << pw->pw_name << "." << endl 
     << " block CalibParameters = {" << endl
     << "    untracked string  Parametrization    = " << '\"' << cal.p->name() << '\"' <<  endl
     << "    untracked  int32  NTowerParamsPerBin = " << cal.GetNumberOfTowerParametersPerBin() << endl
     << "    untracked int32  NJetParamsPerBin   = " << cal.GetNumberOfJetParametersPerBin() << endl
     << "    untracked int32  NEtaBins           = " << cal.eta_granularity << endl
     << "    untracked int32  NPhiBins           = " << cal.phi_granularity << endl
     << "    untracked bool   EtaSymmetryUsed    = " << cal.eta_symmetry << endl
     << "    untracked double FitChi2            = " << cal.GetFitChi2() << endl
     << " }";
  //--------------------------------------------------------------------
  os << endl
     << " block TowerCalibConstants = {" << endl
     << "    untracked vint32 TowMapEta       = { ";
  //1. ieta
  for (int ieta= -41; ieta<=41; ++ieta){
    if (ieta==0) continue;
    for (unsigned int iphi=1; iphi<=cal.phi_ntwr; ++iphi){
      if (ieta!=-41 || iphi!=1)
        os << ", " << ieta;
      else
        os << ieta;	
    }
  }
  os << " }" << endl << "    untracked vint32 TowMapPhi       = { ";

  //2. iphi
  for (int ieta= -41; ieta<=41; ++ieta){
    if (ieta==0) continue;
    for (unsigned int iphi=1; iphi<=cal.phi_ntwr; ++iphi){
      if (ieta!= -41 || iphi!=1)
        os << ", " << iphi;
      else
        os << iphi;
    }
  }
  os << " }"<< endl;
  os << endl;

  //3. tower calibration constants
  os << "    untracked  int32  TowerParam  = " << cal.GetNumberOfTowerParametersPerBin() << endl;
  for (unsigned int n=0; n<cal.p->nTowerPars(); ++n) {
    os << "    untracked vdouble TowerParam"<< n <<" = { ";
    for (int ieta=-41; ieta<=41; ++ieta){
      if (ieta==0) continue;
      for (unsigned int iphi=1; iphi<=cal.phi_ntwr; ++iphi){
        int index = cal.GetBin(cal.GetEtaBin(ieta),cal.GetPhiBin(iphi));
	if (index<0) continue;
	if (ieta!=-41 || iphi!=1)
          os << ", " << cal.k[index*cal.p->nTowerPars()+n];
	else
          os << cal.k[index*cal.p->nTowerPars()+n];;
      }
    }
    os << " }" << endl; 
  }

  //4. calibration constants errors
  for (unsigned int n=0; n<cal.p->nTowerPars(); ++n) {
    os << "    untracked vdouble TowerError"<<n<<" = { ";
    for (int ieta=-41; ieta<=41; ++ieta){
      if (ieta==0) continue;
      for (unsigned int iphi=1; iphi<=cal.phi_ntwr; ++iphi){
        int index = cal.GetBin(cal.GetEtaBin(ieta),cal.GetPhiBin(iphi));
	if (index<0) continue;
	if (ieta!=-41 || iphi!=1)
          os << ", " << cal.e[index*cal.p->nTowerPars()+n];
	else
          os << cal.e[index*cal.p->nTowerPars()+n];;
      }
    }
    os << " }" << endl; 
  }
  os << " }" << endl; 
  //--------------------------------------------------------------------
  os << endl
     << " block JetCalibConstants = {" << endl
     << "    InputTag Jets    = MyFavoriteJetAlgorithm" << endl
     << "    string CalibJets = \"\" " << endl
     << endl
     << "    untracked vint32 JetMapEta     = { ";
  
  //5. ieta
  for (int ieta= -41; ieta<=41; ++ieta){
    if (ieta==0) continue;
    for (unsigned int iphi=1; iphi<=cal.phi_ntwr; ++iphi){
      if (ieta!=-41 || iphi!=1)
        os << ", " << ieta;
      else
        os << ieta;	
    }
  }
  os << " }" << endl << "    untracked vint32 JetMapPhi     = { ";

  //6. iphi
  for (int ieta= -41; ieta<=41; ++ieta){
    if (ieta==0) continue;
    for (unsigned int iphi=1; iphi<=cal.phi_ntwr; ++iphi){
      if (ieta!= -41 || iphi!=1)
        os << ", " << iphi;
      else
        os << iphi;
    }
  }
  os << " }" << endl;
  os << endl;

  //7. jet calibration constants
  os << "    untracked  int32  JetParam  = " << cal.GetNumberOfJetParametersPerBin() << endl;
  for (unsigned int n=0; n<cal.p->nJetPars(); ++n) {
    os << "    untracked vdouble JetParam"<<n<<" = { ";
    for (int ieta=-41; ieta<=41; ++ieta){
      if (ieta==0) continue;
      for (unsigned int iphi=1; iphi<=cal.phi_ntwr; ++iphi){
	int index = cal.GetJetBin(cal.GetJetEtaBin(ieta),cal.GetJetPhiBin(iphi));
	if (index<0) continue;
	if (ieta!=-41 || iphi!=1)
          os << ", " << cal.k[cal.GetNumberOfTowerParameters() + index*cal.p->nJetPars()+n];
	else
          os << cal.k[cal.GetNumberOfTowerParameters() + index+n];
      }
    }
    os << " }" << endl; 
  }
  //8. calibration constants errors
  for (unsigned int n=0; n<cal.p->nJetPars(); ++n) {
    os << "    untracked vdouble JetError"<<n<<" = { ";
    for (int ieta=-41; ieta<=41; ++ieta){
      if (ieta==0) continue;
      for (unsigned int iphi=1; iphi<=cal.phi_ntwr; ++iphi){
	int index = cal.GetBin(cal.GetJetEtaBin(ieta),cal.GetJetPhiBin(iphi));
	if (index<0) continue;
	if (ieta!=-41 || iphi!=1)
          os << ", " << cal.e[cal.GetNumberOfTowerParameters() + index*cal.p->nJetPars()+n];
	else
          os << cal.e[cal.GetNumberOfTowerParameters() + index+n];
      }
    }
    os << " }" << endl; 
  }
  os << " }" << endl; 
  return os;
}
