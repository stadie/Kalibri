#include <fstream>
#include <cassert>
#include <pwd.h>
#include <unistd.h>
#include <ctime>
#include <iomanip>

#include "Parameters.h"
#include "ConfigFile.h"


using namespace std;

TParameters* TParameters::instance = 0;

Parametrization* TParameters::CreateParametrization(const std::string& name) {
  if(name == "StepParametrization") {
    return new StepParametrization();
  } else if(name == "MyParametrization") {
    return new MyParametrization();
  } else if(name == "StepParametrizationEnergy") {
    return new StepParametrizationEnergy();
  } else if(name == "StepEfracParametrization") {
    return new StepEfracParametrization();
  }  else if(name == "JetMETParametrization") {
    return new JetMETParametrization();
  }  else if(name == "TSimpleParametrization") {
    return new SimpleParametrization();
  }  else if(name == "ToyParametrization") {
    return new ToyParametrization();
  } else if(name == "TStepParametrizationEnergy") {
    return new StepParametrizationEnergy();
  } else if(name == "ToyStepParametrizationEnergy") {
    return new ToyStepParametrizationEnergy();
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
  } else if(parclass == "TStepParametersEnergy") {
    parclass = "StepParametrizationEnergy";
  } else if(parclass == "TStepEfracParameters") {
    parclass = "StepEfracParametrization";
  } else if(parclass == "TJetMETParameters") {
    parclass = "JetMETParametrization";
  }  else if(parclass == "TSimpleParameters") {
    parclass = "SimpleParametrization";
  }  else if(parclass == "TToyParameters") {
    parclass = "ToyParametrization";
  }  else if(parclass == "TToyStepParametersEnergy") {
    parclass = "ToyStepParametrizationEnergy";
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
      int i = GetNumberOfTowerParameters() + bin*p->nJetPars() + jp;   
      k[i] = jet_start_values[jp];
      e[i] = 0.0;
    }
  }

  // read predefined calibration contants from cfi
  // or txt file depending on the ending of the name
  cout << "Reading calibration from file '" << input_calibration << endl;
  if(!input_calibration.empty()){
    if( !input_calibration.substr(input_calibration.rfind(".")+1).compare("cfi") ){
      Read_CalibrationCfi(input_calibration); 
    } 
    else if( !input_calibration.substr(input_calibration.rfind(".")+1).compare("txt") ){
      cout << "call function" << endl;
      Read_CalibrationTxt(input_calibration); 
    }
    else{
      cout << "Error: unknown file format: '" 
	   << input_calibration.substr(input_calibration.rfind("."))
	   << "'"<< endl;  
    }
  }
}

std::string TParameters::trim(std::string const& source, char const* delims) 
{
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

void TParameters::Read_CalibrationTxt(std::string const& configFile)
{
  // ---------------------------------------------------------------
  // fills start parameters for fit when read from txt file; expects 
  // 72 lines for 72 bins in phi for each eta bin ranging from -41
  // to 41 (skipping the 0) and the following parameter format:
  // maxEta minEta nPar towerParameters jetParameters separated by
  // blanks
  // ---------------------------------------------------------------
  std::ifstream file(configFile.c_str());
  std::string line; // buffer line

  int      etaBin=-41;
  unsigned phiBin=  0;
  unsigned iLines=  0;
  while( std::getline(file,line) ){
    // determine phi bin on the fly
    phiBin=(iLines%72)+1;      // phi counts from 1...72 for each eta bin
    ++iLines;
    // determine eta bin on the fly
    if(iLines%72==0) ++etaBin; // increas etaValue by 1 each 72 lines...
    if(etaBin   ==0) ++etaBin; // and don't forget to skip the 0

    //cout << "etaBin: " << etaBin << " :: " << "phiBin: " << phiBin << endl;

    // buffers for input parameters
    unsigned nPar=0; //this is not needed but read out for control reasons 
    double etaMax=0; //this is not needed but read out for control reasons  
    double etaMin=0; //this is not needed but read out for control reasons 
    std::vector<double> twrPars, jetPars;
    unsigned entry=0; // controls which parameter is to filled
    while( line.length()>line.substr(0, line.find(" ")).size() ){
      if( 0<line.find(" ")){
   	// extract value
	switch(++entry){
	case 1 : etaMin = std::atof( line.substr(0, line.find(" ")).c_str() ); 
	  break;
	case 2 : etaMax = std::atof( line.substr(0, line.find(" ")).c_str() ); 
	  break; 
	case 3 : nPar   = std::atoi( line.substr(0, line.find(" ")).c_str() ); 
	  break; 
	default:
	  if((entry-3)<=p->nTowerPars()){
	    twrPars.push_back(std::atof( line.substr(0, line.find(" ")).c_str() ));
	  }
	  else{
	    jetPars.push_back(std::atof( line.substr(0, line.find(" ")).c_str() ));
	  }
	  break;
	}
	// cut string
	line = line.substr(line.find(" "));
      }
      else{
	//cut string
	if(line.find(" ")<std::string::npos){
	  line = line.substr(line.find(" ")+1);
	}
      }
    }
    // catch last character
    jetPars.push_back(std::atof( line.substr(0, line.find(" ")).c_str() ));

    // fill parameters
    int towerIdx = GetBin(GetEtaBin(etaBin),GetPhiBin(phiBin));
    if( towerIdx<0 ) continue;
    for (unsigned n=0; n< p->nTowerPars(); ++n) {
      k[towerIdx*p->nTowerPars()+n] = twrPars[n];
      //e[towerIdx*p->nTowerPars()+n] = NOT_READ_OUT;
    }
    int jetIdx = GetJetBin(GetJetEtaBin(etaBin),GetJetPhiBin(phiBin));
    if( jetIdx<0 ) continue;
    for (unsigned n=0; n<p->nJetPars(); ++n) {
      k[GetNumberOfTowerParameters()+jetIdx*p->nJetPars()+n] = jetPars[n];
      //e[GetNumberOfTowerParameters()+jetIdx*p->nJetPars()+n] = NOT_READ_OUT;
    }
  }
}

void TParameters::Read_CalibrationCfi(std::string const& configFile)
{
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
	//std::cout << name << ".\n";
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

int TParameters::GetEtaBin(int eta_id, int etagranu, int phigranu, bool etasym) const
{
//This function knows the number of wanted eta-bins and returns 
//in which eta-bin the tower with eta-ID "eta_id" is located.
  //Case 1 bin:
//cout << "eta="<<eta_id<<", etagranu:"<< etagranu<< ", eta_ntwr_used:"<< eta_ntwr_used<<endl;
  if (etagranu<=1) return 0;
  if (etagranu==2) return (eta_id < 0) ? 0 : 1;

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
  if (!etasym){
    if (etagranu==82) return index;
    else if (etagranu==42) return ta_42[index];
    else if (etagranu==22) return ta_22[index];
    else if (etagranu==10) return ta_10[index];
    else if (etagranu==6) return ta_6[index];
  } else {
    if (etagranu==41) return abs(eta_id)-1;
    else if (etagranu==21) return ts_21[abs(eta_id)-1];
    else if (etagranu==11) return ts_11[abs(eta_id)-1];
    else if (etagranu== 5) return ts_5[ abs(eta_id)-1];
    else if (etagranu== 3) return ts_3[ abs(eta_id)-1];
  }
  
  //Default value, should never be returned!
  return -4;
}

float TParameters::EtaEdge(int const etaBin, bool lowerEdge)
{
  // return eta bin - eta edge mappting
  switch(etaBin){
  case -41: return (lowerEdge ? -5.191 : -4.889); break;
  case -40: return (lowerEdge ? -4.889 : -4.716); break;
  case -39: return (lowerEdge ? -4.716 : -4.538); break;
  case -38: return (lowerEdge ? -4.538 : -4.363); break;
  case -37: return (lowerEdge ? -4.363 : -4.191); break;
  case -36: return (lowerEdge ? -4.191 : -4.013); break;
  case -35: return (lowerEdge ? -4.013 : -3.839); break;
  case -34: return (lowerEdge ? -3.839 : -3.664); break;
  case -33: return (lowerEdge ? -3.664 : -3.489); break;
  case -32: return (lowerEdge ? -3.489 : -3.314); break;
  case -31: return (lowerEdge ? -3.314 : -3.139); break;
  case -30: return (lowerEdge ? -3.139 : -2.964); break;
  case -29: return (lowerEdge ? -2.964 : -2.853); break; 
  case -28: return (lowerEdge ? -2.853 :  -2.65); break;
  case -27: return (lowerEdge ?  -2.65 :   -2.5); break;
  case -26: return (lowerEdge ?   -2.5 : -2.322); break;
  case -25: return (lowerEdge ? -2.322 : -2.172); break;
  case -24: return (lowerEdge ? -2.172 : -2.043); break;
  case -23: return (lowerEdge ? -2.043 :  -1.93); break;
  case -22: return (lowerEdge ?  -1.93 :  -1.83); break;
  case -21: return (lowerEdge ?  -1.83 :  -1.74); break;
  case -20: return (lowerEdge ?  -1.74 : -1.653); break;
  case -19: return (lowerEdge ? -1.653 : -1.566); break;
  case -18: return (lowerEdge ? -1.566 : -1.479); break;
  case -17: return (lowerEdge ? -1.479 : -1.392); break;
  case -16: return (lowerEdge ? -1.392 : -1.305); break;
  case -15: return (lowerEdge ? -1.305 : -1.218); break;
  case -14: return (lowerEdge ? -1.218 : -1.131); break;
  case -13: return (lowerEdge ? -1.131 : -1.044); break;
  case -12: return (lowerEdge ? -1.044 : -0.957); break;
  case -11: return (lowerEdge ? -0.957 : -0.879); break;
  case -10: return (lowerEdge ? -0.879 : -0.783); break;
  case  -9: return (lowerEdge ? -0.783 : -0.696); break;
  case  -8: return (lowerEdge ? -0.696 : -0.609); break;
  case  -7: return (lowerEdge ? -0.609 : -0.522); break;
  case  -6: return (lowerEdge ? -0.522 : -0.435); break;
  case  -5: return (lowerEdge ? -0.435 : -0.348); break;
  case  -4: return (lowerEdge ? -0.348 : -0.261); break;
  case  -3: return (lowerEdge ? -0.261 : -0.174); break;
  case  -2: return (lowerEdge ? -0.174 : -0.087); break;
  case  -1: return (lowerEdge ? -0.087 :      0); break;
  case  +1: return (lowerEdge ?      0 :  0.087); break;
  case  +2: return (lowerEdge ?  0.087 :  0.174); break;
  case  +3: return (lowerEdge ?  0.174 :  0.261); break;
  case  +4: return (lowerEdge ?  0.261 :  0.348); break;
  case  +5: return (lowerEdge ?  0.348 :  0.435); break;
  case  +6: return (lowerEdge ?  0.435 :  0.522); break;
  case  +7: return (lowerEdge ?  0.522 :  0.609); break;
  case  +8: return (lowerEdge ?  0.609 :  0.696); break;
  case  +9: return (lowerEdge ?  0.696 :  0.783); break;
  case +10: return (lowerEdge ?  0.783 :  0.879); break;
  case +11: return (lowerEdge ?  0.879 :  0.957); break;
  case +12: return (lowerEdge ?  0.957 :  1.044); break;
  case +13: return (lowerEdge ?  1.044 :  1.131); break;
  case +14: return (lowerEdge ?  1.131 :  1.218); break;
  case +15: return (lowerEdge ?  1.218 :  1.305); break;
  case +16: return (lowerEdge ?  1.305 :  1.392); break;
  case +17: return (lowerEdge ?  1.392 :  1.479); break;
  case +18: return (lowerEdge ?  1.479 :  1.566); break;
  case +19: return (lowerEdge ?  1.566 :  1.653); break;
  case +20: return (lowerEdge ?  1.653 :   1.74); break;
  case +21: return (lowerEdge ?   1.74 :   1.83); break;
  case +22: return (lowerEdge ?   1.83 :   1.93); break;
  case +23: return (lowerEdge ?   1.93 :  2.043); break;
  case +24: return (lowerEdge ?  2.043 :  2.172); break;
  case +25: return (lowerEdge ?  2.172 :  2.322); break;
  case +26: return (lowerEdge ?  2.322 :    2.5); break;
  case +27: return (lowerEdge ?    2.5 :   2.65); break;
  case +28: return (lowerEdge ?   2.65 :  2.853); break;
  case +29: return (lowerEdge ?  2.853 :  2.964); break;
  case +30: return (lowerEdge ?  2.964 :  3.139); break;
  case +31: return (lowerEdge ?  3.139 :  3.314); break;
  case +32: return (lowerEdge ?  3.314 :  3.489); break;
  case +33: return (lowerEdge ?  3.489 :  3.664); break;
  case +34: return (lowerEdge ?  3.664 :  3.839); break;
  case +35: return (lowerEdge ?  3.839 :  4.013); break;
  case +36: return (lowerEdge ?  4.013 :  4.191); break;
  case +37: return (lowerEdge ?  4.191 :  4.363); break;
  case +38: return (lowerEdge ?  4.363 :  4.538); break;
  case +39: return (lowerEdge ?  4.538 :  4.716); break;
  case +40: return (lowerEdge ?  4.716 :  4.889); break;
  case +41: return (lowerEdge ?  4.889 :  5.191); break;
    //something went wrong;
  default : return -1; break;
  }
}

int TParameters::GetPhiBin(int phi_id, int phigranu) const
//This function knows the number of wanted phi-bins and returns 
//in which phi-bin the tower with eta-ID "phi_id" is located.
{
  return (phi_id-1)*phigranu/phi_ntwr;
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

void TParameters::Write_CalibrationTxt(const char* name)
{
  cout << "Writing calibration to file '" << name << "'" << endl;

  ofstream file(name, ofstream::binary);
  for (int ieta= -41; ieta<=41; ++ieta){
    if (ieta==0) continue;
    for (unsigned iphi=1; iphi<=phi_ntwr; ++iphi){
      int towerIdx = GetBin(GetEtaBin(ieta),GetPhiBin(iphi));
      int jetIdx = GetJetBin(GetJetEtaBin(ieta),GetJetPhiBin(iphi));
      if(towerIdx<0 || jetIdx<0) continue;
      // write: lower eta | upper eta | nparameters, for
      // each eta id of the tower and n times for n phi bins
      file << std::setw(10) << EtaLowerEdge(ieta) 
	   << std::setw(10) << EtaUpperEdge(ieta)  
	   << std::setw(10) << p->nTowerPars()+p->nJetPars();
      // write: each tower parameter
      for(unsigned itower=0; itower<p->nTowerPars(); ++itower){
	file << std::setw(8) << std::setprecision(4) << k[towerIdx*p->nTowerPars()+itower];
      }
      // write: each jet parameter
      for(unsigned ijet=0; ijet<p->nJetPars(); ++ijet){
	file << std::setw(8) << std::setprecision(4) << k[GetNumberOfTowerParameters()+jetIdx*p->nJetPars()+ijet];
      }
      // complete line
      file << std::endl;
    }
  }
  file.close();
}

void TParameters::Write_CalibrationCfi(const char* name)
{
  cout << "Writing calibration to file '" << name << "'" << endl;
  ofstream file(name, ofstream::binary);
  
  time_t rawtime = time(0);
  struct tm * timeinfo;
  char buffer [80];
  timeinfo = localtime ( &rawtime );
  strftime (buffer,80," # Hamburg Calorimeter Calibration Tool, created %c",timeinfo);
  struct passwd* pw = getpwuid(getuid());	
  
  //double aux[10000], fsum;
  //int npar = GetNumberOfParameters(), iflag=0;
  //fitfunction(npar, aux, fsum, k, iflag);  

  file << buffer << " by " << pw->pw_name << "." << endl 
       << " block CalibParameters = {" << endl
       << "    untracked string  Parametrization    = " << '\"' << p->name() << '\"' <<  endl
       << "   untracked  int32  NTowerParamsPerBin = " << GetNumberOfTowerParametersPerBin() << endl
       << "    untracked int32  NJetParamsPerBin   = " << GetNumberOfJetParametersPerBin() << endl
       << "    untracked int32  NEtaBins           = " << eta_granularity << endl
       << "    untracked int32  NPhiBins           = " << phi_granularity << endl
       << "    untracked bool   EtaSymmetryUsed    = " << eta_symmetry << endl
       << "    untracked double FitChi2            = " << GetFitChi2() << endl
       << " }";
  //--------------------------------------------------------------------
  file << endl
       << " block TowerCalibConstants = {" << endl
       << "    untracked vint32 TowMapEta       = { ";
  //1. ieta
  for (int ieta= -41; ieta<=41; ++ieta){
    if (ieta==0) continue;
    for (unsigned int iphi=1; iphi<=phi_ntwr; ++iphi){
      if (ieta!=-41 || iphi!=1)
        file << ", " << ieta;
      else
        file << ieta;	
    }
  }
  file << " }" << endl << "    untracked vint32 TowMapPhi       = { ";

  //2. iphi
  for (int ieta= -41; ieta<=41; ++ieta){
    if (ieta==0) continue;
    for (unsigned int iphi=1; iphi<=phi_ntwr; ++iphi){
      if (ieta!= -41 || iphi!=1)
        file << ", " << iphi;
      else
        file << iphi;
    }
  }
  file << " }"<< endl;
  file << endl;

  //3. tower calibration constants
  file << "    untracked  int32  TowerParam  = " << GetNumberOfTowerParametersPerBin() << endl;
  for (unsigned int n=0; n<p->nTowerPars(); ++n) {
    file << "    untracked vdouble TowerParam"<< n <<" = { ";
    for (int ieta=-41; ieta<=41; ++ieta){
      if (ieta==0) continue;
      for (unsigned int iphi=1; iphi<=phi_ntwr; ++iphi){
        int index = GetBin(GetEtaBin(ieta),GetPhiBin(iphi));
	if (index<0) continue;
	if (ieta!=-41 || iphi!=1)
          file << ", " << k[index*p->nTowerPars()+n];
	else
          file << k[index*p->nTowerPars()+n];;
      }
    }
    file << " }" << endl; 
  }

  //4. calibration constants errors
  for (unsigned int n=0; n<p->nTowerPars(); ++n) {
    file << "    untracked vdouble TowerError"<<n<<" = { ";
    for (int ieta=-41; ieta<=41; ++ieta){
      if (ieta==0) continue;
      for (unsigned int iphi=1; iphi<=phi_ntwr; ++iphi){
        int index = GetBin(GetEtaBin(ieta),GetPhiBin(iphi));
	if (index<0) continue;
	if (ieta!=-41 || iphi!=1)
          file << ", " << e[index*p->nTowerPars()+n];
	else
          file << e[index*p->nTowerPars()+n];;
      }
    }
    file << " }" << endl; 
  }
  file << " }" << endl; 
  //--------------------------------------------------------------------
  file << endl
       << " block JetCalibConstants = {" << endl
       << "    InputTag Jets    = MyFavoriteJetAlgorithm" << endl
       << "    string CalibJets = \"\" " << endl
       << endl
       << "    untracked vint32 JetMapEta     = { ";
  
  //5. ieta
  for (int ieta= -41; ieta<=41; ++ieta){
    if (ieta==0) continue;
    for (unsigned int iphi=1; iphi<=phi_ntwr; ++iphi){
      if (ieta!=-41 || iphi!=1)
        file << ", " << ieta;
      else
        file << ieta;	
    }
  }
  file << " }" << endl << "    untracked vint32 JetMapPhi     = { ";
  
  //6. iphi
  for (int ieta= -41; ieta<=41; ++ieta){
    if (ieta==0) continue;
    for (unsigned int iphi=1; iphi<=phi_ntwr; ++iphi){
      if (ieta!= -41 || iphi!=1)
        file << ", " << iphi;
      else
        file << iphi;
    }
  }
  file << " }" << endl;
  file << endl;

  //7. jet calibration constants
  file << "    untracked  int32  JetParam  = " << GetNumberOfJetParametersPerBin() << endl;
  for (unsigned int n=0; n<p->nJetPars(); ++n) {
    file << "    untracked vdouble JetParam"<<n<<" = { ";
    for (int ieta=-41; ieta<=41; ++ieta){
      if (ieta==0) continue;
      for (unsigned int iphi=1; iphi<=phi_ntwr; ++iphi){
	int index = GetJetBin(GetJetEtaBin(ieta),GetJetPhiBin(iphi));
	if (index<0) continue;
	if (ieta!=-41 || iphi!=1)
          file << ", " << k[GetNumberOfTowerParameters()+index*p->nJetPars()+n];
	else
          file << k[GetNumberOfTowerParameters()+index*p->nJetPars()+n];
      }
    }
    file << " }" << endl; 
  }
  //8. calibration constants errors
  for (unsigned int n=0; n<p->nJetPars(); ++n) {
    file << "    untracked vdouble JetError"<<n<<" = { ";
    for (int ieta=-41; ieta<=41; ++ieta){
      if (ieta==0) continue;
      for (unsigned int iphi=1; iphi<=phi_ntwr; ++iphi){
	int index = GetBin(GetJetEtaBin(ieta),GetJetPhiBin(iphi));
	if (index<0) continue;
	if (ieta!=-41 || iphi!=1)
          file << ", " << e[GetNumberOfTowerParameters() + index*p->nJetPars()+n];
	else
          file << e[GetNumberOfTowerParameters() + index*p->nJetPars()+n];
      }
    }
    file << " }" << endl; 
  }
  file << " }" << endl; 
}
