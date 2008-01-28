#include <fstream>
#include <pwd.h>
#include <unistd.h>
#include <time.h>

#include "Parameters.h"
#include "ConfigFile.h"


using namespace std;

void TParameters::ReadConfigFile(string const file)
{
  ConfigFile config( file.c_str() );
  
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
  if ( start_values.size()<free_pars_per_bin){
    cerr<< "ERROR: Number of start values and free parameters does not match!"<<endl
        << "       There must be at least " << free_pars_per_bin << " parameters!" << endl;
    exit(2);    
  }
  jet_start_values = bag_of<double>(config.read<string>("jet start values","3.8 -0.7  0.04")); 
  if ( jet_start_values.size()<free_pars_per_bin_jet){
    cerr<< "ERROR: Number of jet start values and free jet parameters does not match!"<<endl
        << "       There must be at least " << free_pars_per_bin_jet << " parameters!" << endl;
    exit(3);
  }
  
  k = new double[GetNumberOfParameters()];
  e = new double[GetNumberOfParameters()];

  for (unsigned int bin=0; bin<eta_granularity*phi_granularity; ++bin){
    for (unsigned int tp=0; tp<free_pars_per_bin; ++tp){
      //step[ bin*free_pars_per_bin + tp ]   = step_sizes[ tp ];
      k[ bin*free_pars_per_bin + tp ] = start_values[ tp ];
      e[ bin*free_pars_per_bin + tp ] = 0.0;
    }
  }
  for (unsigned int bin=0; bin<eta_granularity_jet*phi_granularity_jet; ++bin){
    for (unsigned int jp=0; jp<free_pars_per_bin_jet; ++jp){
      int i = GetNumberOfTowerParameters() + bin*free_pars_per_bin_jet + jp;   
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
  unsigned int pos = result.find(",");
  while (pos!=string::npos) {
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
  std::vector<double> param[free_pars_per_bin], error[free_pars_per_bin];
  std::vector<int> eta_jet, phi_jet;
  std::vector<double> param_jet[free_pars_per_bin_jet], error_jet[free_pars_per_bin_jet];
  int posEqual;
  while (std::getline(file,line)) {

    if (! line.length()) continue;
    if( line.find("#") != string::npos) continue;

    //Read Tower Calibration: ---------------------------------------------------
    if ( line.find("module ccctm = CalibratedCaloTowerMaker") != string::npos ) {
      while (std::getline(file,line)) {
        if( line.find("module") != string::npos) break;
	posEqual=line.find('=');
	name  = trim(line.substr(0,posEqual));
        if( name.find("mapEta") != string::npos) 
  	  eta = bag_of<int>(trim(line.substr(posEqual+1)));
        if( name.find("mapPhi") != string::npos) 
  	  phi = bag_of<int>(trim(line.substr(posEqual+1)));
        for (unsigned i=0; i<free_pars_per_bin; ++i){
	  sprintf(dummy,"TowerParam%d",i);
          if( name.find(dummy) != string::npos) 
  	    param[i] = bag_of<double>(trim(line.substr(posEqual+1)));
	  sprintf(dummy,"TowerError%d",i);
          if( name.find(dummy) != string::npos) 
  	    error[i] = bag_of<double>(trim(line.substr(posEqual+1)));
	}
      }
    }
    //Read Jet Calibration --------------------------------------------------------
    if ( line.find("module cccjm = CalibratedJetMaker") != string::npos ) {
      while (std::getline(file,line)) {
        if( line.find("module") != string::npos) break;
	posEqual=line.find('=');
	name  = trim(line.substr(0,posEqual));
        if( name.find("mapEta") != string::npos) 
  	  eta_jet = bag_of<int>(trim(line.substr(posEqual+1)));
        if( name.find("mapPhi") != string::npos) 
  	  phi_jet = bag_of<int>(trim(line.substr(posEqual+1)));
        for (unsigned i=0; i<free_pars_per_bin_jet; ++i){
	  sprintf(dummy,"JetParam%d",i);
          if( name.find(dummy) != string::npos) 
  	    param_jet[i] = bag_of<double>(trim(line.substr(posEqual+1)));
	  sprintf(dummy,"JetError%d",i);
          if( name.find(dummy) != string::npos) 
  	    error_jet[i] = bag_of<double>(trim(line.substr(posEqual+1)));
	}
      }
    }
  }
  //check if the read calibration is ok:
  bool ok=eta.size()==phi.size();
  for (unsigned i=0; i<free_pars_per_bin; ++i){
    ok *= eta.size()==param[i].size();
    ok *= eta.size()==error[i].size();
  }
  //fill tower parameters and errors:  
  if (ok) {
    for (unsigned i=0; i<eta.size(); ++i){
      int index = GetBin(GetEtaBin(eta[i]),GetPhiBin(phi[i]));
      if (index<0) continue;
      for (unsigned n=0; n<free_pars_per_bin; ++n) {
        k[index*free_pars_per_bin+n] = param[n][i];
        e[index*free_pars_per_bin+n] = error[n][i];
      }
    }
  }

  //check if the read calibration is ok:
  ok=eta_jet.size()==phi_jet.size();
  for (unsigned i=0; i<free_pars_per_bin_jet; ++i){
    ok *= eta_jet.size()==param_jet[i].size();
    ok *= eta_jet.size()==error_jet[i].size();
  }
  //fill Jet parameters and errors:  
  if (ok) {
    for (unsigned i=0; i<eta_jet.size(); ++i){
      int index = GetJetBin(GetJetEtaBin(eta_jet[i]),GetJetPhiBin(phi_jet[i]));
      if (index<0) continue;
      for (unsigned n=0; n<free_pars_per_bin_jet; ++n) {
        k[GetNumberOfTowerParameters() + index*free_pars_per_bin_jet+n] = param_jet[n][i];
        e[GetNumberOfTowerParameters() + index*free_pars_per_bin_jet+n] = error_jet[n][i];
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



std::ostream& operator<<( std::ostream& os, const TParameters& cal )
{
  time_t rawtime = time(0);
  struct tm * timeinfo;
  char buffer [80];
  timeinfo = localtime ( &rawtime );
  strftime (buffer,80," # Hamburg Calorimeter Calibration Tool, created %c",timeinfo);
  struct passwd* pw = getpwuid(getuid());	
  os << buffer << " by " << pw->pw_name << "." << endl 
     << " module calibTowerMaker = CalibTowerMaker {" << endl
     << "    untracked vint32 mapEta       = { ";
  
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
  os << " }" << endl << "    untracked vint32 mapPhi       = { ";

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
  os << "    untracked  int32  TowerParam       = " << cal.GetNumberOfTowerParametersPerBin() << endl;
  for (unsigned int n=0; n<cal.free_pars_per_bin; ++n) {
    os << "    untracked vdouble TowerParam"<< n <<" = { ";
    for (int ieta=-41; ieta<=41; ++ieta){
      if (ieta==0) continue;
      for (unsigned int iphi=1; iphi<=cal.phi_ntwr; ++iphi){
        int index = cal.GetBin(cal.GetEtaBin(ieta),cal.GetPhiBin(iphi));
	if (index<0) continue;
	if (ieta!=-41 || iphi!=1)
          os << ", " << cal.k[index*cal.free_pars_per_bin+n];
	else
          os << cal.k[index*cal.free_pars_per_bin+n];;
      }
    }
    os << " }" << endl; 
  }

  //4. calibration constants errors
  for (unsigned int n=0; n<cal.free_pars_per_bin; ++n) {
    os << "    untracked vdouble TowerError"<<n<<" = { ";
    for (int ieta=-41; ieta<=41; ++ieta){
      if (ieta==0) continue;
      for (unsigned int iphi=1; iphi<=cal.phi_ntwr; ++iphi){
        int index = cal.GetBin(cal.GetEtaBin(ieta),cal.GetPhiBin(iphi));
	if (index<0) continue;
	if (ieta!=-41 || iphi!=1)
          os << ", " << cal.e[index*cal.free_pars_per_bin+n];
	else
          os << cal.e[index*cal.free_pars_per_bin+n];;
      }
    }
    os << " }" << endl; 
  }
  os << " }" << endl; 
  //--------------------------------------------------------------------
  os << endl
     << " module calibJetMaker = CalibJetMaker {" << endl
     << "    InputTag Jets    = " << endl
     << "    string CalibJets = \"\" " << endl
     << endl
     << "    untracked vint32 mapEta     = { ";
  
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
  os << " }" << endl << "    untracked vint32 mapPhi     = { ";

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
  os << "    untracked  int32 JetParam   = " << cal.GetNumberOfJetParametersPerBin() << endl;
  for (unsigned int n=0; n<cal.free_pars_per_bin_jet; ++n) {
    os << "    untracked vdouble  JetParam"<<n<<" = { ";
    for (int ieta=-41; ieta<=41; ++ieta){
      if (ieta==0) continue;
      for (unsigned int iphi=1; iphi<=cal.phi_ntwr; ++iphi){
	int index = cal.GetJetBin(cal.GetJetEtaBin(ieta),cal.GetJetPhiBin(iphi));
	if (index<0) continue;
	if (ieta!=-41 || iphi!=1)
          os << ", " << cal.k[cal.GetNumberOfTowerParameters() + index*cal.free_pars_per_bin_jet+n];
	else
          os << cal.k[cal.GetNumberOfTowerParameters() + index+n];
      }
    }
    os << " }" << endl; 
  }
  //8. calibration constants errors
  for (unsigned int n=0; n<cal.free_pars_per_bin_jet; ++n) {
    os << "    untracked vdouble JetError"<<n<<" = { ";
    for (int ieta=-41; ieta<=41; ++ieta){
      if (ieta==0) continue;
      for (unsigned int iphi=1; iphi<=cal.phi_ntwr; ++iphi){
	int index = cal.GetBin(cal.GetJetEtaBin(ieta),cal.GetJetPhiBin(iphi));
	if (index<0) continue;
	if (ieta!=-41 || iphi!=1)
          os << ", " << cal.e[cal.GetNumberOfTowerParameters() + index*cal.free_pars_per_bin_jet+n];
	else
          os << cal.e[cal.GetNumberOfTowerParameters() + index+n];
      }
    }
    os << " }" << endl; 
  }
  os << " }" << endl; 
  return os;
}
