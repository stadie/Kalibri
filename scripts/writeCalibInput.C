// $Id: $
//
// This script reads L2L3 calibration constants
// from a txt file with the same format as in
// CondFormats database. It writes the constants
// to a txt file in a format readable to caliber.

#include <cassert>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>


typedef std::vector<double> ParsL3;
typedef std::vector<double>::iterator ParsL3It;
typedef std::vector< std::vector<double> > ParsL2;
typedef std::vector< std::vector<double> >::iterator ParsL2It;


ParsL3 readL3Correction(const std::string& filename);
ParsL2 readL2Correction(const std::string& filename);
double etaEdge(int const etaBin, bool lowerEdge);



// ----------------------------------------------------------------  
void writeCalibInput(const std::string& corrl2filename,
		     const std::string& corrl3filename,
		     const std::string& outfilename) {

  ParsL2 parL2 = readL2Correction(corrl2filename);
  ParsL3 parL3 = readL3Correction(corrl3filename);

  std::cout << "Writing parameters into kalibri format\n";

  ofstream file(outfilename.c_str(), ofstream::binary);
  int iEta = -41;
  ParsL2It bin = parL2.begin();
  for(; bin != parL2.end(); bin++, iEta++) {
    if( iEta == 0 ) iEta++;

    file << std::setw(8) << std::setprecision(4) << etaEdge(iEta,true);
    file << std::setw(8) << std::setprecision(4) << etaEdge(iEta,false);
    file << std::setw(8) << std::setprecision(0) << 9;
    file << std::setw(8) << std::setprecision(1) << 4;
    file << std::setw(8) << std::setprecision(1) << 2000;
    for(size_t i = 0; i < bin->size(); i++) {
      file << std::setw(12) << std::setprecision(6) << bin->at(i);
    }
    for(size_t i = 0; i < parL3.size(); i++) {
      file << std::setw(14) << std::setprecision(6) << parL3.at(i);
    }      
    file << std::endl;
  }

  std::cout << "Done. Wrote file " << outfilename << std::endl;
}



//  Read parameters of L3 correction from
//  txt file in CondDB format i.e.
//  etaMin etaMax nPar EtMin EtMax Par1 Par2 Par3 Par4
// ----------------------------------------------------------------  
ParsL3 readL3Correction(const std::string& filename) {
  std::cout << "Reading parameters for L3 correction\n";

  std::ifstream file;
  file.open(filename.c_str());

  double val  = -1.;

  file >> val;                       // Eta min
  file >> val;                       // Eta max
  file >> val;                       // Number of values following
  int n = static_cast<int>(val - 2); // Number of L3 parameters
  ParsL3 par(n,0.);                  // The L3 parameters
  file >> val;                       // Et min
  file >> val;                       // Et max
  for(int i = 0; i < n; i++) {       // Store L3 parameters
    file >> val;
    par.at(i) = val;
  }
  file.close();

  return par;
} 



//  Read parameters of L2 correction from
//  txt file in CondDB format i.e.
//  etaMin etaMax nPar EtMin EtMax Par1 Par2 Par3 0 0 0
// ----------------------------------------------------------------  
ParsL2 readL2Correction(const std::string& filename) {
  std::cout << "Reading parameters for L2 correction\n";

  ParsL2 parL2;

  std::ifstream file;
  file.open(filename.c_str());

  int    etaBin = -41;
  double etaMin = 0.;
  double etaMax = 0.;
  int      nPar = 0;
  double    val = -1.;
  while( !file.eof() && etaBin < 42 ) {
    if( etaBin == 0 ) etaBin++;

    file >> etaMin;                      // Eta min
    file >> etaMax;                      // Eta max
    val = 0.;
    file >> val;                         // Number of values following
    if( val != 0 ) {                     // Avoid reading of empty last line
      nPar = static_cast<int>(val - 5);  // Number of L2 parameters
      std::vector<double> par(nPar,0.);  // The L2 parameters
      file >> val;                       // Et min
      file >> val;                       // Et max
      for(int i = 0; i < nPar; i++) {    // Store L2 parameters
	file >> val;
	par.at(i) = val;
      }
      for(int i = 0; i < 3; i++) {       // Last 3 parameters are 0
	file >> val;
      }
      
      // In case some eta bin is missing,
      // add default parameters
      while( etaBin < 42 && 
	     etaMin != etaEdge(etaBin,true) && 
	     etaMax != etaEdge(etaBin,false)    ) {
	std::cout << "  Warning: No parameters for eta bin " << etaBin;
	std::cout << "; using default parameters instead.\n";
	
	std::vector<double> pardef(nPar,0.);
	for(int i = 0; i < nPar; i++) {
	  if( i == 0 ) pardef.at(i) = 1.;
	  else         pardef.at(i) = 0.;
	}
	parL2.push_back(pardef);
	
	etaBin++;
      }
      
      if( etaBin < 42 ) {
	parL2.push_back(par);
      }
     etaBin++;
    }
  }
  file.close();

  // In case last eta bins are missing,
  // add default parameters  
  while( etaBin < 42 ) {
    std::cout << "  Warning: No parameters for eta bin " << etaBin;
    std::cout << "; using default parameters instead.\n";
	
    nPar = 3;
    std::vector<double> pardef(nPar,0.);
    for(int i = 0; i < nPar; i++) {
      if( i == 0 ) pardef.at(i) = 1.;
      else         pardef.at(i) = 0.;
    }
    parL2.push_back(pardef);
    
    etaBin++;
  }
  

  return parL2;
} 



// ----------------------------------------------------------------  
double etaEdge(int const etaBin, bool lowerEdge) {
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
