// $Id: SmearPhotonJet.h,v 1.1 2009/06/10 14:21:13 mschrode Exp $

#ifndef SmearPhotonJet_h
#define SmearPhotonJet_h

#include "CalibData.h"
#include "SmearData.h"


//!  \brief Photon-jet data for jetsmearing method
//!  \author Matthias Schroeder
//!  \date Tue Jun  9 18:23:44 CEST 2009
//!  $Id: SmearPhotonJet.h,v 1.1 2009/06/10 14:21:13 mschrode Exp $
// --------------------------------------------------
class SmearPhotonJet : public SmearData {
 public:
  SmearPhotonJet(TMeasurement * mess, double photonpt, double weight, const Function& respPDF)
    : SmearData(TypeSmearPhotonJet,mess,photonpt,weight,respPDF) {};
  ~SmearPhotonJet() {};

  virtual double chi2() const;
  virtual double chi2_fast(double * temp_derivative1, double * temp_derivative2, double const epsilon) const;
  virtual void PrintInitStats() const {};  //!< No functionality yet
};

#endif
