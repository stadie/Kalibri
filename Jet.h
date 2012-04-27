//!    \brief Class for basic jets 
//!
//!    \author Hartmut Stadie
//!
//!    \date 2008/12/14
//!
//!    $Id: Jet.h,v 1.49 2012/02/06 22:29:37 kirschen Exp $
#ifndef JET_H
#define JET_H

#include "CalibData.h"
#include "Function.h"
#include "Parameters.h"

#include "gsl/gsl_errno.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_roots.h"
   
#include <cstdlib>
  
class CorFactors;


class Jet : public Measurement
{
 public:
  class JetIndex {
  public:
    JetIndex(unsigned int idx, float pt) : idx_(idx), pt_(pt) {};
    const unsigned int idx_;
    const float pt_;
    // For sorting jets in pt
    static bool ptGreaterThan(const JetIndex *idx1, const JetIndex *idx2) {
      // check for 0
      if(idx1 == 0) {
	return idx2 != 0;
      } else if (idx2 == 0) {
	return false;
      } else {
	return idx1->pt_ > idx2->pt_;
      }
    }
  };

  //! For sorting jets in calo pt
  static bool caloPtGreaterThan(const Jet *j1, const Jet *j2) {
    // check for 0
    if (j1 == 0) {
      return j2 != 0;
    } else if (j2 == 0) {
      return false;
    } else {
      return j1->pt() > j2->pt();
    }
  }


  //!  \brief Jet flavor
  //!
  //!  The possible flavors are
  //!  - 0: Gluon
  //!  - 1: u, d, or s quark
  //!  - 2: c quark
  //!  - 3: b quark
  enum Flavor{ unknown = -1, gluon=0, uds=1, c=2, b=3 };

  //! return flavor for pdg id
  static Flavor flavorFromPDG(int pdg) {
    if(pdg == 21) return gluon;
    if(pdg == 0) return unknown;
    unsigned int id = std::abs(pdg);
    if(id < 4) return uds;
    if(id == 4) return c;
    if(id == 5) return b;
    return unknown;
  }

 public:
  Jet(float Et, float EmEt, float HadEt ,float OutEt, float E,
      float eta,float phi, float phiphi, float etaeta, Flavor flavor, 
      float fCH, float fNH, float fPH, float fEL, float fHFEm, float fHFHad, 
      float genPt, float dR, CorFactors* corFactors, const Function& f,
      float (*errfunc)(const float *x, const Measurement *xorig, float err), 
      const Function& gf); 
  virtual ~Jet();

  float Et()     const {return Measurement::pt;}                 //!< Return transverse energy Et
  float pt()     const {return Measurement::pt;}                 //!< Return transverse energy Et
  float EmEt()   const {return EMF;}                //!< Return Et from the ECAL part of the towers
  float HadEt()  const {return HadF;}               //!< Return Et from the HCAL part of the towers
  float OutEt()  const {return OutF;}               //!< Return Et from the HOut part of the towers
  float E()      const {return Measurement::E;}    //!< Return energy
  float emf()    const {return EMF/(EMF+HadF);}    //!< Return fraction of ECAL energy
  float hadf()    const {return HadF/(EMF+HadF);}  //!< Return fraction of HCAL energy
  float eta()    const {return Measurement::eta;}  //!< Return pseudorapidity
  float phi()    const {return Measurement::phi;}  //!< Return azimuthal angle
  float momentPhiPhi() const {return Measurement::phiphi;}  //!< Return phi-phi moment (width of jet in phi)
  float momentEtaEta() const {return Measurement::etaeta;}  //!< Return eta-eta moment (width of jet in eta)
  float meanMoment() const {return 0.5 * (Measurement::phiphi + Measurement::etaeta);}  //!< Return mean moment (width of jet)
  float fCH()    const {return fCH_;}                //!< Return charged hadron fraction
  float fNH()    const {return fNH_;}                //!< Return neutral hadron fraction
  float fPH()    const {return fPH_;}                //!< Return photon fraction
  float fEL()    const {return fEL_;}                //!< Return electron fraction
  float fHFEm()  const {return fHFEm_;}              //!< Return fraction of HF em
  float fHFHad() const {return fHFHad_;}	     //!< Return fraction of HF had
  Flavor flavor() const {return flavor_;}       //!< Return jet flavor
  float genPt()  const {return genPt_;}        //!< Return Pt for corresponding GenJet 
  float dR() const {return dR_;}               //!< \f$ \Delta R \f$ between jet and genjet
  const CorFactors& corFactors() const { return *corFactors_;}
  void updateCorFactors(CorFactors *cor);  
  //! Correct measurement by \p L1
  void correctL1();
  //! Correct measurement by product \p L1*L2*L3
  void correctToL3();
  //! Correct measurement by product \p L2*L3
  void correctL2L3();
  //! Correct measurement by L1*l2*L3*Residual
  void correctToLRes();

  //!  \brief Change address of parameters covered by this jet
  //!  \sa Parameters
  // ---------------------------------------------------------
  virtual void setParameters(Parameters* param); 
  virtual float correctedEt() const { return correctedEt(Et()); }
  virtual float correctedEt(float Et, bool fast = false) const;
  float expectedEt(float truth, float start, float& error,
		    bool fast = false);

  //!  \brief Calculate error from original measurement
  //!
  //!  The error is calculated using the error function
  //!  errf and Measurement::pt, the pt of the original
  //!  measurement.
  //!  
  //!  \return Error of original measurement
  // ---------------------------------------------------------
  virtual float error() const {return errf_(&(Measurement::pt),this,error_);}


  void setError(float error) { error_ = error;}

  //!  \brief Calculate error from given Et
  //!
  //!  The error is calculated using the error function
  //!  errf and a given Et.
  //!
  //!  \param et Jet Et
  //!  \return Error from jet Et
  // ---------------------------------------------------------
  virtual float expectedError(float et) const { return errf_(&et,this,error_);}

  //!  \brief Get number of parameters of this jet
  //!
  //!  Returns the number of parameters of the
  //!  local and global jet correction functions
  //!  (see Parametrization).
  //!
  //!  \return Number of parameters of this jet
  // ---------------------------------------------------------
  virtual int nPar() const {return f_->nPars() + gf_->nPars();}

  //!  \brief Represents a corrected jet if one parameter of
  //!         the correction function is varied (for derivative
  //!         calculation)
  //!
  //!  Stores the corrected Et and error if a parameter is
  //!  varied by +/- eps for derivative calculation.
  // ---------------------------------------------------------
  virtual const Parameters::VariationColl& varyPars(const double* eps, float Et, float start);
  virtual const Parameters::VariationColl& varyParsDirectly(const double* eps, bool computeDeriv = true, float Et = 0);
  
  void print();                       //!< Print some jet members
  static void printInversionStats();  //!< Print some info on inversion

  int parIndex() const { return f_->parIndex(); }

  virtual Jet* clone() const { return new Jet(*this);} //!< Clone this jet
  void setGlobalFunction(const Function& ngf) { gf_ = &ngf;} //!< Set global correction function, needed for constraints

  const Function* f() const { return f_;}
  const Function* gf() const { return gf_;}
  float (*errFunc())(const float *x, const Measurement *xorig, float err) {
    return errf_;}
 protected:
  virtual float expectedEt(float truth, float start, bool fast = false);
  Jet(const Jet&j); //!< disallow copies!

 private: 
  Flavor flavor_;           //!< The jet's Flavor
  float fCH_;               //!< Fraction of charged hadrons 
  float fNH_;		    //!< Fraction of neutral hadrons
  float fPH_;		    //!< Fraction of photons
  float fEL_;               //!< Fraction of electrons
  float fHFEm_;             //!< Fraction of HF em
  float fHFHad_;	    //!< Fraction of HF had

  float genPt_;            //!< The genjet pt
  float dR_;               //!< \f$ \Delta R \f$ between jet and genjet 
  float error_;                //!< Stores error for constant error mode
  const CorFactors* corFactors_;   //!< The correction factors
  const Function*  f_;            //!< Jet correction function
  const Function*  gf_;           //!< Global jet correction function
  float    (*errf_)(const float *x, const Measurement *xorig, float err);   //!< Error function

  bool      secant(double truth, double& x1, double& x2, double eps);
  bool      falseposition(double truth, double& x1, double& x2, double eps);

  static long long ncalls_;        //!< Number of calls of inversion methods secant and falseposition
  static long long ntries_;        //!< Number of tries in iteration during inversion
  static long long nfails_;        //!< Number of failed tries during inversion
  static long long nwarns_;        //!< Number of warnings during inversion

  struct rf_par {
    double y_;
    const Jet* jet_;
  rf_par(double y, const Jet *jet) : y_(y), jet_(jet) {}
  };
  static double rf(double x, void* params) {
    rf_par* p = (rf_par*)params;
    return p->y_ - p->jet_->correctedEt(x,true);
  };
 protected:
  Parameters* parameters_;
};

#endif
