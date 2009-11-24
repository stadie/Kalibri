#ifndef JET_H
#define JET_H

#include"CalibData.h"
#include "Function.h"


//!    \brief Class for basic jets 
//!
//!    \author Hartmut Stadie
//!
//!    \date 2008/12/14
//!
//!    $Id: Jet.h,v 1.23 2009/10/26 21:00:36 mschrode Exp $


#include "gsl/gsl_errno.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_roots.h"
     
class Jet : public Measurement
{
 public:
  //!  \brief Jet flavor
  //!
  //!  The possible flavors are
  //!  - 0: Gluon
  //!  - 1: u, d, or s quark
  //!  - 2: c quark
  //!  - 3: b quark
  enum Flavor{ gluon=0, uds=1, c=2, b=3 };

  //!  \brief   Container class for jet correction factors
  class CorFactors
  {
  public :
  CorFactors(double L1=1.0, double L2=1.0, double L3=1.0, double L4=1.0, 
	     double L5=1.0, double JPT=1.0, double JPTL2L3=1.0) 
    : l1_(L1),l2_(L2),l3_(L3),l4_(L4),l5_(L5),jpt_(JPT),jptL2L3_(JPTL2L3) {};
    double getL1()  const { return l1_; }    //!< Return L1 correction factor (zero-suppression)
    double getL2()  const { return l2_; }    //!< Return L2 correction factor (relative, in eta)
    double getL3()  const { return l3_; }    //!< Return L3 correction factor (absolute, in pt)
    double getL4()  const { return l4_; }    //!< Return L4 correction factor (electromagnetic fraction)
    double getL5()  const { return l5_; }    //!< Return L5 correction factor (flavor)
    double getJPT() const { return jpt_; }   //!< Return Jet+Track correction factor
    double getL2L3() const { return l2_*l3_; }   //!< Return product of L2 and L3 correction factors
    double getJPTL2L3() const { return jptL2L3_; }   //!< Return product of L2 and L3 correction factors for Jet+Track
    double getToL2() const { return l1_*l2_; }         //!< Return factor needed to get L2 corrected from raw jets: L1*L2
    double getToL3() const { return getToL2()*l3_; }   //!< Return factor needed to get L3 corrected from raw jets: L1*L2*L3
    double getToL4() const { return getToL3()*l4_; }   //!< Return factor needed to get L4 corrected from raw jets: L1*L2*L3*L4
    double getToL5() const { return getToL4()*l5_; }   //!< Return factor needed to get L5 corrected from raw jets: L1*L2*L3*L4*L5
    double getToJPTL3() const { return jpt_*l1_*jptL2L3_; }   //!< Return factor needed to get L3 corrected from raw jets for JPT: JPT*L1*JPTL2L3
  private :
    double l1_;      //!< Level 1 correction factor (zero-suppression)
    double l2_;      //!< Level 2 correction factor (relative, in eta)
    double l3_;      //!< Level 3 correction factor (absolute, in pt)
    double l4_;      //!< Level 4 correction factor (electromagnetic fraction)
    double l5_;      //!< Level 5 correction factor (flavor)
    double jpt_;     //!< Jet+Track correction factor
    double jptL2L3_; //!< Product of level 2 and level 3 correction factors for Jet+Track
  };
  
 public:
  Jet(double Et, double EmEt, double HadEt ,double OutEt, double E,
      double eta,double phi, Flavor flavor, double genPt, double dR,
      CorFactors corFactors, const Function& f,
      double (*errfunc)(const double *x, const Measurement *xorig, double err), 
      const Function& gf, double Etmin = 0); 
  virtual ~Jet() {};

  double Et()     const {return Measurement::pt;}                 //!< Return transverse energy Et
  double pt()     const {return Measurement::pt;}                 //!< Return transverse energy Et
  double EmEt()   const {return EMF;}                //!< Return Et from the ECAL part of the towers
  double HadEt()  const {return HadF;}               //!< Return Et from the HCAL part of the towers
  double OutEt()  const {return OutF;}               //!< Return Et from the HOut part of the towers
  double E()      const {return Measurement::E;}    //!< Return energy
  double eta()    const {return Measurement::eta;}  //!< Return pseudorapidity
  double phi()    const {return Measurement::phi;}  //!< Return azimuthal angle
  Flavor flavor() const {return flavor_;}       //!< Return jet flavor
  double genPt()  const {return genPt_;}        //!< Return Pt for corresponding GenJet 
  double ptHat()  const {return ptHat_;}     //!< \f$ \hat{p}_{T} \f$ of the event
  double dR() const {return dR_;}               //!< \f$ \Delta R \f$ between jet and genjet
  const CorFactors& corFactors() const { return corFactors_;}
  //!  \brief Change address of parameters covered by this jet
  //!  \sa Parameters
  // ---------------------------------------------------------
  virtual void ChangeParAddress(double* oldpar, double* newpar) {
    f.changeParBase(oldpar,newpar);
    gf.changeParBase(oldpar,newpar);
  }
  virtual double correctedEt() const { return correctedEt(Et()); }
  virtual double correctedEt(double Et, bool fast = false) const;
  double expectedEt(double truth, double start, double& error,
		    bool fast = false);

  //!  \brief Calculate error from original measurement
  //!
  //!  The error is calculated using the error function
  //!  errf and Measurement::pt, the pt of the original
  //!  measurement.
  //!  
  //!  \return Error of original measurement
  // ---------------------------------------------------------
  virtual double Error() const {return errf(&(Measurement::pt),this,0);}

  //!  \brief Calculate error from given Et
  //!
  //!  The error is calculated using the error function
  //!  errf and a given Et.
  //!
  //!  \param et Jet Et
  //!  \return Error from jet Et
  // ---------------------------------------------------------
  virtual double expectedError(double et) const { return errf(&et,this,0);}

  //!  \brief Get number of parameters of this jet
  //!
  //!  Returns the number of parameters of the
  //!  local and global jet correction functions
  //!  (see Parametrization).
  //!
  //!  \return Number of parameters of this jet
  // ---------------------------------------------------------
  virtual int nPar() const {return f.nPars() + gf.nPars();}

  //!  \brief Represents a corrected jet if one parameter of
  //!         the correction function is varied (for derivative
  //!         calculation)
  //!
  //!  Stores the corrected Et and error if a parameter is
  //!  varied by +/- eps for derivative calculation.
  // ---------------------------------------------------------
  struct ParameterVariation {
    int    parid;        //!< Id of varied parameter
    double upperEt;      //!< Expected Et if parameter is varied by +eps
    double lowerEt;      //!< Expected Et if parameter is varied by -eps
    double upperError;   //!< Expected error if parameter is varied by +eps
    double lowerError;   //!< Expected error if parameter is varied by -eps
    double upperEtDeriv; //!< Derivative of Et if parameter is  varied by +eps
    double lowerEtDeriv; //!< Derivative of Et if parameter is  varied by +eps
    bool operator==(int b) const { return parid == b;} //!< Two ParameterVariation are the same if they have the same parid
  };
  typedef std::vector<ParameterVariation> VariationColl;
  typedef std::vector<ParameterVariation>::const_iterator VariationCollIter;
  virtual const VariationColl& varyPars(double eps, double Et, double start);
  virtual const VariationColl& varyParsDirectly(double eps);

  void print();                       //!< Print some jet members
  static void printInversionStats();  //!< Print some info on inversion

  int parIndex() const { return f.parIndex(); }
 protected:
  mutable VariationColl varcoll;
  virtual double expectedEt(double truth, double start, bool fast = false);

 private: 
  Flavor flavor_;           //!< The jet's Flavor
  double genPt_;            //!< The genjet pt
  double dR_;               //!< \f$ \Delta R \f$ between jet and genjet
  double ptHat_;            //!< \f$ \hat{p}_{T} \f$ of the event
  CorFactors corFactors_;   //!< The correction factors
  double    error;                //!< Stores error for constant error mode
  Function  f;                    //!< Jet correction function
  Function  gf;                   //!< Global jet correction function
  double    (*errf)(const double *x, const Measurement *xorig, double err);   //!< Error function
  double    etmin;                //!< Lower cut on measured Et

  bool      secant(double truth, double& x1, double& x2, double eps);
  bool      falseposition(double truth, double& x1, double& x2, double eps);

  static long long ncalls;        //!< Number of calls of inversion methods secant and falseposition
  static long long ntries;        //!< Number of tries in iteration during inversion
  static long long nfails;        //!< Number of failed tries during inversion
  static long long nwarns;        //!< Number of warnings during inversion

  mutable Measurement temp;
  const double EoverPt;
  class GslImplementation {
    struct rf_par {
      double y;
      const Jet* jet;
      rf_par(double y, const Jet *jet) : y(y), jet(jet) {}
    } par;
    gsl_root_fsolver *s;
    gsl_function F;
    static double rf(double x, void* params) {
      rf_par* p = (rf_par*)params;
      return p->y - p->jet->correctedEt(x,true);
    };
  public:
    GslImplementation(const Jet* jet);
    ~GslImplementation();
    bool root(double truth, double& x1, double& x2, double eps);
  } gsl_impl;
};

#endif
