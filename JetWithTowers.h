#ifndef JETWITHTOWERS_H
#define JETWITHTOWERS_H

#include"Jet.h"

#include <vector>
#include <map>


//!
//!    \brief Class for jets with towers 
//!
//!    \author Hartmut Stadie
//!    \date 2008/12/25
//!    $Id: JetWithTowers.h,v 1.22 2010/05/19 13:34:49 stadie Exp $
// ----------------------------------------------------------------   
class JetWithTowers : public Jet
{
 public:
  JetWithTowers(double Et, double EmEt, double HadEt ,double OutEt, double E,
		double eta,double phi, double phiphi, double etaeta, 
		Flavor flavor, double genPt, double dR, CorFactors* corFactors,
		const Function& f,
		double (*errfunc)(const double *x, const Measurement *xorig, double err), 
		const Function& gf, double Etmin = 0); 
  virtual ~JetWithTowers(); 
  virtual int nPar() const {return Jet::nPar() + towerpars_.size() * ntowerpars_;}
  virtual void changeParAddress(double* oldpar, double* newpar);
  virtual double correctedEt(double Et,bool fast = false) const; 
  virtual double error() const;
  virtual double expectedError(double et) const;
  // varies all parameters for this jet by eps and returns a vector of the
  // parameter id and the Et for the par + eps and par - eps variation
  virtual const VariationColl& varyPars(const double* eps, double Et, double start);
  virtual const VariationColl& varyParsDirectly(const double* eps, bool computeDeriv = false);

  void addTower(double Et, double EmEt, double HadEt ,double OutEt, double E,
		double eta,double phi,const Function& f,
		double (*errfunc)(const double *x, const Measurement *xorig, double err));
  virtual Jet* clone() const { return new JetWithTowers(*this);} //!< Clone this jet
 private:
  JetWithTowers(const JetWithTowers& j); //!< disallow copies!
  class Tower : public Measurement {
  public:
    Tower(double Et, double EmEt, double HadEt ,double OutEt, double E,
	  double eta,double phi, double alpha, const Function& func,
	  double (*errfunc)(const double *x, const Measurement *xorig, double err));
    Tower(double Et, double EmEt, double HadEt ,double OutEt, double E,
	  double EmEttrue, double HadEttrue, double OutEttrue,
	  double eta,double phi, double alpha, const Function& func,
	  double (*errfunc)(const double *x, const Measurement *xorig, double err));
    virtual ~Tower() {}
    double Et()     const {return pt;}
    double EmEt()   const {return EMF;}
    double HadEt()  const {return HadF;}
    double OutEt()  const {return OutF;}
    double E()      const {return Measurement::E;}
    double eta()    const {return Measurement::eta;}
    double phi()    const {return Measurement::phi;}
    double projectionToJetAxis() const {return alpha_;}
    double fractionOfJetHadEt() const { return fraction_;}
    void setFractionOfJetHadEt(double frac) const { fraction_ = frac;}
    void changeParAddress(double* oldpar, double* newpar) {f_.changeParBase(oldpar,newpar);}
    double correctedHadEt(double HadEt) const;
    double lastCorrectedHadEt() const { return lastCorHadEt_;}  
    double error() const {return errf_(&(Measurement::pt),this,error_);}
    double expectedError(double et) const { return  errf_(&et,this,error_);}
    int nPar() const {return f_.nPars();}
    int firstPar() const {return f_.parIndex();}
    double *par() const {return f_.firstPar();}
  private:
    double alpha_;                //!< Projection factor onto jet axis
    double error_;                //!< Error for constant error mode
    double mEttrue_;              //!< True total transverse energy
    double mEmEttrue_;            //!< True Et from the ECAL part of the tower		
    double mHadEttrue_;           //!< True Et from the HCAL part of the towers
    double mOutEttrue_;           //!< True Et from the HO part of the tower
    mutable Measurement temp_;
    mutable double lastCorHadEt_;
    mutable double fraction_;
    Function f_;
    double (*errf_)(const double *x, const Measurement *xorig, double err);

    friend class JetWithTowers;
  };
  typedef std::vector<Tower*> TowerColl;
  typedef TowerColl::iterator TowerCollIter;
  typedef TowerColl::const_iterator TowerCollConstIter;
  TowerColl towers_;
  int ntowerpars_;
  std::map<int,double*> towerpars_;
};

#endif
