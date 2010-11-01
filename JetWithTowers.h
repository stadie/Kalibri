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
//!    $Id: JetWithTowers.h,v 1.24 2010/10/20 11:28:19 stadie Exp $
// ----------------------------------------------------------------   
class JetWithTowers : public Jet
{
 public:
  JetWithTowers(float Et, float EmEt, float HadEt ,float OutEt, float E,
		float eta,float phi, float phiphi, float etaeta, 
		Flavor flavor, float genPt, float dR, CorFactors* corFactors,
		const Function& f,
		float (*errfunc)(const float *x, const Measurement *xorig, float err), 
		const Function& gf); 
  virtual ~JetWithTowers(); 
  virtual int nPar() const {return Jet::nPar() + towerpars_.size() * ntowerpars_;}
  virtual void setParameters(Parameters* param);
  virtual float correctedEt(float Et,bool fast = false) const; 
  virtual float error() const;
  virtual float expectedError(float et) const;
  // varies all parameters for this jet by eps and returns a vector of the
  // parameter id and the Et for the par + eps and par - eps variation
  virtual const Parameters::VariationColl& varyPars(const double* eps, float Et, float start);
  virtual const Parameters::VariationColl& varyParsDirectly(const double* eps, bool computeDeriv = false);

  void addTower(float Et, float EmEt, float HadEt ,float OutEt, float E,
		float eta,float phi, const Function& f,
		float (*errfunc)(const float *x, const Measurement *xorig, float err));
  virtual Jet* clone() const { return new JetWithTowers(*this);} //!< Clone this jet
 private:
  JetWithTowers(const JetWithTowers& j); //!< disallow copies!
  class Tower : public Measurement {
  public:
    Tower(float Et, float EmEt, float HadEt ,float OutEt, float E,
	  float eta,float phi, float alpha, const Function& func,
	  float (*errfunc)(const float *x, const Measurement *xorig, float err));
    Tower(float Et, float EmEt, float HadEt ,float OutEt, float E,
	  float EmEttrue, float HadEttrue, float OutEttrue,
	  float eta,float phi, float alpha, const Function& func,
	  float (*errfunc)(const float *x, const Measurement *xorig, float err));
    virtual ~Tower() {}
    float Et()     const {return pt;}
    float EmEt()   const {return EMF;}
    float HadEt()  const {return HadF;}
    float OutEt()  const {return OutF;}
    float E()      const {return Measurement::E;}
    float eta()    const {return Measurement::eta;}
    float phi()    const {return Measurement::phi;}
    float projectionToJetAxis() const {return alpha_;}
    float fractionOfJetHadEt() const { return fraction_;}
    void setFractionOfJetHadEt(float frac) const { fraction_ = frac;}
    const Function& setParameters(Parameters* param);
    float correctedHadEt(float HadEt) const;
    float lastCorrectedHadEt() const { return lastCorHadEt_;}  
    float error() const {return errf_(&(Measurement::pt),this,error_);}
    float expectedError(float et) const { return  errf_(&et,this,error_);}
    int nPar() const {return f_->nPars();}
    int firstPar() const {return f_->parIndex();}
    double *par() const {return f_->firstPar();}
  private:
    float alpha_;                //!< Projection factor onto jet axis
    float error_;                //!< Error for constant error mode
    float mEttrue_;              //!< True total transverse energy
    float mEmEttrue_;            //!< True Et from the ECAL part of the tower		
    float mHadEttrue_;           //!< True Et from the HCAL part of the towers
    float mOutEttrue_;           //!< True Et from the HO part of the tower
    mutable Measurement temp_;
    mutable float lastCorHadEt_;
    mutable float fraction_;
    const Function* f_;
    float (*errf_)(const float *x, const Measurement *xorig, float err);

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
