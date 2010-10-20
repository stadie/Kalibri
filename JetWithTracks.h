#ifndef JETWITHTRACKS_H
#define JETWITHTRACKS_H

#include"Jet.h"

#include <vector>
#include <map>


//!
//!    \brief Class for jets with towers 
//!
//!    \author Hartmut Stadie
//!    \date 2008/12/25
//!    $Id: JetWithTracks.h,v 1.12 2010/05/19 16:01:42 stadie Exp $
// ---------------------------------------------------------------   
class JetWithTracks : public Jet
{
 public:
  JetWithTracks(float Et, float EmEt, float HadEt ,float OutEt, float E,
		float eta,float phi, float phiphi, float etaeta, 
		Flavor flavor,float genPt, float dR, CorFactors* corFactors,
		const Function& f,
		float (*errfunc)(const float *x, const Measurement *xorig, float err), 
		const Function& gf); 
  virtual ~JetWithTracks(); 
  virtual int nPar() const {return Jet::nPar() + trackpars_.size() * ntrackpars_;}
  virtual void changeParAddress(double* oldpar, double* newpar);
  virtual float correctedEt(float Et,bool fast = false) const; 
  virtual float error() const;
  virtual float expectedError(float et) const;
  // varies all parameters for this jet by eps and returns a vector of the
  // parameter id and the Et for the par + eps and par - eps variation
  virtual const VariationColl& varyPars(const double* eps, float Et, float start);
  virtual const VariationColl& varyParsDirectly(const double* eps, bool computeDeriv);

  void addTrack(float Et, float EmEt, float HadEt ,float OutEt, float E,
		float eta,float phi,int TrackId, int TowerId, float DR, float DRout,
		float etaOut, float phiOut, float EM1, float EM5, float Had1, 
		float Had5, float TrackChi2, int NValidHits, bool TrackQualityT, 
		float MuDR, float MuDE, float Efficiency, const Function& f,
		float (*errfunc)(const float *x, const Measurement *xorig, float err));
  virtual Jet* clone() const { return new JetWithTracks(*this);} //!< Clone this jet
 protected:
  virtual float expectedEt(float truth, float start, bool fast = false);
 private:
  JetWithTracks(const JetWithTracks& j); //!< disallow copies!
  class Track : public TTrack {
  public:
    Track(float Et, float EmEt, float HadEt ,float OutEt, float E,
	  float eta,float phi, int TrackId, int TowerId, float DR, float DRout,
	  float etaOut, float phiOut, float EM1, float EM5, float Had1, float Had5,
	  float TrackChi2, int NValidHits, bool TrackQualityT, float MuDR, float MuDE,
	  float Efficiency, const Function& func,
	  float (*errfunc)(const float *x, const Measurement *xorig, float err));
    virtual ~Track() {}
    float Et()     const {return Measurement::pt;}
    float EmEt()   const {return Measurement::EMF;}
    float HadEt()  const {return Measurement::HadF;}
    float OutEt()  const {return Measurement::OutF;}
    float E()      const {return Measurement::E;}
    float eta()    const {return Measurement::eta;}
    float phi()    const {return Measurement::phi;}
    int trackId() const { return TTrack::TrackId;}
    int towerId() const { return TTrack::TowerId;}
    float dR() const { return TTrack::DR;}
    float dRout() const { return TTrack::DRout;}
    bool goodTrack() const { return TTrack::TrackQualityT;}
    void changeParAddress(double* oldpar, double* newpar) {f_.changeParBase(oldpar,newpar);}
    float expectedEt() const;
    float error() const {return 0;}
    int nPar() const {return f_.nPars();}
    int firstPar() const {return f_.parIndex();}
    double *par() const {return f_.firstPar();}
  private:
    Function f_;
    float (*errf_)(const float *x, const Measurement *xorig, float err);

    friend class JetWithTracks;
  };
  typedef std::vector<Track*> TrackColl;
  typedef TrackColl::iterator TrackCollIter;
  typedef TrackColl::const_iterator TrackCollConstIter;
  TrackColl tracks_;
  int ntrackpars_;
  std::map<int,double*> trackpars_;
  mutable float expectedCaloEt_;
  mutable float trackPt_;
  
};

#endif
