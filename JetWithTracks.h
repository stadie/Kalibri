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
//!    $Id: JetWithTracks.h,v 1.11 2010/05/19 13:34:49 stadie Exp $
// ---------------------------------------------------------------   
class JetWithTracks : public Jet
{
 public:
  JetWithTracks(double Et, double EmEt, double HadEt ,double OutEt, double E,
		double eta,double phi, double phiphi, double etaeta, 
		Flavor flavor,double genPt, double dR, CorFactors* corFactors,
		const Function& f,
		double (*errfunc)(const double *x, const Measurement *xorig, double err), 
		const Function& gf, double Etmin = 0); 
  virtual ~JetWithTracks(); 
  virtual int nPar() const {return Jet::nPar() + trackpars_.size() * ntrackpars_;}
  virtual void changeParAddress(double* oldpar, double* newpar);
  virtual double correctedEt(double Et,bool fast = false) const; 
  virtual double error() const;
  virtual double expectedError(double et) const;
  // varies all parameters for this jet by eps and returns a vector of the
  // parameter id and the Et for the par + eps and par - eps variation
  virtual const VariationColl& varyPars(const double* eps, double Et, double start);
  virtual const VariationColl& varyParsDirectly(const double* eps, bool computeDeriv);

  void addTrack(double Et, double EmEt, double HadEt ,double OutEt, double E,
		double eta,double phi,int TrackId, int TowerId, double DR, double DRout,
		double etaOut, double phiOut, double EM1, double EM5, double Had1, 
		double Had5, double TrackChi2, int NValidHits, bool TrackQualityT, 
		double MuDR, double MuDE, double Efficiency, const Function& f,
		double (*errfunc)(const double *x, const Measurement *xorig, double err));
  virtual Jet* clone() const { return new JetWithTracks(*this);} //!< Clone this jet
 protected:
  virtual double expectedEt(double truth, double start, bool fast = false);
 private:
  JetWithTracks(const JetWithTracks& j); //!< disallow copies!
  class Track : public TTrack {
  public:
    Track(double Et, double EmEt, double HadEt ,double OutEt, double E,
	  double eta,double phi, int TrackId, int TowerId, double DR, double DRout,
	  double etaOut, double phiOut, double EM1, double EM5, double Had1, double Had5,
	  double TrackChi2, int NValidHits, bool TrackQualityT, double MuDR, double MuDE,
	  double Efficiency, const Function& func,
	  double (*errfunc)(const double *x, const Measurement *xorig, double err));
    virtual ~Track() {}
    double Et()     const {return Measurement::pt;}
    double EmEt()   const {return Measurement::EMF;}
    double HadEt()  const {return Measurement::HadF;}
    double OutEt()  const {return Measurement::OutF;}
    double E()      const {return Measurement::E;}
    double eta()    const {return Measurement::eta;}
    double phi()    const {return Measurement::phi;}
    int trackId() const { return TTrack::TrackId;}
    int towerId() const { return TTrack::TowerId;}
    double dR() const { return TTrack::DR;}
    double dRout() const { return TTrack::DRout;}
    bool goodTrack() const { return TTrack::TrackQualityT;}
    void changeParAddress(double* oldpar, double* newpar) {f_.changeParBase(oldpar,newpar);}
    double expectedEt() const;
    double error() const {return 0;}
    int nPar() const {return f_.nPars();}
    int firstPar() const {return f_.parIndex();}
    double *par() const {return f_.firstPar();}
  private:
    Function f_;
    double (*errf_)(const double *x, const Measurement *xorig, double err);

    friend class JetWithTracks;
  };
  typedef std::vector<Track*> TrackColl;
  typedef TrackColl::iterator TrackCollIter;
  typedef TrackColl::const_iterator TrackCollConstIter;
  TrackColl tracks_;
  int ntrackpars_;
  std::map<int,double*> trackpars_;
  mutable double expectedCaloEt_;
  mutable double trackPt_;
  
};

#endif
