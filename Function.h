#ifndef FUNCTION_H
#define FUNCTION_H
 
class TMeasurement;

class Function {
 public:
  Function(double (*func)(const TMeasurement*, const double*),
	   double (*invfunc)(const TMeasurement*, const double*), 
	   double *firstpar, int parindex, int npars) 
    : func(func),invfunc(invfunc),firstpar(firstpar),parindex(parindex),npars(npars)
    {}
  //double (*f)(const TMeasurement*, const double*) const { return func;}
  double* firstPar() const { return firstpar;}
  int parIndex() const { return parindex;}
  int nPars() const { return npars;}
  double operator()(const TMeasurement* x) const { return func(x,firstpar);}
  void changeParBase(double* oldpar, double* newpar) { firstpar += newpar - oldpar;}

  bool hasInverse() { return invfunc;}
  double inverse(const TMeasurement* x) const { return invfunc ? invfunc(x,firstpar) : 0;}
 private:
    double (*func)(const TMeasurement*, const double*);
    double (*invfunc)(const TMeasurement*, const double*);
    double *firstpar;
    int parindex;
    int npars;
};
#endif
