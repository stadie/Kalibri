#ifndef DefaultStyles_h
#define DefaultStyles_h

#include "TString.h"

class DefaultStyles {
 public:

  DefaultStyles();
  int getColor(int color_i) {return colors_.at(color_i);}
  int getMarker(int marker_i) {return markers_.at(marker_i);}
  void setStyle(TString stylelabel = "PFComp");
  
 private:
  void init();
  std::vector <int> colors_;
  std::vector <int> markers_;
};



#endif 


