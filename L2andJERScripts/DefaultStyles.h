#ifndef DefaultStyles_h
#define DefaultStyles_h

#include "TString.h"

class DefaultStyles {
 public:

  DefaultStyles();
  int getColor(int color_i) {
    if(color_i>colors_.size())std::cout << "color index out of range" <<std::endl;
    assert(color_i<colors_.size());
    return colors_.at(color_i);
  }
  int getMarker(int marker_i) {
    if(marker_i>markers_.size())std::cout << "marker index out of range" <<std::endl;
    assert(marker_i<markers_.size());
    return markers_.at(marker_i);
}
  void setStyle(TString stylelabel = "PFComp");
  
 private:
  void init();
  std::vector <int> colors_;
  std::vector <int> markers_;
};



#endif 


