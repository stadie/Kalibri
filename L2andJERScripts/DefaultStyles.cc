#include "DefaultStyles.h"

DefaultStyles::DefaultStyles() {


  init();
}

void DefaultStyles::init() {

   markers_.push_back(20);
   markers_.push_back(24);
   markers_.push_back(21);
   markers_.push_back(25);
   markers_.push_back(22);
   markers_.push_back(26);
   markers_.push_back(29);
   markers_.push_back(30);
   markers_.push_back(23);
   markers_.push_back(28);
   markers_.push_back(34);

  colors_.push_back(1);
  colors_.push_back(2);
  colors_.push_back(4);
  colors_.push_back(3);
  colors_.push_back(46);
  colors_.push_back(9);
  colors_.push_back(12);

}


void DefaultStyles::setStyle (TString stylelabel){
  if(stylelabel=="PFComp"){
    colors_.clear();
    markers_.clear();

   markers_.push_back(20);
   markers_.push_back(24);
   markers_.push_back(21);
   markers_.push_back(25);
   markers_.push_back(22);
   markers_.push_back(26);
   markers_.push_back(29);
   markers_.push_back(30);
   markers_.push_back(23);
   markers_.push_back(28);
   markers_.push_back(34);

  colors_.push_back(2);
  colors_.push_back(4);
  colors_.push_back(3);
  colors_.push_back(46);
  colors_.push_back(9);
  colors_.push_back(12);
  colors_.push_back(1);

  }


  if(stylelabel=="Confidence"){
    colors_.clear();
    markers_.clear();

   markers_.push_back(1);
   markers_.push_back(1);
   markers_.push_back(1);
   markers_.push_back(1);
   markers_.push_back(1);
   markers_.push_back(1);
   markers_.push_back(1);
   markers_.push_back(1);
   markers_.push_back(1);
   markers_.push_back(1);
   markers_.push_back(1);

  colors_.push_back(kGray);
  colors_.push_back(kOrange-2);
  colors_.push_back(3);
  colors_.push_back(46);
  colors_.push_back(9);
  colors_.push_back(12);
  colors_.push_back(1);

  }

}
