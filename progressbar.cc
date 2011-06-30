//
// taken from http://teknut.blogspot.com/2009/08/c-console-progress-bar.html
//

#include <iostream>
#include <sstream>

void progressbar(int percent)
{
  static int x = 0;
  char slash[] = {'\\', '|', '/','-'};

  std::string bars;
  for(int i = 5 ; i <= percent ; i+= 5) {
    bars +='#';
  }
  std::cout << '\r'; // carriage return back to beginning of line
  std::cout << bars << " " << slash[x] << " " << percent << " %" << std::flush; // print the bars and percentage
  x++; // increment to make the slash appear to rotate
  if(x == 4) x = 0; // reset slash animation
} 
