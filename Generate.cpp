#include "Util.hpp"

int main(int argc, char *argv[])
{
  srand( time(NULL) );
  
  ofstream outputFile;
  outputFile.open(PHIDATA);
  
  for (int i = 0; i < NUMPOINTS; ++i) {
    // get random number x,y,z between 1 and 100
    double x = get_random();
    double y = get_random();
    double z = get_random();
    
    outputFile << x << "\t" << y << "\t" << z << "\n";
  }
  
  outputFile.close();
  outputFile.clear();

  outputFile.open(SIGMADATA);

 for (int i = 0; i < NUMPOINTS; ++i) {
    // get random number x,y,z between 1 and 100
    double sigma = get_random();
    
    outputFile << sigma << "\n";
  }
  
  outputFile.close();
  return 0;
}
