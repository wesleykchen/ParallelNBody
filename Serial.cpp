#include "Util.hpp"

// Serial version of n-body algorithm, with 3 points

int main(int argc, char *argv[])
{
  string x, y, z;
  ifstream inputFile(PHIDATA);
  double data[NUMPOINTS][DATADIM];

  if (inputFile.is_open()) {
    int index = 0;
    while ( true ) {
      getline( inputFile, x, '\t');
      if (inputFile.eof()) {
	break;
      }
      getline( inputFile, y, '\t');
      getline( inputFile, z, '\n');

      //cout << x << "." << y << "." << z << endl;

      data[index][0] = atof(x.c_str());
      data[index][1] = atof(y.c_str());
      data[index][2] = atof(z.c_str());

      ++index;
    }
    inputFile.close();
    inputFile.clear();
  } else {
    cout << "Unable to open phi file";
  }

  double sigma[NUMPOINTS];
  inputFile.open(SIGMADATA);
  if (inputFile.is_open()) {
    int index = 0;
    while ( true ) {
      getline( inputFile, x, '\n');
      if (inputFile.eof()) {
	break;
      }
      sigma[index] = atof(x.c_str());

      ++index;
    }
    inputFile.close();
  } else {
    cout << "Unable to open sigma file";
  }

  // cout << sigma[1] << endl;

  Clock timer;

  double phi[NUMPOINTS];
  // clear initial values

  for (int i = 0; i < NUMPOINTS; ++i) {
    phi[i] = 0;
  }

  timer.start();
  for (int i = 0; i < NUMPOINTS; ++i) {
    for (int j = 0; j < NUMPOINTS; ++j) {
      phi[i] += evaluate(data[i][0], data[i][1], data[i][2],
			 data[j][0], data[j][1], data[j][2]) * sigma[j];
    }
  }
  double time = timer.elapsed();

  double checkSum = 0;
  for (int i = 0; i < NUMPOINTS; ++i) {
  	checkSum += phi[i];
  }

  //printf("The last xI computed answer is: %e\n", phi[99]);
  printPhi(phi);
  printf("Serial - checksum answer is: %e\n", checkSum);
  printf("Computed in %e seconds\n", time);
  return 0;
}
