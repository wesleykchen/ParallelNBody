#if defined(_WIN32)
#include <time.h>
#else
#include <sys/time.h>
#endif

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cassert>

#include <mpi.h>

using namespace std;

/** Clock class, useful when timing code.
 */
struct Clock {
  /** Construct a Clock and start timing. */
  Clock() {
    start();
  }
  /** Start the clock. */
  inline void start() {
    time_ = now();
  }
  /** Return the seconds elapsed since the last start. */
  inline double elapsed() const {
    return sec_diff(time_, now());
  }
  /** Return the seconds difference between two Clocks */
  inline double operator-(const Clock& clock) const {
    return sec_diff(time_, clock.time_);
  }
 private:
  timeval time_;
  inline static timeval now() {
    timeval tv;
    gettimeofday(&tv, 0);
    return tv;
  }
  // Return the time difference (t2 - t1) in seconds
  inline static double sec_diff(const timeval& t1, const timeval& t2) {
    timeval dt;
    timersub(&t2, &t1, &dt);
    return dt.tv_sec + 1e-6 * dt.tv_usec;
  }
};

/*
// TODO: maybe use?
struct Point {
  double x, y, z;

  Point(double _x, double _y, double _z)
      : x(_x), y(_y), z(_z) {
  }

  template <typename Stream>
  friend Stream& operator<<(Stream& s, const Point& p) {
    return s << p.x << '\t' << p.y << '\t' << p.z << std::endl;
  }
  template <typename Stream>
  friend Stream& operator>>(Stream& s, Poi) {
    return s >> p.x  >> p.y >> p.z;
  }
};


*/

// Random number in (0,1)
inline double get_random() {
  return drand48();
}
// Random number in (A,B)
inline double get_random(double A, double B) {
  return A + (B-A)*get_random();
}

// Problem specific
#define NUMPOINTS 100000
#define PHIDATA "phiData100k.txt"
#define SIGMADATA "sigmaData100k.txt"
#define MASTER 0
#define DATADIM 3


// tool to help print out phi
void printPhi (double phi[NUMPOINTS]) {
    ofstream printFile;
    printFile.open("phi.txt");
    for (int i = 0; i < NUMPOINTS; ++i) {
        printFile << phi[i] << "\n";
    }
    printFile.close();
}

/*
double* readData(double phiData[][3]) {
  string x, y, z;
  ifstream phiFile("phi.txt");
  if (phiFile.is_open()) {
    int index = 0;
    while ( true ) {

      getline( phiFile, x, '\t');
      if (phiFile.eof()) {
	break;
      }
      getline( phiFile, y, '\t');
      getline( phiFile, z, '\n');

      //cout << x << "." << y << "." << z << endl;

      phiData[index][0] = atof(x.c_str());
      phiData[index][1] = atof(y.c_str());
      phiData[index][2] = atof(z.c_str());

      ++index;
    }
    phiFile.close();
  }
  else {
    printf("Unable to open file");
  }
  double** dataPtr = phiData;
  return dataPtr;
}
*/

// 1/R kernel
double evaluate(double x1, double x2, double x3,
		double y1, double y2, double y3) {
  double dx = x1 - y1;
  double dy = x2 - y2;
  double dz = x3 - y3;
  double R = sqrt(dx*dx + dy*dy + dz*dz);
  double invR = 1.0 / R;
  if (R < 1e-10) invR = 0;
  return invR;
}

#if 0
// 1/R^2 kernel
double evaluate(double x1, double x2, double x3,
		double y1, double y2, double y3) {
  double dx = x1 - y1;
  double dy = x2 - y2;
  double dz = x3 - y3;
  double R2 = (dx*dx + dy*dy + dz*dz);
  double invR2 = 1.0 / R2;
  if (R2 < 1e-20) invR2 = 0;
  return invR2;
}
#endif
