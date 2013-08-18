#if defined(_WIN32)
#include <time.h>
#else
#include <sys/time.h>
#endif

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <cassert>

#include <vector>
#include <iterator>
#include <numeric>
#include <algorithm>

#include <mpi.h>

// Random number in (0,1)
inline double get_random() {
  return drand48();
}
// Random number in (A,B)
inline double get_random(double A, double B) {
  return A + (B-A)*get_random();
}

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

/** Read a line from @a s, parse it as type T, and store it in @a value.
 * @param[in]   s      input stream
 * @param[out]  value  value returned if the line in @a s doesn't parse
 *
 * If the line doesn't parse correctly, then @a s is set to the "failed"
 * state. Ignores blank lines and lines that start with '#'.
 */
template <typename T>
std::istream& getline_parsed(std::istream& s, T& value) {
  std::string str;
  do {
    getline(s, str);
  } while (s && (str.empty() || str[0] == '#'));
  std::istringstream is(str);
  is >> value;
  if (is.fail())
    s.setstate(std::istream::failbit);
  return s;
}

template <typename T>
std::istream& operator>>(std::istream& s, std::vector<T>& v) {
  T temp;
  while (getline_parsed(s, temp))
    v.push_back(temp);
  return s;
}

template <typename T>
std::ostream& operator<<(std::ostream& s, const std::vector<T>& v) {
  std::copy(v.begin(), v.end(), std::ostream_iterator<T>(s,"\n"));
  return s;
}

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


// Problem specific -- XXX: NOT NEEDED
#define NUMPOINTS 100
#define PHIDATA "data/phiDataSmall.txt"
#define SIGMADATA "data/sigmaDataSmall.txt"
#define MASTER 0
#define DATADIM 3

// tool to help print out phi
void printPhi (double phi[NUMPOINTS]) {
  std::ofstream printFile;
  printFile.open("phi.txt");
  for (int i = 0; i < NUMPOINTS; ++i) {
    printFile << phi[i] << "\n";
  }
  printFile.close();
}
