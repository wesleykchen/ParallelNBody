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

/** Integer divide, rounded up
 * @param[in] a Numerator
 * @param[in] b Denominator
 * @returns If b divides a, then a/b
 *                          else a/b + 1
 */
inline unsigned idiv_up(unsigned a, unsigned b) {
  return (a+b-1)/b;
}

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

template <class T>
inline std::string to_string(const T& t) {
  std::stringstream ss;
  ss << t;
  return ss.str();
}

template <class T>
inline T string_to_(const std::string& s) {
  T val;
  std::istringstream(s) >> val;
  return val;
}

// 1/R kernel
template <typename P>
double evaluate(const P& x, const P& y) {
  double dx = x[0] - y[0];
  double dy = x[1] - y[1];
  double dz = x[2] - y[2];
  double R = sqrt(dx*dx + dy*dy + dz*dz);
  double invR = 1.0 / R;
  if (R < 1e-10) invR = 0;
  return invR;
}

#if 0
// 1/R^2 kernel
template <typename P>
double evaluate(const P& x, const P& y) {
  double dx = x[0] - y[0];
  double dy = x[1] - y[1];
  double dz = x[2] - y[2];
  double R2 = (dx*dx + dy*dy + dz*dz);
  double invR2 = 1.0 / R2;
  if (R2 < 1e-20) invR2 = 0;
  return invR2;
}
#endif


/** Asymmetric block P2P using the evaluation operator
 * r_i += sum_j K(t_i, s_j) * c_j
 *
 * @param[in] s_first,s_last  Iterator range to sources (s_j)
 * @param[in] c_first         Associated iterator to charges (c_j)
 * @param[in] t_first,t_last  Iterator range to targets (t_i)
 * @param[in,out] r_first     Associated iterator to results (r_i)
 */
template <typename SourceIter, typename ChargeIter,
          typename TargetIter, typename ResultIter>
inline void
block_eval(SourceIter s_first, SourceIter s_last, ChargeIter c_first,
           TargetIter t_first, TargetIter t_last, ResultIter r_first)
{
  typedef typename std::iterator_traits<SourceIter>::value_type source_type;
  typedef typename std::iterator_traits<ChargeIter>::value_type charge_type;
  typedef typename std::iterator_traits<TargetIter>::value_type target_type;
  typedef typename std::iterator_traits<ResultIter>::value_type result_type;

  for ( ; t_first != t_last; ++t_first, ++r_first) {
    const target_type& t = *t_first;
    result_type& r       = *r_first;

    SourceIter s = s_first;
    ChargeIter c = c_first;
    for ( ; s != s_last; ++s, ++c)
      r += evaluate(t,*s) * (*c);
  }
}


// Problem specific -- XXX: NOT NEEDED
#define NUMPOINTS 100
#define PHIDATA "data/phiDataSmall.txt"
#define SIGMADATA "data/sigmaDataSmall.txt"
#define MASTER 0
#define DATADIM 3

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

// tool to help print out phi
void printPhi (double phi[NUMPOINTS]) {
  std::ofstream printFile;
  printFile.open("phi.txt");
  for (int i = 0; i < NUMPOINTS; ++i) {
    printFile << phi[i] << "\n";
  }
  printFile.close();
}
