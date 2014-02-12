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

#include "P2P.hpp"

/** Integer divide, rounded up
 * @param[in] a Numerator
 * @param[in] b Denominator
 * @returns If b divides a, then a/b, else a/b + 1
 */
inline unsigned idiv_up(unsigned a, unsigned b) {
  return (a+b-1)/b;
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

// Problem specific -- XXX: NOT NEEDED
#define NUMPOINTS 1024
#define SOURCE_DATA "data/sourceDataSmall.txt"
#define CHARGE_DATA "data/chargeDataSmall.txt"
#define MASTER 0
