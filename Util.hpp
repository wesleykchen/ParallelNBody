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
#include <tuple>
#include <iterator>
#include <numeric>
#include <algorithm>

#include <mpi.h>

#include "P2P.hpp"
#include "numeric/Norm.hpp"

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

template<typename Type, unsigned N, unsigned Last>
struct tuple_printer {
  static void print(std::ostream& out, const Type& value) {
    out << std::get<N>(value) << ", ";
    tuple_printer<Type, N + 1, Last>::print(out, value);
  }
};
template<typename Type, unsigned N>
struct tuple_printer<Type, N, N> {
  static void print(std::ostream& out, const Type& value) {
    out << std::get<N>(value);
  }
};
template<typename... Types>
std::ostream& operator<<(std::ostream& out, const std::tuple<Types...>& value) {
  out << "(";
  tuple_printer<std::tuple<Types...>, 0, sizeof...(Types) - 1>::print(out, value);
  out << ")";
  return out;
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

template <typename result_type>
void print_error(const std::vector<result_type>& exact,
                 const std::vector<result_type>& result) {
  double tot_error_sq = 0;
  double tot_norm_sq = 0;
  double tot_ind_rel_err = 0;
  double max_ind_rel_err = 0;
  for (unsigned k = 0; k < result.size(); ++k) {
    // Individual relative error
    double rel_error = norm(exact[k] - result[k]) / norm(exact[k]);
    tot_ind_rel_err += rel_error;
    // Maximum relative error
    max_ind_rel_err  = std::max(max_ind_rel_err, rel_error);

    // Total relative error
    tot_error_sq += normSq(exact[k] - result[k]);
    tot_norm_sq  += normSq(exact[k]);
  }
  double tot_rel_err = sqrt(tot_error_sq/tot_norm_sq);
  std::cout << "Vector  relative error: " << tot_rel_err << std::endl;

  double ave_rel_err = tot_ind_rel_err / result.size();
  std::cout << "Average relative error: " << ave_rel_err << std::endl;

  std::cout << "Maximum relative error: " << max_ind_rel_err << std::endl;
}

// Problem specific -- XXX: NOT NEEDED
#define NUMPOINTS 1024
#define SOURCE_DATA "data/sourceDataSmall.txt"
#define CHARGE_DATA "data/chargeDataSmall.txt"
#define MASTER 0
