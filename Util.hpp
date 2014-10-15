#pragma once

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
#include <chrono>

#include <mpi.h>

#include "P2P.hpp"
#include "numeric/Norm.hpp"

/** Integer divide, rounded up
 * @param[in] a Numerator
 * @param[in] b Denominator
 * @returns If b divides a, then a/b, else a/b + 1
 */
inline constexpr unsigned idiv_up(unsigned a, unsigned b) {
  return (a+b-1)/b;
}

/** Positive modulus
 * @param[in] i Dividend
 * @param[in] n Divisor
 * @returns r such that 0 <= r < n and i == a*n+r for some integer a.
 */
inline constexpr int pos_mod(int i, unsigned n) {
  return ((i % n) + n) % n;
}

/** Clock class, useful when timing code.
 */
class Clock {
 public:
  typedef std::chrono::high_resolution_clock clock;
  typedef typename clock::time_point         time_point;
  typedef typename clock::duration           duration_type;
  typedef typename duration_type::rep        tick_type;

  // Default constructor
  Clock() : starttime_(clock::now()) {}
  // Restart the timer on this Clock
  void start() {
    starttime_ = clock::now();
  }
  // Get the duration on this Clock
  duration_type duration() const {
    return clock::now() - starttime_;
  }
  // Get the seconds on this Clock
  double elapsed() const {
    typedef std::chrono::duration<double> units;
    return std::chrono::duration_cast<units>(duration()).count();
  }
 private:
  time_point starttime_;
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
  std::streamsize p = s.precision();
  s.precision(20);
  std::copy(v.begin(), v.end(), std::ostream_iterator<T>(s,"\n"));
  s.precision(p);
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
  assert(exact.size() == result.size());

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

    //if (rel_error > 1e-10) std::cout << k << ": " << result[k] << "\t" << exact[k] << std::endl;
  }
  double tot_rel_err = sqrt(tot_error_sq/tot_norm_sq);
  std::cout << "Vector  relative error: " << tot_rel_err << std::endl;

  double ave_rel_err = tot_ind_rel_err / result.size();
  std::cout << "Average relative error: " << ave_rel_err << std::endl;

  std::cout << "Maximum relative error: " << max_ind_rel_err << std::endl;
}

#define MASTER 0
