#pragma once

#include <random>
#include <limits>

namespace meta {

static std::mt19937 default_generator;

template <typename T>
struct random;

template <>
struct random<double> {
  static double get(double a, double b) {
    std::uniform_real_distribution<double> dist(a, b);
    return dist(default_generator);
  }
  static double get() {
    return get(0,1);
  }
};

template <>
struct random<float> {
  static float get(float a, float b) {
    std::uniform_real_distribution<float> dist(a, b);
    return dist(default_generator);
  }
  static float get() {
    return get(0,1);
  }
};

template <>
struct random<unsigned> {
  static unsigned get(unsigned a, unsigned b) {
    std::uniform_int_distribution<unsigned> dist(a, b);
    return dist(default_generator);
  }
  static unsigned get() {
    return get(0, std::numeric_limits<unsigned>::max());
  }
};

template <>
struct random<int> {
  static int get(int a, int b) {
    std::uniform_int_distribution<int> dist(a, b);
    return dist(default_generator);
  }
  static int get() {
    return get(0, std::numeric_limits<int>::max());
  }
};

} // end namepsace fmmtl


#include <complex>
using std::complex;

namespace meta {
template <typename T>
struct random<complex<T> > {
  static complex<T> get(T a, T b) {
    return complex<T>(random<T>::get(a,b), random<T>::get(a,b));
  }
  static complex<T> get() {
    return get(T(0), T(1));
  }
};
} // end namespace fmmtl


#include "numeric/Vec.hpp"

namespace meta {
template <std::size_t N, typename T>
struct random<Vec<N,T> > {
  static Vec<N,T> get(T a, T b) {
    Vec<N,T> v;
    for (std::size_t i = 0; i != N; ++i)
      v[i] = random<T>::get(a, b);
    return v;
  }
  static Vec<N,T> get() {
    return get(T(0), T(1));
  }
};
} // end namespace fmmtl
