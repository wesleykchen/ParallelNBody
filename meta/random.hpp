#pragma once

namespace meta {

template <typename T>
struct random;

template <>
struct random<double> {
  static double get() {
    return ::drand48();
  }
};

template <>
struct random<float> {
  static float get() {
    return (float) random<double>::get();
  }
};

template <>
struct random<unsigned> {
  static unsigned get() {
    return rand();
  }
};

template <>
struct random<int> {
  static int get() {
    return rand();
  }
};

} // end namespace meta
