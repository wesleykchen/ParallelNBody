#pragma once
/** @file Vec.hpp
 * @brief A small N-dimensional numerical vector type that works on CPU and GPU.
 */

#include <iostream>
#include <cmath>

#define for_i for(std::size_t i = 0; i != N; ++i)

/** @class Vec
 * @brief Class representing ND points and vectors.
 */
template <std::size_t N, typename T = double>
struct Vec {
  T elem[N];

  typedef T               value_type;
  typedef T&              reference;
  typedef const T&        const_reference;
  typedef T*              iterator;
  typedef const T*        const_iterator;
  typedef std::size_t     size_type;
  typedef std::ptrdiff_t  difference_type;

  // CONSTRUCTORS

  inline Vec() {
    for_i elem[i] = value_type();
  }
  // TODO: Force 0-initialization to get POD and trivial semantics?
  //inline Vec() = default;
  inline explicit Vec(value_type b) {
    for_i elem[i] = b;
  }
  inline Vec(value_type b0, value_type b1) {
    static_assert(N >= 2, "Too many arguments in Vec constructor");
    elem[0] = b0; elem[1] = b1;
    for(std::size_t i = 2; i != N; ++i) elem[i] = value_type();
  }
  inline Vec(value_type b0, value_type b1, value_type b2) {
    static_assert(N >= 3, "Too many arguments in Vec constructor");
    elem[0] = b0; elem[1] = b1; elem[2] = b2;
    for(std::size_t i = 3; i != N; ++i) elem[i] = value_type();
  }
  inline Vec(value_type b0, value_type b1, value_type b2, value_type b3) {
    static_assert(N >= 4, "Too many arguments in Vec constructor");
    elem[0] = b0; elem[1] = b1; elem[2] = b2; elem[3] = b3;
    for(std::size_t i = 4; i != N; ++i) elem[i] = value_type();
  }

  // COMPARATORS

  inline bool operator==(const Vec& b) const {
    for_i if (elem[i] != b[i]) return false;
    return true;
  }
  inline bool operator!=(const Vec& b) const {
    return !(*this == b);
  }

  // MODIFIERS

  /** Add scalar @a b to this Vec */
  template <typename D>
  inline Vec& operator+=(const D& b) {
    for_i elem[i] += b;
    return *this;
  }
  /** Subtract scalar @a b from this Vec */
  template <typename D>
  inline Vec& operator-=(const D& b) {
    for_i elem[i] -= b;
    return *this;
  }
  /** Scale this Vec up by scalar @a b */
  template <typename D>
  inline Vec& operator*=(const D& b) {
    for_i elem[i] *= b;
    return *this;
  }
  /** Scale this Vec down by scalar @a b */
  template <typename D>
  inline Vec& operator/=(const D& b) {
    for_i elem[i] /= b;
    return *this;
  }
  /** Add Vec @a b to this Vec */
  inline Vec& operator+=(const Vec& b) {
    for_i elem[i] += b[i];
    return *this;
  }
  /** Subtract Vec @a b from this Vec */
  inline Vec& operator-=(const Vec& b) {
    for_i elem[i] -= b[i];
    return *this;
  }
  /** Scale this Vec up by factors in @a b */
  inline Vec& operator*=(const Vec& b) {
    for_i elem[i] *= b[i];
    return *this;
  }
  /** Scale this Vec down by factors in @a b */
  inline Vec& operator/=(const Vec& b) {
    for_i elem[i] /= b[i];
    return *this;
  }

  // ACCESSORS

  inline reference       operator[](size_type i)       { return elem[i]; }
  inline const_reference operator[](size_type i) const { return elem[i]; }

  inline T*       data()       { return elem; }
  inline const T* data() const { return elem; }

  inline reference       front()       { return elem[0]; }
  inline const_reference front() const { return elem[0]; }
  inline reference        back()       { return elem[N-1]; }
  inline const_reference  back() const { return elem[N-1]; }

  inline static size_type     size() { return N; }
  inline static size_type max_size() { return N; }
  inline static bool         empty() { return false; }

  // ITERATORS

  inline iterator        begin()       { return elem; }
  inline const_iterator  begin() const { return elem; }
  inline const_iterator cbegin() const { return elem; }

  inline iterator          end()       { return elem+N; }
  inline const_iterator    end() const { return elem+N; }
  inline const_iterator   cend() const { return elem+N; }
};

// OPERATORS

/** Write a Vec to an output stream */
template <std::size_t N, typename T>
inline std::ostream& operator<<(std::ostream& s, const Vec<N,T>& a) {
  for_i s << a[i] << " ";
  return s;
}
/** Read a Vec from an input stream */
template <std::size_t N, typename T>
inline std::istream& operator>>(std::istream& s, Vec<N,T>& a) {
  for_i s >> a[i];
  return s;
}

/** Compute cross product of two 3D Vecs */
template <typename T>
inline Vec<3,T> cross(const Vec<3,T>& a, const Vec<3,T>& b) {
  return Vec<3,T>(a[1]*b[2] - a[2]*b[1],
                  a[2]*b[0] - a[0]*b[2],
                  a[0]*b[1] - a[1]*b[0]);
}

// ARITHMETIC OPERATORS

/** Unary negation: Return -@a a */
template <std::size_t N, typename T>
inline Vec<N,T> operator-(Vec<N,T> a) {
  for_i a[i] = -a[i];
  return a;
}
/** Unary plus: Return @a a. ("+a" should work if "-a" works.) */
template <std::size_t N, typename T>
inline Vec<N,T> operator+(const Vec<N,T>& a) {
  return a;
}
template <std::size_t N, typename T>
inline Vec<N,T> operator+(Vec<N,T> a, const Vec<N,T>& b) {
  return a += b;
}
template <std::size_t N, typename T, typename D>
inline Vec<N,T> operator+(Vec<N,T> a, const D& b) {
  return a += b;
}
template <std::size_t N, typename T, typename D>
inline Vec<N,T> operator+(const D& b, Vec<N,T> a) {
  return a += b;
}
template <std::size_t N, typename T>
inline Vec<N,T> operator-(Vec<N,T> a, const Vec<N,T>& b) {
  return a -= b;
}
template <std::size_t N, typename T, typename D>
inline Vec<N,T> operator-(Vec<N,T> a, const D& b) {
  return a -= b;
}
template <std::size_t N, typename T, typename D>
inline Vec<N,T> operator-(const D& b, const Vec<N,T>& a) {
  return (-a) += b;
}
template <std::size_t N, typename T>
inline Vec<N,T> operator*(Vec<N,T> a, const Vec<N,T>& b) {
  return a *= b;
}
template <std::size_t N, typename T, typename D>
inline Vec<N,T> operator*(Vec<N,T> a, const D& b) {
  return a *= b;
}
template <std::size_t N, typename T, typename D>
inline Vec<N,T> operator*(const D& b, Vec<N,T> a) {
  return a *= b;
}
template <std::size_t N, typename T>
inline Vec<N,T> operator/(Vec<N,T> a, const Vec<N,T>& b) {
  return a /= b;
}
template <std::size_t N, typename T, typename D>
inline Vec<N,T> operator/(Vec<N,T> a, const D& b) {
  return a /= b;
}

// ELEMENTWISE OPERATORS

template <std::size_t N, typename T>
inline Vec<N,T> abs(Vec<N,T> a) {
  using std::abs;
  for_i a[i] = abs(a[i]);
  return a;
}
template <std::size_t N, typename T>
inline Vec<N,T> sqrt(Vec<N,T> a) {
  using std::sqrt;
  for_i a[i] = sqrt(a[i]);
  return a;
}

#undef for_i


#include "numeric/Norm.hpp"
