#pragma once
/** @file Vec.hpp
 * @brief A wrapper class to provide simple vector operations for
 * primitive types and classes that only require op[].
 */

#include <iostream>
#include <cmath>

#define for_i for(unsigned i=0; i!=N; ++i)

/** @class Vec
 * @brief Class representing ND points and vectors.
 */
template <unsigned N, typename T = double>
class Vec {
  T a[N];

 public:
  typedef T         value_type;
  typedef T*        pointer;
  typedef const T*  const_pointer;
  typedef T&        reference;
  typedef const T&  const_reference;

  typedef T*        iterator;
  typedef const T*  const_iterator;

  // CONSTRUCTORS

  inline Vec() {
    for_i a[i] = value_type();
  }
  inline explicit Vec(value_type b) {
    for_i a[i] = b;
  }
  inline Vec(value_type b0, value_type b1) {
    assert(N == 2);
    a[0] = b0; a[1] = b1;
  }
  inline Vec(value_type b0, value_type b1, value_type b2) {
    assert(N == 3);
    a[0] = b0; a[1] = b1; a[2] = b2;
  }
  inline Vec(value_type b0, value_type b1, value_type b2, value_type b3) {
    assert(N == 4);
    a[0] = b0; a[1] = b1; a[2] = b2; a[3] = b3;
  }

  // COMPARATORS

  inline bool operator==(const Vec& b) const {
    for_i if (a[i] != b[i]) return false;
    return true;
  }
  inline bool operator!=(const Vec& b) const {
    return !(*this == b);
  }

  // MODIFIERS

  /** Return a negated version of @a p. */
  inline Vec operator-() const {
    Vec<N,T> m = *this;
    for_i m[i] = -m[i];
    return m;
  }
  /** Add scalar @a b to this Vec */
  template <typename D>
  inline Vec& operator+=(const D& b) {
    for_i a[i] += b;
    return *this;
  }
  /** Subtract scalar @a b from this Vec */
  template <typename D>
  inline Vec& operator-=(const D& b) {
    for_i a[i] -= b;
    return *this;
  }
  /** Scale this Vec up by scalar @a b */
  template <typename D>
  inline Vec& operator*=(const D& b) {
    for_i a[i] *= b;
    return *this;
  }
  /** Scale this Vec down by scalar @a b */
  template <typename D>
  inline Vec& operator/=(const D& b) {
    for_i a[i] /= b;
    return *this;
  }
  /** Add Vec @a b to this Vec */
  inline Vec& operator+=(const Vec& b) {
    for_i a[i] += b[i];
    return *this;
  }
  /** Subtract Vec @a b from this Vec */
  inline Vec& operator-=(const Vec& b) {
    for_i a[i] -= b[i];
    return *this;
  }
  /** Scale this Vec up by factors in @a b */
  inline Vec& operator*=(const Vec& b) {
    for_i a[i] *= b[i];
    return *this;
  }
  /** Scale this Vec down by factors in @a b */
  inline Vec& operator/=(const Vec& b) {
    for_i a[i] /= b[i];
    return *this;
  }
  /** Compute the dot product of this Vec with another Vec */
  inline value_type dot(const Vec& b) const {
    value_type d = value_type();
    for_i d += a[i]*b[i];
    return d;
  }

  // ACCESSORS

  /** Access the @a i th (lvalue) element of this Vec
   * @pre i < size() */
  inline value_type& operator[](unsigned i) {
    return a[i];
  }
  /** Access the @a i th (rvalue) element of this Vec
   * @pre i < size() */
  inline const value_type& operator[](unsigned i) const {
    return a[i];
  }
  inline static unsigned size() {
    return N;
  }
  inline iterator begin() {
    return a;
  }
  inline const_iterator begin() const {
    return a;
  }
  inline iterator end() {
    return a + N;
  }
  inline const_iterator end() const {
    return a + N;
  }
};


// OPERATORS

/** Write a Vec to an output stream */
template <unsigned N, typename P>
inline std::ostream& operator<<(std::ostream& s, const Vec<N,P>& a) {
  for_i s << a[i] << " ";
  return s;
}
/** Read a Vec from an input stream */
template <unsigned N, typename P>
inline std::istream& operator>>(std::istream& s, Vec<N,P>& a) {
  for_i s >> a[i];
  return s;
}

/** Compute the dot product of two Vecs */
template <unsigned N, typename P>
inline typename Vec<N,P>::value_type dot(const Vec<N,P>& a,
                                               const Vec<N,P>& b) {
  return a.dot(b);
}
/** Compute the dot product of two Vecs */
template <unsigned N, typename P>
inline typename Vec<N,P>::value_type inner_prod(const Vec<N,P>& a,
                                                      const Vec<N,P>& b) {
  return a.dot(b);
}
/** Compute cross product of two 3D Vecs */
template <typename P>
inline Vec<3,P> cross(const Vec<3,P>& a, const Vec<3,P>& b) {
  return Vec<3,P>(a[1]*b[2] - a[2]*b[1],
                  a[2]*b[0] - a[0]*b[2],
                  a[0]*b[1] - a[1]*b[0]);
}
/** Compute the squared L2 norm */
template <unsigned N, typename P>
inline typename Vec<N,P>::value_type normSq(const Vec<N,P>& a) {
  return a.dot(a);
}
/** Compute the L2 norm */
template <unsigned N, typename P>
inline typename Vec<N,P>::value_type norm(const Vec<N,P>& a) {
  using std::sqrt;
  return sqrt(normSq(a));
}
/** Compute the L2 norm */
template <unsigned N, typename P>
inline typename Vec<N,P>::value_type norm_2(const Vec<N,P>& a) {
  return norm(a);
}
/** Compute the L1 norm */
template <unsigned N, typename P>
inline typename Vec<N,P>::value_type norm_1(const Vec<N,P>& a) {
  using std::abs;
  typename Vec<N,P>::value_type r = typename Vec<N,P>::value_type();
  for_i r += abs(a[i]);
  return r;
}
/** Compute the L-infinity norm */
template <unsigned N, typename P>
inline typename Vec<N,P>::value_type norm_inf(const Vec<N,P>& a) {
  using std::abs;
  using std::max;
  typename Vec<N,P>::value_type a_max = typename Vec<N,P>::value_type();
  for_i a_max = max(a_max, abs(a[i]));
  return a_max;
}

// ARITHMETIC

/** Unary plus: Return @a p. ("+p" should work if "-p" works.) */
template <unsigned N, typename P>
inline Vec<N,P> operator+(const Vec<N,P>& a) {
  return a;
}
template <unsigned N, typename P>
inline Vec<N,P> operator+(Vec<N,P> a, const Vec<N,P>& b) {
  return a += b;
}
template <unsigned N, typename P, typename D>
inline Vec<N,P> operator+(Vec<N,P> a, const D& b) {
  return a += b;
}
template <unsigned N, typename P, typename D>
inline Vec<N,P> operator+(const D& b, Vec<N,P> a) {
  return a += b;
}
template <unsigned N, typename P>
inline Vec<N,P> operator-(Vec<N,P> a, const Vec<N,P>& b) {
  return a -= b;
}
template <unsigned N, typename P, typename D>
inline Vec<N,P> operator-(Vec<N,P> a, const D& b) {
  return a -= b;
}
template <unsigned N, typename P, typename D>
inline Vec<N,P> operator-(const D& b, const Vec<N,P>& a) {
  return (-a) += b;
}
template <unsigned N, typename P>
inline Vec<N,P> operator*(Vec<N,P> a, const Vec<N,P>& b) {
  return a *= b;
}
template <unsigned N, typename P, typename D>
inline Vec<N,P> operator*(Vec<N,P> a, const D& b) {
  return a *= b;
}
template <unsigned N, typename P, typename D>
inline Vec<N,P> operator*(const D& b, Vec<N,P> a) {
  return a *= b;
}
template <unsigned N, typename P>
inline Vec<N,P> operator/(Vec<N,P> a, const Vec<N,P>& b) {
  return a /= b;
}
template <unsigned N, typename P, typename D>
inline Vec<N,P> operator/(Vec<N,P> a, const D& b) {
  return a /= b;
}

// ADDITIONAL OPERATORS (Helps with Vec<N,Vec<M,T>> norms, etc)
template <unsigned N, typename P>
inline Vec<N,P> abs(Vec<N,P> a) {
  using std::abs;
  for_i a[i] = abs(a[i]);
  return a;
}
template <unsigned N, typename P>
inline Vec<N,P> sqrt(Vec<N,P> a) {
  using std::sqrt;
  for_i a[i] = sqrt(a[i]);
  return a;
}

#undef for_i
