#pragma once

#include <vector>
#include <valarray>

#include <type_traits>

#if !defined(P2P_DECAY_ITERATOR)
#  define P2P_DECAY_ITERATOR 1
#endif

// If Iterator has a base, return it
// Otherwise, Iterator is returned untouched.
template <typename Iterator, bool decay>
struct _Iter_base {
  typedef Iterator iterator_type;
  static inline iterator_type base(Iterator it) {
    return it;
  }
};


#if __GNUC__
// forward declaration of gnu's __normal_iterator
namespace __gnu_cxx
{
template <typename Iterator, typename Container>
class __normal_iterator;
} // end __gnu_cxx

template <typename T>
struct is_gnu_normal_iterator
    : std::false_type {};

// catch gnu __normal_iterators
template <typename T, typename Container>
struct is_gnu_normal_iterator<__gnu_cxx::__normal_iterator<T, Container> >
    : std::true_type {};

template <typename Iterator>
struct _Iter_base<Iterator, true /* is_gnu_normal_iterator */> {
  typedef typename Iterator::iterator_type iterator_type;
  static inline iterator_type base(Iterator it) {
    return it.base();
  }
};
#endif // __GNUC__


#if _MSC_VER
// forward declaration of MSVC's "normal iterators"
namespace std
{
template <typename Value, typename Diff, typename Ptr, typename Ref>
struct _Ranit;
} // end std

// catch msvc _Ranit
// TODO: This needs te be is_same, not is_convertible?
template <typename Iterator>
struct is_convertible_to_msvc_Ranit
    : std::is_convertible<
        Iterator,
        std::_Ranit<
          typename std::iterator_value<Iterator>::type,
          typename std::iterator_difference<Iterator>::type,
          typename std::iterator_pointer<Iterator>::type,
          typename std::iterator_reference<Iterator>::type
        >
      > {};

// TODO: MSC _Iter_base<Iterator, true>
#endif // _MSC_VER


template <typename T>
struct is_normal_iterator
    : std::integral_constant<bool,
#if __GNUC__
                             is_gnu_normal_iterator<T>::value
#endif // __GNUC__
#ifdef _MSC_VER
                             is_convertible_to_msvc_Ranit<T>::value
#endif // _MSC_VER
                             > {};

template <typename T>
struct is_trivial_iterator
    : std::integral_constant<bool,
                             std::is_pointer<T>::value
                             | is_normal_iterator<T>::value
                             > {};


// If Iterator is a __normal_iterator return its base (a plain pointer,
// normally) otherwise return it untouched.  See copy, fill, ...
template <typename Iterator>
struct _iter_base
    : _Iter_base<Iterator, is_normal_iterator<Iterator>::value> {};


#if P2P_DECAY_ITERATOR == 1     // Decay trivial iterators to pointers
template <typename Iterator>
inline typename _iter_base<Iterator>::iterator_type
iter_base(Iterator it) {
  return _iter_base<Iterator>::base(it);
}
#else                           // Return the same iterator
template <typename Iterator>
inline Iterator
iter_base(Iterator it) {
  return it;
}
#endif
