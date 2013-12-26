#pragma once

/* @class      : HAS_TYPEDEF
 * @brief      : This macro may be used to check if a class has a public typedef
 * @param NAME        : Name of struct this macro creates.
 * @param TYPEDEF     : Name of typedef to test for.
 * @note This macro is C++03 compatible
 */
#define HAS_TYPEDEF(NAME, TYPEDEF)                                      \
  template <typename CLASS>                                             \
  struct NAME {                                                         \
    template <typename U> static char chk(U::TYPEDEF*);                 \
    template <typename  > static long chk(...);                         \
    static const bool value = sizeof(chk<CLASS>(0)) == sizeof(char);    \
  }

/* @class      : HAS_MEM_FUNC
 * @brief      : This macro may be used to check if a class has a public
 *               const member function with particular signature.
 * @param NAME        : Name of struct this macro creates.
 * @param RETURN_TYPE : Return type of the member function to test for.
 * @param FUNC        : Name of member function to test for.
 * @param (...)       : The argument types of the member function to test for.
 *                      These complete the signature of the member funtion.
 * @note This macro is C++03 compatible
 */
#define HAS_MEM_FUNC(NAME, RETURN_TYPE, FUNC, ...)                      \
  template <typename CLASS>                                             \
  struct NAME {                                                         \
    typedef RETURN_TYPE (CLASS::*A)(__VA_ARGS__) const;                 \
    template <typename U, U> struct type_check;                         \
    template <typename U> static char chk(type_check<A, &U::FUNC>*);    \
    template <typename  > static long chk(...);                         \
    static bool const value = sizeof(chk<CLASS>(0)) == sizeof(char);    \
  }


#include <tuple>

template <class F>
struct function_traits;

// function pointer
template <class R, class... Args>
struct function_traits<R(*)(Args...)>
    : public function_traits<R(Args...)> {};

template <class R, class... Args>
struct function_traits<R(Args...)> {
  using return_type = R;

  static constexpr std::size_t arity = sizeof...(Args);

  template <std::size_t N>
  struct argument {
    static_assert(N < arity, "error: invalid parameter index.");
    using type = typename std::tuple_element<N,std::tuple<Args...>>::type;
  };
};

// member function pointer
template <class C, class R, class... Args>
struct function_traits<R(C::*)(Args...)>
    : public function_traits<R(C&,Args...)> {};

// const member function pointer
template <class C, class R, class... Args>
struct function_traits<R(C::*)(Args...) const>
    : public function_traits<R(C&,Args...)> {};

// member object pointer
template <class C, class R>
struct function_traits<R(C::*)>
    : public function_traits<R(C&)> {};

// functor
template <class F>
class function_traits {
  using call_type = function_traits<decltype(&F::operator())>;
 public:
  using return_type = typename call_type::return_type;

  static constexpr std::size_t arity = call_type::arity - 1;

  template <std::size_t N>
  struct argument {
    static_assert(N < arity, "error: invalid parameter index.");
    using type = typename call_type::template argument<N+1>::type;
  };
};

template <class F>
struct function_traits<F&>
    : public function_traits<F> {};

template <class F>
struct function_traits<F&&>
    : public function_traits<F> {};
