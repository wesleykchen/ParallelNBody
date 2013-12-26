#pragma once

#include <iostream>

#include "meta/func_traits.hpp"

template <typename T>
struct Void {
  typedef void type;
};

/** @brief Get the target type of a kernel function.
 *
 * If typename K::target_type exists, use it
 * Else, use the type of the first argument to k(t,s)
 */
template <typename K, typename _ = void>
struct target {
  typedef typename function_traits<K>::template argument<0>::type type;
};
template <typename K>
struct target<K, typename Void<typename K::target_type>::type> {
  typedef typename K::target_type type;
};

/** @brief Get the source type of a kernel function.
 *
 * If typename K::source_type exists, use it
 * Else, use the type of the second argument to k(t,s)
 */
template <typename K, typename _ = void>
struct source {
  typedef typename function_traits<K>::template argument<1>::type type;
};
template <typename K>
struct source<K, typename Void<typename K::source_type>::type> {
  typedef typename K::source_type type;
};


/** @brief Get the value type of a kernel function
 *
 * If typename K::kernel_value_type exists, use it
 * Else, use the return type of K(t,s)
 */
template <typename K, typename _ = void>
struct kernel_value {
  typedef typename function_traits<K>::return_type type;
};
template <typename K>
struct kernel_value<K, typename Void<typename K::kernel_value_type>::type> {
  typedef typename K::kernel_value_type type;
};


template <typename Kernel>
struct KernelTraits {
  typedef Kernel kernel_type;

  // Define the source_type, target_type, and kernel_value_type
  typedef typename source<kernel_type>::type source_type;
  typedef typename target<kernel_type>::type target_type;
  typedef typename kernel_value<kernel_type>::type kernel_value_type;

  // Kernel evaluation operator, K(t,s)
  HAS_MEM_FUNC(HasEvalOp,
               kernel_value_type, operator(),
               const target_type&, const source_type&);
  static const bool has_eval_op = HasEvalOp<kernel_type>::value;
  // Kernel transpose, K.transpose(K(t,s))
  HAS_MEM_FUNC(HasTranspose,
               kernel_value_type, transpose,
               const kernel_value_type&);
  static const bool has_transpose = HasTranspose<kernel_type>::value;

  friend std::ostream& operator<<(std::ostream& s, const KernelTraits& traits) {
    s << "has_eval_op: "           << traits.has_eval_op           << std::endl;
    s << "has_transpose: "         << traits.has_transpose         << std::endl;
    return s;
  }
};
