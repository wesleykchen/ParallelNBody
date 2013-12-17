#pragma once

#include <iostream>

#include "meta/func_traits.hpp"

template <typename Kernel>
struct KernelTraits {
 public:
  typedef KernelTraits<Kernel>                    self_type;
  typedef Kernel                                  kernel_type;

  typedef typename kernel_type::source_type       source_type;
  typedef typename kernel_type::target_type       target_type;
  typedef typename kernel_type::charge_type       charge_type;
  typedef typename kernel_type::kernel_value_type kernel_value_type;
  typedef typename kernel_type::result_type       result_type;

  // Kernel evaluation operator, K(t,s)
  HAS_MEM_FUNC(HasEvalOp,
               kernel_value_type, operator(),
               const target_type&, const source_type&);
  static const bool has_eval_op = HasEvalOp<Kernel>::value;
  // Kernel transpose, K.transpose(K(t,s))
  HAS_MEM_FUNC(HasTranspose,
               kernel_value_type, transpose,
               const kernel_value_type&);
  static const bool has_transpose = HasTranspose<Kernel>::value;

  friend std::ostream& operator<<(std::ostream& s, const self_type& traits) {
    s << "has_eval_op: "           << traits.has_eval_op           << std::endl;
    s << "has_transpose: "         << traits.has_transpose         << std::endl;
    return s;
  }
};

#define IMPORT_KERNEL_TRAITS(K)                                         \
  typedef typename KernelTraits<K>::kernel_type        kernel_type;     \
  typedef typename KernelTraits<K>::kernel_value_type  kernel_value_type; \
  typedef typename KernelTraits<K>::source_type        source_type;     \
  typedef typename KernelTraits<K>::target_type        target_type;     \
  typedef typename KernelTraits<K>::charge_type        charge_type;     \
  typedef typename KernelTraits<K>::result_type        result_type
