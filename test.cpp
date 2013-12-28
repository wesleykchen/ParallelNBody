#include "P2P.hpp"

#include "kernel/Laplace.kern"

int main() {
  typedef LaplacePotential kernel_type;
  kernel_type K;

  typedef kernel_type::source_type source_type;
  typedef kernel_type::charge_type charge_type;
  typedef kernel_type::target_type target_type;
  typedef kernel_type::result_type result_type;

  std::vector<source_type> s(3);
  std::vector<target_type> t(3);
  std::vector<charge_type> c(3);
  std::vector<result_type> r(3);

  p2p(K,
      s.begin(), s.end(), c.begin(),
      t.begin(), t.end(), r.begin());

  p2p(K,
      s.cbegin(), s.cend(), c.begin(),
      t.begin(), t.end(), r.begin());

  p2p([](const target_type&, const source_type&) { return 1; },
      s.begin(), s.end(), c.begin(), r.begin());
}
