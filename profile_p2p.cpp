#include "P2P.hpp"
#include "Util.hpp"
#include "meta/random.hpp"

#include "kernel/NormSq.kern"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <iterator>


template <typename T>
std::vector<T> generate(unsigned N) {
  std::vector<T> a;
  a.reserve(N);
  for ( ; N != 0; --N)
    a.push_back(meta::random<T>::get());
  return a;
}


int main(int argv, char** arg) {
  typedef NormSq kernel_type;
  kernel_type K;

  typedef kernel_type::source_type source_type;
  typedef kernel_type::charge_type charge_type;
  typedef kernel_type::target_type target_type;
  typedef kernel_type::result_type result_type;

  unsigned N = 1 << 17;
  char run_type = string_to_<char>(arg[1]);

  Clock timer;

  switch (run_type) {
    case 'D':

  //std::cout << "Symmetric Diagonal" << std::endl;
  for (unsigned n = 64; n < N; n *= 1.1) {
    std::vector<source_type> s = generate<source_type>(n);
    std::vector<charge_type> c = generate<charge_type>(n);
    std::vector<result_type> r = generate<result_type>(n);

    // Copy to prevent warm cache
    std::vector<source_type> s1 = s;
    std::vector<charge_type> c1 = c;
    std::vector<result_type> r1 = r;

    timer.start();
    detail::block_eval(K, s1.begin(), s1.end(), c1.begin(), r1.begin());
    double old_time = timer.elapsed();

    // Copy to prevent warm cache
    std::vector<source_type> s2 = s;
    std::vector<charge_type> c2 = c;
    std::vector<result_type> r2 = r;

    timer.start();
    p2p(K, s2.begin(), s2.end(), c2.begin(), r2.begin());
    double new_time = timer.elapsed();

    double error = 0;
    for (unsigned i = 0; i < n; ++i) {
      error += normSq(r2[i] - r1[i]) / normSq(r1[i]);
    }
    error = std::sqrt(error);

    std::cout << std::setw(10) << P2P_BLOCK_SIZE << "\t"
              << std::setw(10) << n << "\t"
              << std::setw(10) << error << "\t"
              << std::setw(10) << old_time << "\t"
              << std::setw(10) << new_time << "\t"
              << std::endl;

    //if (std::max(old_time, new_time) > 10) break;
  }

    break;

    case 'O':

  //std::cout << "Symmetric Off-Diagonal" << std::endl;
  for (unsigned n = 64; n < N; n *= 1.1) {
    std::vector<source_type> s = generate<source_type>(n);
    std::vector<charge_type> c = generate<charge_type>(n);
    std::vector<result_type> r = generate<result_type>(n);

    // Copy to prevent warm cache
    std::vector<source_type> s1 = s;
    std::vector<charge_type> c1 = c;
    std::vector<result_type> r1 = r;

    timer.start();
    detail::block_eval(K, s1.begin(), s1.begin()+n/2, c1.begin(), r1.begin(),
                          s1.begin()+n/2, s1.end(), c1.begin()+n/2, r1.begin()+n/2);
    double old_time = timer.elapsed();

    // Copy to prevent warm cache
    std::vector<source_type> s2 = s;
    std::vector<charge_type> c2 = c;
    std::vector<result_type> r2 = r;

    timer.start();
    p2p(K, s2.begin(), s2.begin()+n/2, c2.begin(), r2.begin(),
           s2.begin()+n/2, s2.end(), c2.begin()+n/2, r2.begin()+n/2);
    double new_time = timer.elapsed();

    double error = 0;
    for (unsigned i = 0; i < n; ++i) {
      error += normSq(r2[i] - r1[i]) / normSq(r1[i]);
    }
    error = std::sqrt(error);

    std::cout << std::setw(10) << P2P_BLOCK_SIZE << "\t"
              << std::setw(10) << n << "\t"
              << std::setw(10) << error << "\t"
              << std::setw(10) << old_time << "\t"
              << std::setw(10) << new_time << "\t"
              << std::endl;

    //if (std::max(old_time, new_time) > 10) break;
  }

  break;
    case 'A':

  //std::cout << "Asymmetric off-diagonal" << std::endl;
  for (unsigned n = 64; n < N; n *= 1.1) {
    std::vector<source_type> s = generate<source_type>(n);
    std::vector<target_type> t = generate<target_type>(n);
    std::vector<charge_type> c = generate<charge_type>(n);
    std::vector<result_type> r = generate<result_type>(n);

    // Copy to prevent warm cache
    std::vector<source_type> s1 = s;
    std::vector<source_type> t1 = t;
    std::vector<charge_type> c1 = c;
    std::vector<result_type> r1 = r;

    timer.start();
    detail::block_eval(K, s1.begin(), s1.end(), c1.begin(),
                          t1.begin(), t1.end(), r1.begin());
    double old_time = timer.elapsed();

    // Copy to prevent warm cache
    std::vector<source_type> s2 = s;
    std::vector<source_type> t2 = t;
    std::vector<charge_type> c2 = c;
    std::vector<result_type> r2 = r;

    timer.start();
    p2p(K, s2.begin(), s2.end(), c2.begin(),
           t2.begin(), t2.end(), r2.begin());
    double new_time = timer.elapsed();

    double error = 0;
    for (unsigned i = 0; i < n; ++i) {
      error += normSq(r2[i] - r1[i]) / normSq(r1[i]);
    }
    error = std::sqrt(error);

    std::cout << std::setw(10) << P2P_BLOCK_SIZE << "\t"
              << std::setw(10) << n << "\t"
              << std::setw(10) << error << "\t"
              << std::setw(10) << old_time << "\t"
              << std::setw(10) << new_time << "\t"
              << std::endl;

    //if (std::max(old_time, new_time) > 10) break;
  }

  }
}
