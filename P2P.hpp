#pragma once
/** @file P2P.hpp
 * @brief Dispatch methods for the P2P stage
 *
 */

#include <iterator>
#include <type_traits>
#include <thread>

#include "meta/kernel_traits.hpp"
#include "meta/trivial_iterator.hpp"

#if !defined(P2P_BLOCK_SIZE)
#  define P2P_BLOCK_SIZE 128
#endif

namespace detail {

/** Dual-Evaluation dispatch when K.transpose does not exist */
template <typename Kernel,
          typename Source, typename Charge,
          typename Target, typename Result>
inline
typename std::enable_if<!KernelTraits<Kernel>::has_transpose>::type
symm_eval(const Kernel& K,
          const Source& p1, const Charge& c1, Result& r1,
          const Target& p2, const Charge& c2, Result& r2)
{
  r1 += K(p1,p2) * c2;
  r2 += K(p2,p1) * c1;
}

/** Dual-Evaluation dispatch when K.transpose does exist */
template <typename Kernel,
          typename Source, typename Charge,
          typename Target, typename Result>
inline
typename std::enable_if<KernelTraits<Kernel>::has_transpose>::type
symm_eval(const Kernel& K,
          const Source& p1, const Charge& c1, Result& r1,
          const Target& p2, const Charge& c2, Result& r2)
{
  typedef typename kernel_value<Kernel>::type kernel_value_type;

  kernel_value_type k12 = K(p1,p2);
  r1 += k12 * c2;
  r2 += K.transpose(k12) * c1;
}

/** Asymmetric block P2P evaluation */
template <typename Kernel,
          typename SourceIter, typename ChargeIter,
          typename TargetIter, typename ResultIter>
inline void
block_eval(const Kernel& K,
           SourceIter s_first, SourceIter s_last, ChargeIter c_first,
           TargetIter t_first, TargetIter t_last, ResultIter r_first)
{
  typedef typename std::iterator_traits<SourceIter>::value_type source_type;
  typedef typename std::iterator_traits<TargetIter>::value_type target_type;
  typedef typename std::iterator_traits<ChargeIter>::value_type charge_type;
  typedef typename std::iterator_traits<ResultIter>::value_type result_type;

  for ( ; t_first != t_last; ++t_first, ++r_first) {
    const target_type& t = *t_first;
    result_type& r       = *r_first;

    SourceIter si = s_first;
    ChargeIter ci = c_first;
    for ( ; si != s_last; ++si, ++ci)
      r += K(t,*si) * (*ci);
  }
}

/** Symmetric off-diagonal block P2P evaluation */
template <typename Kernel,
          typename SourceIter, typename ChargeIter,
          typename TargetIter, typename ResultIter>
inline void
block_eval(const Kernel& K,
           SourceIter p1_first, SourceIter p1_last,
           ChargeIter c1_first, ResultIter r1_first,
           TargetIter p2_first, TargetIter p2_last,
           ChargeIter c2_first, ResultIter r2_first)
{
  typedef typename std::iterator_traits<SourceIter>::value_type source_type;
  typedef typename std::iterator_traits<TargetIter>::value_type target_type;
  typedef typename std::iterator_traits<ChargeIter>::value_type charge_type;
  typedef typename std::iterator_traits<ResultIter>::value_type result_type;

  for ( ; p1_first != p1_last; ++p1_first, ++c1_first, ++r1_first) {
    const source_type& pi = *p1_first;
    const charge_type& ci = *c1_first;
    result_type& ri       = *r1_first;

    TargetIter p2i = p2_first;
    ChargeIter c2i = c2_first;
    ResultIter r2i = r2_first;
    for ( ; p2i != p2_last; ++p2i, ++c2i, ++r2i) {
      const target_type& pj = *p2i;
      const charge_type& cj = *c2i;
      result_type& rj       = *r2i;

      symm_eval(K, pi, ci, ri, pj, cj, rj);
    }
  }
}

/** Symmetric diagonal block P2P evaluation */
template <typename Kernel,
          typename SourceIter, typename ChargeIter, typename ResultIter>
inline static void
block_eval(const Kernel& K,
           SourceIter p_first, SourceIter p_last,
           ChargeIter c_first, ResultIter r_first)
{
  typedef typename std::iterator_traits<SourceIter>::value_type source_type;
  typedef typename std::iterator_traits<ChargeIter>::value_type charge_type;
  typedef typename std::iterator_traits<ResultIter>::value_type result_type;

  SourceIter ipi = p_first;
  ChargeIter ici = c_first;
  ResultIter iri = r_first;

  for ( ; ipi != p_last; ++ipi, ++ici, ++iri) {
    const source_type& pi = *ipi;
    const charge_type& ci = *ici;
    result_type& ri       = *iri;

    // The diagonal element
    ri += K(pi,pi) * ci;

    // The off-diagonal elements
    SourceIter ipj = p_first;
    ChargeIter icj = c_first;
    ResultIter irj = r_first;
    for ( ; ipj != ipi; ++ipj, ++icj, ++irj) {
      const source_type& pj = *ipj;
      const charge_type& cj = *icj;
      result_type& rj       = *irj;

      symm_eval(K, pi, ci, ri, pj, cj, rj);
    }
  }
}


/*************************************/
/****** Dispatch Methods *************/
/*************************************/


/** Asymmetric block P2P */
template <typename Kernel,
          typename SourceIter, typename ChargeIter,
          typename TargetIter, typename ResultIter>
inline void
p2p(const Kernel& K,
    SourceIter s_first, SourceIter s_last, ChargeIter c_first,
    TargetIter t_first, TargetIter t_last, ResultIter r_first)
{
  return block_eval(K,
                    s_first, s_last, c_first,
                    t_first, t_last, r_first);
}

/** Symmetric off-diagonal block P2P */
template <typename Kernel,
          typename SourceIter, typename ChargeIter,
          typename TargetIter, typename ResultIter>
inline void
p2p(const Kernel& K,
    SourceIter p1_first, SourceIter p1_last,
    ChargeIter c1_first, ResultIter r1_first,
    TargetIter p2_first, TargetIter p2_last,
    ChargeIter c2_first, ResultIter r2_first)
{
  return block_eval(K,
                    p1_first, p1_last,
                    c1_first, r1_first,
                    p2_first, p2_last,
                    c2_first, r2_first);
}

/** Symmetric diagonal block P2P */
template <typename Kernel,
          typename SourceIter, typename ChargeIter, typename ResultIter>
inline void
p2p(const Kernel& K,
    SourceIter p_first, SourceIter p_last,
    ChargeIter c_first, ResultIter r_first)
{
  return block_eval(K,
                    p_first, p_last,
                    c_first, r_first);
}

/** Asymmetric block P2P optimized for pointers-to-data */
template <typename Kernel,
          typename Source, typename Charge,
          typename Target, typename Result>
inline void
p2p(const Kernel& K,
    Source* s_first, Source* s_last, Charge* c_first,
    Target* t_first, Target* t_last, Result* r_first,
    unsigned threads = std::thread::hardware_concurrency())
{
  const int count1 = s_last - s_first;
  const int count2 = t_last - t_first;

  const char flag = ((count1 > P2P_BLOCK_SIZE) << 1) | (count2 > P2P_BLOCK_SIZE);
  switch (flag) {
    case 0: { // Both are small, evaluate
      block_eval(K, s_first, s_last, c_first,
                    t_first, t_last, r_first);
    } break;
    case 1: { // Split the targets
      Target* t_half = t_first + count2/2;
      Result* r_half = r_first + count2/2;

      if (threads > 0) {
        // In parallel
        std::thread thr([=](){
        p2p(K, s_first, s_last, c_first,
               t_first, t_half, r_first, threads-1);
          });
        p2p(K, s_first, s_last, c_first,
               t_half, t_last, r_half, threads-1);
        thr.join();
      } else {
        p2p(K, s_first, s_last, c_first,
               t_first, t_half, r_first, threads);
        p2p(K, s_first, s_last, c_first,
               t_half, t_last, r_half, threads);
      }
    } break;
    case 2: { // Split the sources
      Source* s_half = s_first + count1/2;
      Charge* c_half = c_first + count1/2;
      p2p(K, s_first, s_half, c_first,
             t_first, t_last, r_first, threads);
      p2p(K, s_half,  s_last, c_half,
             t_first, t_last, r_first, threads);
    } break;
    case 3: { // Split both
      Source* s_half = s_first + count1/2;
      Charge* c_half = c_first + count1/2;
      Target* t_half = t_first + count2/2;
      Result* r_half = r_first + count2/2;

      if (threads > 0) {
        // Top and bottom in parallel
        std::thread thr([=](){
        p2p(K, s_first, s_half, c_first,
               t_first, t_half, r_first, threads-1);
        p2p(K, s_half,  s_last, c_half,
               t_first, t_half, r_first, threads-1);
          });
        p2p(K, s_first, s_half, c_first,
               t_half,  t_last, r_half, threads-1);
        p2p(K, s_half, s_last, c_half,
               t_half, t_last, r_half, threads-1);
        thr.join();
      } else {
        p2p(K, s_first, s_half, c_first,
               t_first, t_half, r_first, threads);
        p2p(K, s_half,  s_last, c_half,
               t_first, t_half, r_first, threads);
        p2p(K, s_first, s_half, c_first,
               t_half,  t_last, r_half, threads);
        p2p(K, s_half, s_last, c_half,
               t_half, t_last, r_half, threads);
      }
    } break;
  }
}

/** Symmetric off-diagonal block P2P optimized for pointer-to-data */
template <typename Kernel,
          typename Source, typename Charge,
          typename Target, typename Result>
inline void
p2p(const Kernel& K,
    Source* p1_first, Source* p1_last,
    Charge* c1_first, Result* r1_first,
    Target* p2_first, Target* p2_last,
    Charge* c2_first, Result* r2_first,
    unsigned threads = std::thread::hardware_concurrency())
{
  const int count1 = p1_last - p1_first;
  const int count2 = p2_last - p2_first;

  const char flag = ((count1 > P2P_BLOCK_SIZE) << 1) | (count2 > P2P_BLOCK_SIZE);
  switch (flag) {
    case 0: { // Both are small, evaluate
      block_eval(K, p1_first, p1_last, c1_first, r1_first,
                    p2_first, p2_last, c2_first, r2_first);
    } break;
    case 1: { // Split the p2
      Target* p2_half = p2_first + count2/2;
      Charge* c2_half = c2_first + count2/2;
      Result* r2_half = r2_first + count2/2;
      p2p(K, p1_first, p1_last, c1_first, r1_first,
             p2_first, p2_half, c2_first, r2_first, threads);
      p2p(K, p1_first, p1_last, c1_first, r1_first,
             p2_half,  p2_last, c2_half, r2_half, threads);
    } break;
    case 2: { // Split the p1
      Source* p1_half = p1_first + count2/2;
      Charge* c1_half = c1_first + count2/2;
      Result* r1_half = r1_first + count2/2;
      p2p(K, p1_first, p1_half, c1_first, r1_first,
             p2_first, p2_last, c2_first, r2_first, threads);
      p2p(K, p1_half,  p1_last, c1_half,  r1_half,
             p2_first, p2_last, c2_first, r2_first, threads);
    } break;
    case 3: { // Split both
      Source* p1_half = p1_first + count1/2;
      Charge* c1_half = c1_first + count1/2;
      Result* r1_half = r1_first + count1/2;
      Target* p2_half = p2_first + count2/2;
      Charge* c2_half = c2_first + count2/2;
      Result* r2_half = r2_first + count2/2;

      if (threads > 0) {
        // Upper left and bottom right in parallel
        std::thread thr1([=](){
        p2p(K, p1_first, p1_half, c1_first, r1_first,
               p2_first, p2_half, c2_first, r2_first, threads-1);
          });
        p2p(K, p1_half, p1_last, c1_half, r1_half,
               p2_half, p2_last, c2_half, r2_half, threads-1);
        thr1.join();

        // Bottom left and top right in parallel
        std::thread thr2([=](){
        p2p(K, p1_half,  p1_last, c1_half,  r1_half,
               p2_first, p2_half, c2_first, r2_first, threads-1);
          });
        p2p(K, p1_first, p1_half, c1_first, r1_first,
               p2_half,  p2_last, c2_half,  r2_half, threads-1);
        thr2.join();
      } else {
        p2p(K, p1_first, p1_half, c1_first, r1_first,
               p2_first, p2_half, c2_first, r2_first, threads);
        p2p(K, p1_half,  p1_last, c1_half,  r1_half,
               p2_first, p2_half, c2_first, r2_first, threads);
        p2p(K, p1_first, p1_half, c1_first, r1_first,
               p2_half,  p2_last, c2_half,  r2_half, threads);
        p2p(K, p1_half, p1_last, c1_half, r1_half,
               p2_half, p2_last, c2_half, r2_half, threads);
      }
    } break;
  }
}

/** Symmetric diagonal block P2P optimized for pointer-to-data */
template <typename Kernel,
          typename Source, typename Charge, typename Result>
inline void
p2p(const Kernel& K,
    Source* p_first, Source* p_last,
    Charge* c_first, Result* r_first,
    unsigned threads = std::thread::hardware_concurrency())
{
  const int count = p_last - p_first;
  if (count > P2P_BLOCK_SIZE) {   // TODO: Generalize
    Source* p_half = p_first + count/2;
    Charge* c_half = c_first + count/2;
    Result* r_half = r_first + count/2;

    if (threads > 0) {
      // Two symmetric diagonal blocks in parallel
      std::thread thr([=](){
      p2p(K, p_first, p_half, c_first, r_first, threads-1);
        });
      p2p(K, p_half,  p_last, c_half,  r_half, threads-1);
      thr.join();
      // Symmetric off-diagonal block
      p2p(K, p_first, p_half, c_first, r_first,
             p_half,  p_last, c_half,  r_half, threads);
    } else {
      p2p(K, p_first, p_half, c_first, r_first, threads);
      p2p(K, p_first, p_half, c_first, r_first,
             p_half,  p_last, c_half,  r_half, threads);
      p2p(K, p_half,  p_last, c_half,  r_half, threads);
    }
  } else {
    block_eval(K, p_first, p_last, c_first, r_first);
  }
}

} // end namespace detail




/** Asymmetric block P2P
 * r_i += sum_j K(t_i, s_j) * c_j
 *
 * @param[in] ...
 */
template <typename Kernel,
          typename SourceIter, typename ChargeIter,
          typename TargetIter, typename ResultIter>
inline void
p2p(const Kernel& K,
    SourceIter s_first, SourceIter s_last, ChargeIter c_first,
    TargetIter t_first, TargetIter t_last, ResultIter r_first,
    unsigned threads = std::thread::hardware_concurrency())
{
  return detail::p2p(K,
                     iter_base(s_first), iter_base(s_last),
                     iter_base(c_first),
                     iter_base(t_first), iter_base(t_last),
                     iter_base(r_first),
                     threads);
}

/** Symmetric off-diagonal block P2P
 * r2_i += sum_j K(p2_i, p1_j) * c1_j
 * r1_j += sum_i K(p1_j, p2_i) * c2_i
 *
 * @param[in] ...
 * @pre source_type == target_type
 * @pre For all i,j we have p1_i != p2_j
 */
template <typename Kernel,
          typename SourceIter, typename ChargeIter,
          typename TargetIter, typename ResultIter>
inline void
p2p(const Kernel& K,
    SourceIter p1_first, SourceIter p1_last,
    ChargeIter c1_first, ResultIter r1_first,
    TargetIter p2_first, TargetIter p2_last,
    ChargeIter c2_first, ResultIter r2_first,
    unsigned threads = std::thread::hardware_concurrency())
{
  return detail::p2p(K,
                     iter_base(p1_first), iter_base(p1_last),
                     iter_base(c1_first), iter_base(r1_first),
                     iter_base(p2_first), iter_base(p2_last),
                     iter_base(c2_first), iter_base(r2_first),
                     threads);
}

/** Symmetric diagonal block P2P
 * r_i += sum_j K(p_i, p_j) * c_j
 *
 * @pre source_type == target_type
 */
template <typename Kernel,
          typename SourceIter, typename ChargeIter, typename ResultIter>
inline void
p2p(const Kernel& K,
    SourceIter p_first, SourceIter p_last,
    ChargeIter c_first, ResultIter r_first,
    unsigned threads = std::thread::hardware_concurrency())
{
  return detail::p2p(K,
                     iter_base(p_first), iter_base(p_last),
                     iter_base(c_first), iter_base(r_first),
                     threads);
}
