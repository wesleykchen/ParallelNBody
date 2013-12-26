#pragma once
/** @file P2P.hpp
 * @brief Dispatch methods for the P2P stage
 *
 */

#include <iterator>
#include <type_traits>

#include "meta/kernel_traits.hpp"
#include "meta/trivial_iterator.hpp"

namespace detail {

/** Dual-Evaluation dispatch when K.transpose does not exist
 */
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

/** Dual-Evaluation dispatch when K.transpose does exist
 */
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

/** Asymmetric block P2P
 * r_i += sum_j K(t_i, s_j) * c_j
 *
 * @param[in] ...
 */
template <typename Kernel,
          typename SourceIter, typename ChargeIter,
          typename TargetIter, typename ResultIter>
inline void
block_eval(const Kernel& K,
           SourceIter s_first, SourceIter s_last, ChargeIter c_first,
           TargetIter t_first, TargetIter t_last, ResultIter r_first)
{
  std::cout << "Iterator Version!" << std::endl;

  typedef typename std::iterator_traits<SourceIter>::value_type source_type;
  typedef typename std::iterator_traits<TargetIter>::value_type target_type;
  typedef typename std::iterator_traits<ChargeIter>::value_type charge_type;
  typedef typename std::iterator_traits<ResultIter>::value_type result_type;

  // TODO: Add static_asserts with improved KernelTraits for improved errors

  for ( ; t_first != t_last; ++t_first, ++r_first) {
    const target_type& t = *t_first;
    result_type& r       = *r_first;

    SourceIter si = s_first;
    ChargeIter ci = c_first;
    for ( ; si != s_last; ++si, ++ci)
      r += K(t,*si) * (*ci);
  }
}

/** Asymmetric block P2P optimized for pointers-to-data
 * r_i += sum_j K(t_i, s_j) * c_j
 *
 * @note This method attempts to take advantage of cache optimization (TODO)
 * @param[in] ...
 */
template <typename Kernel,
          typename Source, typename Charge,
          typename Target, typename Result>
inline void
block_eval(const Kernel& K,
           Source* s_first, Source* s_last, Charge* c_first,
           Target* t_first, Target* t_last, Result* r_first) {
  std::cout << "Pointer Version!" << std::endl;

  for ( ; t_first != t_last; ++t_first, ++r_first) {
    const Target& t = *t_first;
    Result& r       = *r_first;

    const Source* si = s_first;
    const Charge* ci = c_first;
    for ( ; si != s_last; ++si, ++ci)
      r += K(t,*si) * (*ci);
  }
}


/** Symmetric off-diagonal block P2P using the evaluation operator
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

/** Symmetric diagonal block P2P using the evaluation operator
 * r_i += sum_j K(p_i, p_j) * c_j
 *
 * @pre source_type == target_type
 */
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

  for (++ipi, ++ici, ++iri; ipi != p_last; ++ipi, ++ici, ++iri) {
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
block_eval(const Kernel& K,
           SourceIter s_first, SourceIter s_last, ChargeIter c_first,
           TargetIter t_first, TargetIter t_last, ResultIter r_first)
{
  return detail::block_eval(K,
                            iter_base(s_first), iter_base(s_last),
                            iter_base(c_first),
                            iter_base(t_first), iter_base(t_last),
                            iter_base(r_first));
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
          typename SourceIter, typename ChargeIter, typename ResultIter>
inline void
block_eval(const Kernel& K,
           SourceIter p1_first, SourceIter p1_last,
           ChargeIter c1_first, ResultIter r1_first,
           SourceIter p2_first, SourceIter p2_last,
           ChargeIter c2_first, ResultIter r2_first)
{
  return detail::block_eval(K,
                            iter_base(p1_first), iter_base(p1_last),
                            iter_base(c1_first), iter_base(r1_first),
                            iter_base(p2_first), iter_base(p2_last),
                            iter_base(c2_first), iter_base(r2_first));
}

/** Symmetric diagonal block P2P
 * r_i += sum_j K(p_i, p_j) * c_j
 *
 * @pre source_type == target_type
 */
template <typename Kernel,
          typename SourceIter, typename ChargeIter, typename ResultIter>
inline static void
block_eval(const Kernel& K,
           SourceIter p_first, SourceIter p_last,
           ChargeIter c_first, ResultIter r_first)
{
  return detail::block_eval(K,
                            iter_base(p_first), iter_base(p_last),
                            iter_base(c_first), iter_base(r_first));
}
