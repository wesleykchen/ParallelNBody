#pragma once
/** @file P2P.hpp
 * @brief Dispatch methods for the P2P stage
 *
 */

#include <iterator>
#include <type_traits>

struct P2P
{
	/** Asymmetric block P2P using the evaluation operator
	 * r_i += sum_j K(t_i, s_j) * c_j
	 *
	 * @param[in] ...
	 */
	template <typename Kernel,
	          typename SourceIter, typename ChargeIter,
	          typename TargetIter, typename ResultIter>
	inline static void
  block_eval(const Kernel& K,
             SourceIter s_first, SourceIter s_last, ChargeIter c_first,
             TargetIter t_first, TargetIter t_last, ResultIter r_first)
  {
    typedef typename Kernel::source_type source_type;
    typedef typename Kernel::charge_type charge_type;
    typedef typename Kernel::target_type target_type;
    typedef typename Kernel::result_type result_type;

    static_assert(std::is_same<source_type,
                  typename std::iterator_traits<SourceIter>::value_type
                  >::value, "SourceIter::value_type != Kernel::source_type");
    static_assert(std::is_same<charge_type,
                  typename std::iterator_traits<ChargeIter>::value_type
                  >::value, "ChargeIter::value_type != Kernel::charge_type");
    static_assert(std::is_same<target_type,
                  typename std::iterator_traits<TargetIter>::value_type
                  >::value, "TargetIter::value_type != Kernel::target_type");
    static_assert(std::is_same<result_type,
                  typename std::iterator_traits<ResultIter>::value_type
                  >::value, "ResultIter::value_type != Kernel::result_type");

    for ( ; t_first != t_last; ++t_first, ++r_first) {
      const target_type& t = *t_first;
      result_type& r       = *r_first;

      SourceIter si = s_first;
      ChargeIter ci = c_first;
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
            typename SourceIter, typename ChargeIter, typename ResultIter>
  inline static void
  block_eval(const Kernel& K,
             SourceIter p1_first, SourceIter p1_last, ChargeIter c1_first,
             ResultIter r1_first,
             SourceIter p2_first, SourceIter p2_last, ChargeIter c2_first,
             ResultIter r2_first)
  {
    typedef typename Kernel::source_type source_type;
    typedef typename Kernel::charge_type charge_type;
    typedef typename Kernel::target_type target_type;
    typedef typename Kernel::result_type result_type;

    typedef typename Kernel::kernel_value_type kernel_value_type;

    static_assert(std::is_same<source_type,
                  typename std::iterator_traits<SourceIter>::value_type
                  >::value, "SourceIter::value_type != Kernel::source_type");
    static_assert(std::is_same<charge_type,
                  typename std::iterator_traits<ChargeIter>::value_type
                  >::value, "ChargeIter::value_type != Kernel::charge_type");
    static_assert(std::is_same<target_type,
                  typename std::iterator_traits<SourceIter>::value_type
                  >::value, "SourceIter::value_type != Kernel::target_type");
    static_assert(std::is_same<result_type,
                  typename std::iterator_traits<ResultIter>::value_type
                  >::value, "ResultIter::value_type != Kernel::result_type");

    for ( ; p1_first != p1_last; ++p1_first, ++c1_first, ++r1_first) {
      const source_type& pi = *p1_first;
      const charge_type& ci = *c1_first;
      result_type& ri       = *r1_first;

      SourceIter p2i = p2_first;
      ChargeIter c2i = c2_first;
      ResultIter r2i = r2_first;
      for ( ; p2i != p2_last; ++p2i, ++c2i, ++r2i) {
        const source_type& pj = *p2i;
        const charge_type& cj = *c2i;
        result_type& rj       = *r2i;

        kernel_value_type kij = K(pi,pj);
        ri += kij * cj;
        kernel_value_type kji = K.transpose(kij);
        rj += kji * ci;
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
    typedef typename Kernel::source_type source_type;
    typedef typename Kernel::charge_type charge_type;
    typedef typename Kernel::target_type target_type;
    typedef typename Kernel::result_type result_type;

    typedef typename Kernel::kernel_value_type kernel_value_type;

    static_assert(std::is_same<source_type, target_type>::value,
                  "source_type != target_type in symmetric P2P");
    static_assert(std::is_same<source_type,
                               typename SourceIter::value_type>::value,
                  "SourceIter::value_type != Kernel::source_type");
    static_assert(std::is_same<charge_type,
                               typename ChargeIter::value_type>::value,
                  "ChargeIter::value_type != Kernel::charge_type");
    static_assert(std::is_same<result_type,
                               typename ResultIter::value_type>::value,
                  "ResultIter::value_type != Kernel::result_type");

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

        kernel_value_type kij = K(pi,pj);
        ri += kij * cj;
        kernel_value_type kji = K.transpose(kij);
        rj += kji * ci;
      }
    }
  }
};

