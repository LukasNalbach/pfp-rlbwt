/**
 * part of LukasNalbach/pfp-rlbwt
 *
 * MIT License
 *
 * Copyright (c) Lukas Nalbach
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef PFP_RLBWT_PFP_RLBWT_PFP_RLBWT_HPP
#define PFP_RLBWT_PFP_RLBWT_PFP_RLBWT_HPP

#include <cstdint>
#include <string>

#include <pfp_rlbwt/text_ranges.hpp>
#include <pfp_rlbwt/types.hpp>
#include <rlbwt/text_sources.hpp>

// the largest length and alphabet size of the inputs of libsais that are suffix-sorted with 32-bit and with 40-bit
// indices; larger inputs are suffix-sorted with 48-bit indices (the tests lower them, so that they also cover the
// wider indices)
#ifndef PFP_RLBWT_MAX_INDEX32
#define PFP_RLBWT_MAX_INDEX32 INT32_MAX
#endif
#ifndef PFP_RLBWT_MAX_INDEX40
#define PFP_RLBWT_MAX_INDEX40 SA_INT40_MAX_N
#endif

namespace pfp_rlbwt {

/**
 * @brief builds the BWT of T$ with the prefix-free parsing of T and writes it run-length encoded to the runs file
 *        (see ext_runs) in the format prm.format, optionally together with the SA values at the starts and at the ends
 *        of the runs and the LCP values at the starts of the runs, and optionally the BWT with one symbol per position,
 *        the suffix array and the LCP array
 * @param text the text T (all symbols >= min_symbol)
 * @param prefix path prefix of the output files (see ext_runs, ext_bwt, ext_sa and ext_lcp)
 * @param prm parameters
 * @return sizes of the parsing and of the output
 */
template <typename sym_t>
result build(const text_source<sym_t>& text, const std::string& prefix, const params& prm);

/**
 * @brief builds the output of the text T that text_ranges reads (see the build above); if the length of T is not
 *        known before T is read, it is checked against the parameters after the parse
 * @param text the text T in ranges (all symbols >= min_symbol, after a renaming by text_ranges::ranges_read)
 * @param prefix path prefix of the output files
 * @param prm parameters
 * @return sizes of the parsing and of the output
 */
template <typename sym_t>
result build(text_ranges<sym_t>& text, const std::string& prefix, const params& prm);

} // namespace pfp_rlbwt

// the instantiations are compiled into the library (see src/instantiations.cpp)
#define PFP_RLBWT_EXTERN(sym_t) \
    extern template class rlbwt::file_source<sym_t>; \
    extern template class rlbwt::memory_source<sym_t>; \
    extern template pfp_rlbwt::result pfp_rlbwt::build<sym_t>(const rlbwt::text_source<sym_t>&, const std::string&, \
        const pfp_rlbwt::params&); \
    extern template pfp_rlbwt::result pfp_rlbwt::build<sym_t>(pfp_rlbwt::text_ranges<sym_t>&, const std::string&, \
        const pfp_rlbwt::params&);

STRING_UTILS_FOR_EACH_SYMBOL_TYPE(PFP_RLBWT_EXTERN)

#undef PFP_RLBWT_EXTERN

#include <pfp_rlbwt/pfp_rlbwt.tpp>

#endif // PFP_RLBWT_PFP_RLBWT_PFP_RLBWT_HPP
