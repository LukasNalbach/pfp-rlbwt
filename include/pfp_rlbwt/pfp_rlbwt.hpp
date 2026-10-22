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

#include <array>
#include <cstdint>
#include <memory>
#include <string>

#include <pfp_rlbwt/types.hpp>

// the largest length and alphabet size of the inputs of libsais that are suffix-sorted with 32-bit and with 40-bit
// indices (the tests lower them, so that they also cover the wider indices)
#ifndef PFP_RLBWT_MAX_INDEX32
#define PFP_RLBWT_MAX_INDEX32 INT32_MAX
#endif
#ifndef PFP_RLBWT_MAX_INDEX40
#define PFP_RLBWT_MAX_INDEX40 SA_INT40_MAX_N
#endif

namespace pfp_rlbwt {

// the readers of a file_source and of a memory_source (see pfp_rlbwt.tpp)
template <typename sym_t>
class file_reader;
template <typename sym_t>
class memory_reader;

/**
 * @brief a text stored in a file (sizeof(sym_t) bytes per symbol, little endian), optionally read in reverse and,
 *        for one-byte symbols, with its bytes mapped
 */
template <typename sym_t = uint8_t>
class file_source : public text_source<sym_t> {
private:
    std::string path; // path of the file
    uint64_t n = 0; // n, the number of symbols in the file
    bool reverse = false; // true <=> the text is the reversed file
    bool mapped = false; // true <=> map is applied to each byte
    // the map of the bytes: [0..255] map[c] = text symbol of byte c (only used if mapped)
    std::array<uint8_t, 256> map {};

public:
    /**
     * @brief opens a text file
     * @param path path of the file
     * @param reverse true <=> the text is the reversed file
     * @param map if not nullptr, the text symbol of byte c is (*map)[c] (only for one-byte symbols)
     */
    file_source(const std::string& path, bool reverse = false, const std::array<uint8_t, 256>* map = nullptr);

    uint64_t size() const override { return n; }

    std::unique_ptr<typename text_source<sym_t>::reader> new_reader() const override
    {
        return std::make_unique<file_reader<sym_t>>(path, n, reverse, mapped ? &map : nullptr);
    }
};

/**
 * @brief a text stored in memory, optionally read in reverse
 */
template <typename sym_t = uint8_t>
class memory_source : public text_source<sym_t> {
private:
    const sym_t* text; // text[0..n-1]
    uint64_t n; // n, the number of symbols
    bool reverse; // true <=> the text is text[0..n-1] reversed

public:
    /**
     * @brief a text in memory
     * @param text the symbols (not copied; they must stay valid until the construction has finished)
     * @param n number of symbols
     * @param reverse true <=> the text is text[0..n-1] reversed
     */
    memory_source(const sym_t* text, uint64_t n, bool reverse = false) : text(text), n(n), reverse(reverse) { }

    uint64_t size() const override { return n; }

    std::unique_ptr<typename text_source<sym_t>::reader> new_reader() const override
    {
        return std::make_unique<memory_reader<sym_t>>(text, n, reverse);
    }
};

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

// the instantiations are compiled into the library (see src/instantiations.cpp)
#define PFP_RLBWT_EXTERN(sym_t) \
    extern template class file_source<sym_t>; \
    extern template class memory_source<sym_t>; \
    extern template result build<sym_t>(const text_source<sym_t>&, const std::string&, const params&);

PFP_RLBWT_EXTERN(uint8_t)
PFP_RLBWT_EXTERN(uint16_t)
PFP_RLBWT_EXTERN(sym24_t)
PFP_RLBWT_EXTERN(uint32_t)
PFP_RLBWT_EXTERN(sym40_t)
PFP_RLBWT_EXTERN(sym48_t)

#undef PFP_RLBWT_EXTERN

} // namespace pfp_rlbwt

#include <pfp_rlbwt/pfp_rlbwt.tpp>

#endif // PFP_RLBWT_PFP_RLBWT_PFP_RLBWT_HPP
