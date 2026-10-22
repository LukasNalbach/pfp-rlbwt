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

/**
 * @file run_model.hpp
 * @brief the context models of the symbols and the lengths of the runs of the BWT
 */

#ifndef PFP_RLBWT_PFP_RLBWT_COMP_RUN_MODEL_HPP
#define PFP_RLBWT_PFP_RLBWT_COMP_RUN_MODEL_HPP

#include <algorithm>
#include <bit>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include <coding/byte_codes.hpp>
#include <coding/integer_model.hpp>
#include <coding/mixer.hpp>
#include <coding/range_coder.hpp>
#include <misc/bits.hpp>

namespace pfp_rlbwt::comp {

/**
 * @brief codes the runs of the BWT of T$ in the order of the BWT: the symbol c in [0..sigma-1] of a run into one
 *        stream, its length into another one. The symbol is coded with three mixed models with the contexts: none,
 *        the symbol of the previous run, the symbols of the previous two runs. The length l is coded as its bit
 *        width b with three mixed models with the contexts: none, c, (the bit width of the previous run length, the
 *        bit width of the last run length of c); then the modeled_length_bits bits of l below its leading one with
 *        the context b, then its remaining bits with probability 1/2. The run of $ is not coded (see skip_terminator).
 */
class run_model {
private:
    uint64_t sigma = 0; // the number of distinct symbols of the text; sigma is the context value of $
    uint8_t max_width = 0; // the bit width of the largest possible run length n
    mixed_trees<3> symbols; // the models of the symbols
    mixed_trees<3> widths; // the models of the bit widths minus one of the run lengths
    bit_trees top_bits; // the models of the bits below the leading one of a run length
    uint64_t prev = 0; // the symbol of the previous run (sigma: $ or none)
    uint64_t prev2 = 0; // the symbol of the run before the previous run (sigma: $ or none)
    uint64_t prev_width = 0; // the bit width of the length of the previous run (0: none)
    std::vector<uint8_t> last_width; // [0..sigma-1] the bit width of the length of the last run of c (0: none)

    /**
     * @brief codes a run
     * @param symbol_coder the coder of the symbols
     * @param length_coder the coder of the lengths
     * @param c the symbol of the run (set by a decoder)
     * @param length the length of the run (set by a decoder)
     */
    template <typename coder_t>
    void code(coder_t& symbol_coder, coder_t& length_coder, uint64_t& c, uint64_t& length)
    {
        const uint64_t num_widths = uint64_t(max_width) + 1;
        c = symbols.code(symbol_coder, c, { 0, prev, prev2 * (sigma + 1) + prev });
        if (c >= sigma) throw invalid_file();
        const uint64_t width = widths.code(length_coder, std::bit_width(length) - 1,
            { 0, c, prev_width * num_widths + last_width[c] }) + 1;
        if (width > max_width) throw invalid_file();

        const uint8_t below = uint8_t(width - 1); // the number of bits below the leading one
        const uint8_t modeled = std::min(below, modeled_length_bits);
        const uint8_t raw = below - modeled;
        const uint64_t top = top_bits.code(length_coder, (length >> raw) & low_bits_mask(modeled), modeled, width);
        const uint64_t rest = length_coder.code_bits(length & low_bits_mask(raw), raw);
        length = (uint64_t(1) << below) | (top << raw) | rest;

        prev2 = prev;
        prev = c;
        prev_width = width;
        last_width[c] = uint8_t(width);
    }

public:
    /**
     * @brief creates the models
     * @param sigma the number of distinct symbols of the text
     * @param n the length of the text
     */
    run_model(uint64_t sigma, uint64_t n)
        : sigma(sigma), max_width(std::bit_width(n)), prev(sigma), prev2(sigma), last_width(sigma, 0)
    {
        const uint64_t num_widths = uint64_t(max_width) + 1;
        symbols = mixed_trees<3>(std::bit_width(sigma - 1), { 1, sigma + 1, (sigma + 1) * (sigma + 1) });
        widths = mixed_trees<3>(std::bit_width(uint64_t(max_width) - 1), { 1, sigma, num_widths * num_widths });
        top_bits = bit_trees(modeled_length_bits, num_widths);
    }

    /**
     * @brief encodes a run
     * @param symbol_coder the encoder of the symbols
     * @param length_coder the encoder of the lengths
     * @param c the symbol of the run, in [0..sigma-1]
     * @param length the length of the run, in [1..n]
     */
    void encode(range_encoder& symbol_coder, range_encoder& length_coder, uint64_t c, uint64_t length)
    {
        code(symbol_coder, length_coder, c, length);
    }

    /**
     * @brief decodes a run; throws std::runtime_error if the symbol or the bit width of the length is out of range
     * @param symbol_coder the decoder of the symbols
     * @param length_coder the decoder of the lengths
     * @param c is set to the symbol of the run
     * @param length is set to the length of the run
     */
    void decode(range_decoder& symbol_coder, range_decoder& length_coder, uint64_t& c, uint64_t& length)
    {
        code(symbol_coder, length_coder, c, length);
    }

    // updates the contexts after the run of $, which is not coded (its length is 1)
    void skip_terminator()
    {
        prev2 = prev;
        prev = sigma;
        prev_width = 1;
    }
};

} // namespace pfp_rlbwt::comp

#endif // PFP_RLBWT_PFP_RLBWT_COMP_RUN_MODEL_HPP
