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
 * @file blocks.hpp
 * @brief the blocks of the compressed format: a block of runs (Elias-Fano coded starts and Huffman-coded symbols, see
 *        section_run_blocks) and a block of LCP values (RLE-GEF, see section_lcp_blocks); each block starts at a
 *        byte and is decoded from its bytes alone (together with numbers in RAM)
 */

#ifndef PFP_RLBWT_PFP_RLBWT_OUTPUT_BLOCKS_HPP
#define PFP_RLBWT_PFP_RLBWT_OUTPUT_BLOCKS_HPP

#include <algorithm>
#include <bit>
#include <cstdint>
#include <vector>

#include <pfp_rlbwt/output/bits.hpp>
#include <pfp_rlbwt/output/symbol_code.hpp>

namespace pfp_rlbwt {

/**
 * @brief a block of m runs of the compressed format; the starts of its runs relative to the start of its first run are
 *        x_0 = 0 < x_1 < ... < x_{m-1} < span, where span is the distance from its first start to the first start of the
 *        next block (or to n). Layout (bit-packed):
 *        1. sample_bits (width_bits bits): bits per symbol sample;
 *        2. the upper bits of the starts: bit (x_k >> low_bits) + k is set for each k ((span - 1) >> low_bits) + m
 *           bits, low_bits = floor(log2(span / m)) (Elias-Fano);
 *        3. the low_bits low bits of each start;
 *        4. the symbol samples: for s in [1..(m - 1) / period] the offset (in bits) of the code of run s * period
 *           relative to the first code (sample_bits bits each);
 *        5. the codes of the symbols (see symbol_code).
 */
class run_block {
public:
    static constexpr uint64_t upper_begin = width_bits; // offset (in bits) of the upper bits

private:
    const uint8_t* data = nullptr; // the bytes of the block
    const symbol_code* code = nullptr; // the code of the symbols
    uint64_t m = 0; // number of runs
    uint64_t period = 0; // number of runs per symbol sample
    uint8_t low_bits = 0; // number of low bits per start
    uint8_t sample_bits = 0; // bits per symbol sample
    uint64_t upper_bits = 0; // number of upper bits
    uint64_t lows_begin = 0, samples_begin = 0, codes_begin = 0; // offsets (in bits) of the parts

public:
    // returns the number of low bits per start of a block of m runs that span span BWT positions
    static uint8_t low_bits_of(uint64_t m, uint64_t span) { return span > m ? std::bit_width(span / m) - 1 : 0; }

    /**
     * @brief appends a block to bits, padded to a multiple of 8 bits
     * @param bits the bits the block is appended to (it must end at a byte)
     * @param starts [0..m] the starts of the runs of the block relative to the first one, then the span
     * @param symbols [0..m-1] the symbols of the runs
     * @param m number of runs (>= 1)
     * @param period number of runs per symbol sample
     * @param code the code of the symbols
     */
    static void encode(bit_buffer& bits, const uint64_t* starts, const uint64_t* symbols, uint64_t m, uint64_t period,
        const symbol_code& code);

    run_block() = default;

    /**
     * @brief a block stored at data
     * @param data the bytes of the block, followed by read_padding readable bytes
     * @param m number of runs
     * @param span the span of the block
     * @param period number of runs per symbol sample
     * @param code the code of the symbols
     */
    run_block(const uint8_t* data, uint64_t m, uint64_t span, uint64_t period, const symbol_code& code)
        : data(data), code(&code), m(m), period(period), low_bits(low_bits_of(m, span)),
          sample_bits(uint8_t(read_bits(data, 0, width_bits))), upper_bits(((span - 1) >> low_bits) + m),
          lows_begin(upper_begin + upper_bits), samples_begin(lows_begin + m * low_bits),
          codes_begin(samples_begin + ((m - 1) / period) * sample_bits)
    { }

    // returns the start of run k relative to the first start
    uint64_t start(uint64_t k) const;

    // returns the symbol of run k
    uint64_t symbol(uint64_t k) const;

    // returns the largest k such that start(k) <= x (x < span)
    uint64_t predecessor(uint64_t x) const;

    // calls fnc(start, symbol) for the runs in order (starts relative to the first one)
    template <typename fnc_t>
    void for_each(fnc_t fnc) const
    {
        uint64_t pos = 0; // position in the upper bits after the last set bit found
        uint64_t bit = codes_begin; // offset of the next code

        for (uint64_t k = 0; k < m; k++) {
            uint64_t word = read_bits(data, upper_begin + pos, max_packed_bits);

            while (word == 0) {
                pos += max_packed_bits;
                word = read_bits(data, upper_begin + pos, max_packed_bits);
            }

            pos += std::countr_zero(word);
            const uint64_t low = read_bits(data, lows_begin + k * low_bits, low_bits);
            fnc(((pos - k) << low_bits) | low, code->decode(data, bit));
            pos++;
        }
    }
};

/**
 * @brief a block of m LCP values of the compressed format, coded with RLE-GEF (the run-length variant of Generalized
 *        Elias-Fano): the values minus their minimum base are split into high parts (their bits from bit split on) and
 *        low parts (their split low bits); the high part is stored once for each maximal run of equal high parts.
 *        Layout (bit-packed):
 *        1. total (width_bits bits): the bit width of the largest value minus base;
 *        2. split (width_bits bits): the number of low bits; high = total - split;
 *        3. base (position_bits bits);
 *        4. the low parts (split bits each);
 *        5. if high > 0: m bits, bit k set <=> value k starts a run of equal high parts;
 *        6. if high > 0: the high part of each run (high bits each).
 *        split minimizes the size of the block.
 */
class lcp_block {
private:
    const uint8_t* data = nullptr; // the bytes of the block
    uint64_t m = 0; // number of values
    uint8_t split = 0, high = 0; // number of low and high bits
    uint64_t base = 0; // the smallest value
    uint64_t lows_begin = 0, starts_begin = 0, highs_begin = 0; // offsets (in bits) of the parts

public:
    /**
     * @brief appends a block to bits, padded to a multiple of 8 bits
     * @param bits the bits the block is appended to (it must end at a byte)
     * @param values [0..m-1] the values
     * @param m number of values (>= 1)
     * @param position_bits bits of base
     */
    static void encode(bit_buffer& bits, const uint64_t* values, uint64_t m, uint8_t position_bits);

    lcp_block() = default;

    /**
     * @brief a block stored at data
     * @param data the bytes of the block, followed by read_padding readable bytes
     * @param m number of values
     * @param position_bits bits of base
     */
    lcp_block(const uint8_t* data, uint64_t m, uint8_t position_bits)
        : data(data), m(m), split(uint8_t(read_bits(data, width_bits, width_bits))),
          high(uint8_t(read_bits(data, 0, width_bits) - split)),
          base(read_bits(data, 2 * width_bits, position_bits)), lows_begin(2 * width_bits + position_bits),
          starts_begin(lows_begin + m * split), highs_begin(starts_begin + (high > 0 ? m : 0))
    { }

    // returns value k
    uint64_t value(uint64_t k) const;

    // calls fnc(value) for the values in order
    template <typename fnc_t>
    void for_each(fnc_t fnc) const
    {
        uint64_t run = 0; // number of runs of equal high parts started so far
        uint64_t high_part = 0; // the high part of the current run

        for (uint64_t k = 0; k < m; k++) {
            if (high > 0 && read_bits(data, starts_begin + k, 1) != 0) {
                high_part = read_bits(data, highs_begin + run * high, high) << split;
                run++;
            }

            fnc(base + (high_part | read_bits(data, lows_begin + k * split, split)));
        }
    }
};

} // namespace pfp_rlbwt

#include <pfp_rlbwt/output/blocks.tpp>

#endif // PFP_RLBWT_PFP_RLBWT_OUTPUT_BLOCKS_HPP
