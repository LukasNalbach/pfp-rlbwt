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

#ifndef PFP_RLBWT_PFP_RLBWT_OUTPUT_BLOCKS_TPP
#define PFP_RLBWT_PFP_RLBWT_OUTPUT_BLOCKS_TPP

#include <pfp_rlbwt/output/blocks.hpp>

namespace pfp_rlbwt {

// ############################# run_block #############################

inline void run_block::encode(bit_buffer& bits, const uint64_t* starts, const uint64_t* symbols, uint64_t m,
    uint64_t period, const symbol_code& code)
{
    const uint64_t span = starts[m];
    const uint8_t low = low_bits_of(m, span);
    // the offsets of the sampled codes relative to the first code; the last one is the largest
    uint64_t offset = 0, last_sample = 0;

    for (uint64_t k = 0; k < m; k++) {
        if (k % period == 0) last_sample = offset;
        offset += code.code_length(symbols[k]);
    }

    const uint8_t sample_bits = std::bit_width(last_sample);
    bits.put(sample_bits, width_bits);
    uint64_t prev_high = 0;

    // the upper bits: (high part - previous high part) 0-bits, then a 1-bit, for each start
    for (uint64_t k = 0; k < m; k++) {
        const uint64_t high = starts[k] >> low;
        bits.put_zeros(high - prev_high);
        bits.put(1, 1);
        prev_high = high;
    }

    bits.put_zeros(((span - 1) >> low) - prev_high);
    for (uint64_t k = 0; k < m; k++) bits.put(starts[k] & low_mask(low), low);
    offset = 0;

    for (uint64_t k = 0; k < m; k++) {
        if (k % period == 0 && k > 0) bits.put(offset, sample_bits);
        offset += code.code_length(symbols[k]);
    }

    for (uint64_t k = 0; k < m; k++) code.encode(bits, symbols[k]);
    bits.pad_to_byte();
}

inline uint64_t run_block::start(uint64_t k) const
{
    // the position of the k-th set upper bit
    uint64_t pos = 0, rest = k;

    while (true) {
        const uint64_t word = read_bits(data, upper_begin + pos, max_packed_bits);
        const uint64_t ones = std::popcount(word);
        if (rest < ones) break;
        rest -= ones;
        pos += max_packed_bits;
    }

    pos += select_in_word(read_bits(data, upper_begin + pos, max_packed_bits), rest);
    return ((pos - k) << low_bits) | read_bits(data, lows_begin + k * low_bits, low_bits);
}

inline uint64_t run_block::symbol(uint64_t k) const
{
    const uint64_t s = k / period;
    uint64_t bit = codes_begin + (s == 0 ? 0 : read_bits(data, samples_begin + (s - 1) * sample_bits, sample_bits));
    for (uint64_t j = s * period; j < k; j++) code->skip(data, bit);
    return code->decode(data, bit);
}

inline uint64_t run_block::predecessor(uint64_t x) const
{
    const uint64_t h = x >> low_bits;
    uint64_t pos = 0; // position in the upper bits of the first start with high part h (if any)

    if (h > 0) {
        // the (h-1)-th 0-bit ends the starts with high parts < h
        uint64_t rest = h - 1;

        while (true) {
            const uint64_t word = ~read_bits(data, upper_begin + pos, max_packed_bits) & low_mask(max_packed_bits);
            const uint64_t zeros = std::popcount(word);
            if (rest < zeros) break;
            rest -= zeros;
            pos += max_packed_bits;
        }

        pos += select_in_word(~read_bits(data, upper_begin + pos, max_packed_bits) & low_mask(max_packed_bits), rest) + 1;
    }

    // the starts with high part h are at the consecutive 1-bits from pos on
    uint64_t k = pos - h;
    const uint64_t x_low = x & low_mask(low_bits);

    while (k < m && pos < upper_bits && read_bits(data, upper_begin + pos, 1) != 0
        && read_bits(data, lows_begin + k * low_bits, low_bits) <= x_low) {
        k++;
        pos++;
    }

    return k - 1;
}

// ############################# lcp_block #############################

inline void lcp_block::encode(bit_buffer& bits, const uint64_t* values, uint64_t m, uint8_t position_bits)
{
    const auto [min_it, max_it] = std::minmax_element(values, values + m);
    const uint64_t base = *min_it;
    const uint8_t total = std::bit_width(*max_it - base);
    /* [0..total] differ[j] = number of k in [1..m-1] such that the bit width of (values[k] - base) xor
       (values[k-1] - base) is j: the high parts of the values k-1 and k differ for each split < j */
    std::vector<uint64_t> differ(total + 1, 0);
    for (uint64_t k = 1; k < m; k++) differ[std::bit_width((values[k] - base) ^ (values[k - 1] - base))]++;

    // the split with the fewest bits; the split total stores no high parts
    uint8_t split = total;
    uint64_t best = m * total;
    uint64_t runs = 1; // number of runs of equal high parts for the split in the loop

    for (uint8_t b = total; b-- > 0;) {
        runs += differ[b + 1];
        const uint64_t size = m * b + m + runs * (total - b);

        if (size < best) {
            best = size;
            split = b;
        }
    }

    const uint8_t high = total - split;
    bits.put(total, width_bits);
    bits.put(split, width_bits);
    bits.put(base, position_bits);
    for (uint64_t k = 0; k < m; k++) bits.put((values[k] - base) & low_mask(split), split);

    if (high > 0) {
        for (uint64_t k = 0; k < m; k++) {
            bits.put(k == 0 || ((values[k] - base) >> split) != ((values[k - 1] - base) >> split), 1);
        }

        for (uint64_t k = 0; k < m; k++) {
            const uint64_t part = (values[k] - base) >> split;
            if (k == 0 || part != ((values[k - 1] - base) >> split)) bits.put(part, high);
        }
    }

    bits.pad_to_byte();
}

inline uint64_t lcp_block::value(uint64_t k) const
{
    const uint64_t low = read_bits(data, lows_begin + k * split, split);
    if (high == 0) return base + low;
    // the number of runs of equal high parts that start at the values 0..k
    uint64_t runs = 0, pos = 0;

    for (; pos + max_packed_bits <= k + 1; pos += max_packed_bits) {
        runs += std::popcount(read_bits(data, starts_begin + pos, max_packed_bits));
    }

    runs += std::popcount(read_bits(data, starts_begin + pos, uint8_t(k + 1 - pos)));
    return base + ((read_bits(data, highs_begin + (runs - 1) * high, high) << split) | low);
}

} // namespace pfp_rlbwt

#endif // PFP_RLBWT_PFP_RLBWT_OUTPUT_BLOCKS_TPP
