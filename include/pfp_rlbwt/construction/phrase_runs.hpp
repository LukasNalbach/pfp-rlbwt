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

#ifndef PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_PHRASE_RUNS_HPP
#define PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_PHRASE_RUNS_HPP

#include <algorithm>
#include <bit>
#include <cstdint>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>

#include <pfp_rlbwt/construction/common.hpp>

namespace pfp_rlbwt {

/* The suffixes of P that start with a word s: a phrase run i of s with m_i phrases is followed by the suffix Y_i,
   which starts with a word other than s. The suffix starting at the phrase of run i that has j phrases of the run
   left (j in [1..m_i]) is s^j Y_i. If Y_i < s (down), s^j Y_i < s^j' Y_i' for all j' > j, and if Y_i > s (up),
   s^j Y_i > s^j' Y_i' for all j' > j. Hence the suffixes starting with s are ordered: first those of the down runs
   by j in increasing order, then those of the up runs by j in decreasing order; suffixes with the same j are ordered
   by Y_i. Layer j of a class (down or up) is the set of runs with m_i >= j; layer 1 contains all runs of the
   class, and the suffix with j = m_i is the head of run i (the suffix starting at its first phrase). */

/**
 * @brief a set of integers in [0..u-1], as a bit vector with the number of integers in the set before each 64-bit
 *        word of it
 */
class ranked_bits {
private:
    static constexpr uint64_t word_bits = std::numeric_limits<uint64_t>::digits;

    // bit s % 64 of bits[s / 64] = 1 <=> s is in the set
    std::vector<uint64_t> bits;
    // [0..u/64] before[b] = number of integers < 64b in the set
    std::vector<uint64_t> before;

public:
    ranked_bits() = default;

    // the set of the values in [0..u-1] for which value(z) returns true, for z in [0..num_values-1]
    template <typename value_t>
    ranked_bits(value_t value, uint64_t num_values, uint64_t u)
        : bits(u / word_bits + 1, 0), before(u / word_bits + 1, 0)
    {
        for (uint64_t z = 0; z < num_values; z++) {
            const uint64_t s = value(z);
            bits[s / word_bits] |= uint64_t { 1 } << (s % word_bits);
        }

        for (uint64_t b = 0; b + 1 < bits.size(); b++) before[b + 1] = before[b] + std::popcount(bits[b]);
    }

    // returns whether s is in the set
    inline bool contains(uint64_t s) const { return (bits[s / word_bits] >> (s % word_bits)) & 1; }

    // returns the number of integers < s in the set (s <= u)
    inline uint64_t rank(uint64_t s) const
    {
        return before[s / word_bits] + std::popcount(bits[s / word_bits] & ((uint64_t { 1 } << (s % word_bits)) - 1));
    }
};

/**
 * @brief a sequence A[0..t-1] of distinct values in [0..u-1] that answers, for a prefix A[0..len-1], the k-th
 *        smallest value and the number of values smaller than x (a wavelet matrix)
 */
class prefix_wavelet {
private:
    static constexpr uint64_t word_bits = std::numeric_limits<uint64_t>::digits;

    uint64_t t = 0; // the number of values
    uint8_t levels = 0; // the number of bits per value
    // [0..levels-1] bits[l] = bit l (from the most significant one) of the values, in the order of level l
    std::vector<std::vector<uint64_t>> bits;
    // [0..levels-1][0..t/64] ones_before[l][b] = number of 1-bits in bits[l][0..b-1]
    std::vector<std::vector<uint64_t>> ones_before;
    // [0..levels-1] zeros[l] = number of 0-bits on level l
    std::vector<uint64_t> zeros;

    // returns the number of 1-bits in bits[l] before position i
    inline uint64_t rank1(uint8_t l, uint64_t i) const
    {
        const uint64_t b = i / word_bits;
        const uint64_t r = i % word_bits;
        return ones_before[l][b] + (r == 0 ? 0 : std::popcount(bits[l][b] & (~uint64_t { 0 } >> (word_bits - r))));
    }

public:
    prefix_wavelet() = default;

    // stores the values A[0..t-1] = values, each smaller than u
    prefix_wavelet(std::vector<uint64_t> values, uint64_t u) : t(values.size())
    {
        if (t == 0) return;
        levels = std::max<uint8_t>(1, std::bit_width(u - 1));
        bits.resize(levels);
        ones_before.resize(levels);
        zeros.resize(levels);
        std::vector<uint64_t> next(t);

        for (uint8_t l = 0; l < levels; l++) {
            const uint8_t shift = levels - 1 - l;
            bits[l].assign(t / word_bits + 1, 0);
            ones_before[l].assign(t / word_bits + 1, 0);
            uint64_t z = 0;

            for (uint64_t i = 0; i < t; i++) {
                if ((values[i] >> shift) & 1) bits[l][i / word_bits] |= uint64_t { 1 } << (i % word_bits);
                else z++;
            }

            for (uint64_t b = 0; b + 1 < bits[l].size(); b++) {
                ones_before[l][b + 1] = ones_before[l][b] + std::popcount(bits[l][b]);
            }

            zeros[l] = z;
            uint64_t zi = 0, oi = z;
            for (uint64_t i = 0; i < t; i++) next[((values[i] >> shift) & 1) ? oi++ : zi++] = values[i];
            values.swap(next);
        }
    }

    // returns the k-th smallest value (from 0) of A[0..len-1]
    uint64_t kth_smallest(uint64_t len, uint64_t k) const
    {
        uint64_t lo = 0, hi = len, val = 0;

        for (uint8_t l = 0; l < levels; l++) {
            const uint64_t ones_lo = rank1(l, lo), ones_hi = rank1(l, hi);
            const uint64_t num_zeros = (hi - ones_hi) - (lo - ones_lo);

            if (k < num_zeros) {
                lo -= ones_lo;
                hi -= ones_hi;
            } else {
                k -= num_zeros;
                lo = zeros[l] + ones_lo;
                hi = zeros[l] + ones_hi;
                val |= uint64_t { 1 } << (levels - 1 - l);
            }
        }

        return val;
    }

    // returns the number of values in A[0..len-1] that are smaller than x
    uint64_t count_smaller(uint64_t len, uint64_t x) const
    {
        uint64_t lo = 0, hi = len, cnt = 0;

        for (uint8_t l = 0; l < levels && hi > lo; l++) {
            const uint64_t ones_lo = rank1(l, lo), ones_hi = rank1(l, hi);

            if ((x >> (levels - 1 - l)) & 1) {
                cnt += (hi - ones_hi) - (lo - ones_lo);
                lo = zeros[l] + ones_lo;
                hi = zeros[l] + ones_hi;
            } else {
                lo -= ones_lo;
                hi -= ones_hi;
            }
        }

        return cnt;
    }
};

// the phrase runs of one class (down or up) of a word, and the layers of their suffixes
class run_class {
private:
    // the distinct numbers of phrases of the runs in increasing order, each with the number of runs that have it
    std::vector<std::pair<uint64_t, uint64_t>> copies_count;
    // [0..copies_count.size()] copies_prefix[y] = sum of m * count over copies_count[0..y-1]
    std::vector<uint64_t> copies_prefix;
    // [0..copies_count.size()] runs_prefix[y] = sum of count over copies_count[0..y-1]
    std::vector<uint64_t> runs_prefix;
    // the indices x (in runs) of the runs with at least two phrases, in decreasing order of their numbers of phrases
    prefix_wavelet long_runs;

public:
    /* the runs in increasing order of the suffix after them (layer 1), each as its offset in the range of the runs of
       the word in the suffix array of the parse of the phrase runs */
    uint_array runs;
    uint64_t num_runs = 0; // the number of runs
    uint64_t total = 0; // the number of suffixes of the class (the sum of the numbers of phrases)
    uint64_t max_copies = 0; // the largest number of phrases of a run

    run_class() = default;

    /**
     * @brief builds the layers
     * @param runs the runs in increasing order of the suffix after them (see runs)
     * @param num_runs their number
     * @param copies_count the distinct numbers of phrases of the runs in increasing order with their numbers of runs
     * @param long_idx the indices (in runs) of the runs with at least two phrases, in decreasing order of their
     *        numbers of phrases
     */
    run_class(uint_array runs, uint64_t num_runs, std::vector<std::pair<uint64_t, uint64_t>> copies_count,
        std::vector<uint64_t> long_idx)
        : copies_count(std::move(copies_count)), runs(std::move(runs)), num_runs(num_runs)
    {
        const uint64_t c = this->copies_count.size();
        copies_prefix.assign(c + 1, 0);
        runs_prefix.assign(c + 1, 0);

        for (uint64_t y = 0; y < c; y++) {
            const auto [m, cnt] = this->copies_count[y];
            copies_prefix[y + 1] = copies_prefix[y] + m * cnt;
            runs_prefix[y + 1] = runs_prefix[y] + cnt;
        }

        total = copies_prefix[c];
        max_copies = c == 0 ? 0 : this->copies_count.back().first;
        long_runs = prefix_wavelet(std::move(long_idx), num_runs);
    }

    // returns the number of suffixes in the layers 1..j (the sum of min(m, j) over the runs)
    inline uint64_t layers_size(uint64_t j) const
    {
        // the number of distinct numbers of phrases <= j
        const uint64_t y = std::partition_point(copies_count.begin(), copies_count.end(),
            [&](const std::pair<uint64_t, uint64_t>& mc) { return mc.first <= j; }) - copies_count.begin();
        return copies_prefix[y] + j * (num_runs - runs_prefix[y]);
    }

    // returns the layer j (in [1..max_copies]) with layers_size(j-1) <= o < layers_size(j)
    inline uint64_t layer_of(uint64_t o) const
    {
        uint64_t lo = 1, hi = max_copies;

        while (lo < hi) {
            const uint64_t mid = (lo + hi) / 2;
            if (layers_size(mid) > o) hi = mid;
            else lo = mid + 1;
        }

        return lo;
    }

    // returns the number of runs with at least j >= 2 phrases
    inline uint64_t num_long(uint64_t j) const
    {
        const uint64_t y = std::partition_point(copies_count.begin(), copies_count.end(),
            [&](const std::pair<uint64_t, uint64_t>& mc) { return mc.first < j; }) - copies_count.begin();
        return num_runs - runs_prefix[y];
    }

    // returns the index x (in runs) of the idx-th run (from 0) of layer j
    inline uint64_t run_in_layer(uint64_t j, uint64_t idx) const
    {
        return j == 1 ? idx : long_runs.kth_smallest(num_long(j), idx);
    }

    // returns the position of runs[x] in layer j <= the number of phrases of runs[x]
    inline uint64_t index_in_layer(uint64_t x, uint64_t j) const
    {
        return j == 1 ? x : long_runs.count_smaller(num_long(j), x);
    }
};

/**
 * @brief the suffixes of P that start with a word s with a phrase run of at least two phrases: their offsets
 *        q in [0..occ-1] in the range of the suffixes starting with s are the down class (layers 1, 2, ...), then
 *        the up class (layers ..., 2, 1). The heads, in increasing order of their offsets, are the runs of s in the
 *        order of the suffix array of the parse of the phrase runs.
 */
struct repeated_word_ranks {
    uint64_t word = 0; // the word s
    uint64_t occ = 0; // the number of phrases equal to s
    uint64_t advance = 0; // the advance of s
    uint64_t first_run = 0; // the rank of the first run of s in the suffix array of the parse of the phrase runs
    uint64_t first_rank = 0; // the rank of the first suffix of P that starts with s among the suffixes of P
    run_class down, up;
    /* the offsets of the heads in increasing order: [0..R-1] head_offset[b] = the offset of the head of the b-th run of
       s in the order of the suffix array of the parse of the phrase runs */
    uint_array head_offset;
    uint64_t num_runs = 0; // R, the number of runs of s

    // returns the offset of the suffix in layer j of the run with index x in the down (is_down) or the up class
    inline uint64_t offset_of(bool is_down, uint64_t x, uint64_t j) const
    {
        if (is_down) return down.layers_size(j - 1) + down.index_in_layer(x, j);
        return down.total + up.total - up.layers_size(j) + up.index_in_layer(x, j);
    }

    // returns the class (true: down), the layer j and the index x (in the runs of the class) of the suffix of offset q
    inline std::tuple<bool, uint64_t, uint64_t> locate(uint64_t q) const
    {
        if (q < down.total) {
            const uint64_t j = down.layer_of(q);
            return { true, j, down.run_in_layer(j, q - down.layers_size(j - 1)) };
        }

        const uint64_t o = q - down.total;
        const uint64_t j = up.layer_of(up.total - 1 - o);
        return { false, j, up.run_in_layer(j, o - (up.total - up.layers_size(j))) };
    }

    /* returns, for the suffix of offset q, the offset of its phrase run in the range of the runs of s in the suffix
       array of the parse of the phrase runs, and the number j of phrases of the run that the suffix starts with */
    inline std::pair<uint64_t, uint64_t> run_of(uint64_t q) const
    {
        const auto [is_down, j, x] = locate(q);
        return { (is_down ? down : up).runs[x], j };
    }

    // returns the offset of the t-th (from 0) suffix whose offset is not the offset of a head
    inline uint64_t non_head_offset(uint64_t t) const
    {
        // the number of heads b with head_offset[b] - b <= t (the number of the other suffixes before head b)
        uint64_t lo = 0, hi = num_runs;

        while (lo < hi) {
            const uint64_t mid = (lo + hi) / 2;
            if (head_offset[mid] - mid <= t) lo = mid + 1;
            else hi = mid;
        }

        return t + lo;
    }
};

} // namespace pfp_rlbwt

#endif // PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_PHRASE_RUNS_HPP
