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

#ifndef PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_SUFFIX_ARRAYS_HPP
#define PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_SUFFIX_ARRAYS_HPP

#include <algorithm>
#include <bit>
#include <chrono>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
#include <misc/files.hpp>
#include <misc/log.hpp>
#include <misc/sort.hpp>

#include <libsais.h>
#include <libsais40.h>
#include <move_data_structure/move_data_structure.hpp>
#include <omp.h>

#include <pfp_rlbwt/construction/parse.hpp>
#include <pfp_rlbwt/construction/phrase_runs.hpp>
#include <pfp_rlbwt/types.hpp>

namespace pfp_rlbwt {

/**
 * @brief the first positions of the runs of a sequence over the positions [0..u-1], as a bit vector over [0..u] that
 *        stores the number of set bits before each cache line in the line (see cache_line_rank_bit_vector):
 *        is_start(y) and run_of(y) read one cache line
 */
class run_start_bits {
private:
    // bit y = 1 <=> a run starts at position y (y < u); bit u = 1
    cache_line_rank_bit_vector bits;

public:
    run_start_bits() = default;

    // allocates the bits of the positions [0..u] and marks position u
    run_start_bits(uint64_t u, uint16_t threads)
    {
        bits.allocate(u + 1, threads);
        mark(u);
    }

    // marks position y as the first position of a run (several threads can mark positions at once)
    inline void mark(uint64_t y) { bits.set_parallel(y); }

    // stores the number of marks before each cache line; called after all positions are marked
    void finish()
    {
        constexpr uint64_t line_bits = cache_line_rank_bit_vector::bits_per_line;
        uint64_t marks = 0;

        for (uint64_t l = 0; l <= bits.size() / line_bits; l++) {
            bits.set_line_rank(l, marks);
            const uint64_t* line = bits.line(l * line_bits);
            for (uint64_t j = 1; j < cache_line_rank_bit_vector::words_per_line; j++) marks += std::popcount(line[j]);
        }
    }

    // returns whether position y is the first position of a run
    inline bool is_start(uint64_t y) const { return bits.bit(y); }

    // returns the run containing position y
    inline uint64_t run_of(uint64_t y) const { return bits.rank(y + 1) - 1; }

    // fetches the cache line that run_of(y) reads
    inline void prefetch(uint64_t y) const { __builtin_prefetch(bits.line(y + 1)); }

    // returns the first position after position y < u that starts a run, or u if there is none
    inline uint64_t next_start(uint64_t y) const
    {
        constexpr uint64_t line_bits = cache_line_rank_bit_vector::bits_per_line;
        constexpr uint64_t word_bits = std::numeric_limits<uint64_t>::digits;
        /* bit i is bit i % word_bits of the word 1 + (i % line_bits) / word_bits of its line (line_bits is a multiple
           of word_bits) */
        uint64_t i = y + 1;
        uint64_t word = bits.line(i)[1 + (i % line_bits) / word_bits] >> (i % word_bits);

        while (word == 0) {
            i += word_bits - i % word_bits;
            word = bits.line(i)[1 + (i % line_bits) / word_bits];
        }

        return i + std::countr_zero(word);
    }
};

/**
 * @brief Psi of the suffix array of the parse of the phrase runs (see parse_sa::psi) as a move data structure. Each run
 *        [a, a + len) of the BWT of the parse is the image of an interval [y, y + len) of ranks, Psi[y + i] = a + i
 *        (see for_each_run_by_value); these intervals, split by the balancing of move_data_structure, are the input
 *        intervals. Row x stores p(x), q(x) = Psi[p(x)], the input interval idx(x) containing q(x), the word that the
 *        suffixes of the ranks of interval x start with, and the advance of their first phrase runs; the advance of a
 *        phrase run of m phrases of the word s is its length minus w, m (|s| - w). The word and the advance are 0 for
 *        the interval [0, 1) of the end.
 */
class psi_mds {
private:
    // number of bits of idx(x) and of the advance in a row (both are smaller than 2^value_bits)
    static constexpr uint64_t value_bits = 8 * max_value_bytes;
    static constexpr uint64_t value_mask = low_bits_mask(value_bits);
    // number of bits of the word that are stored above idx(x) and above the advance
    static constexpr uint64_t word_part_bits = std::numeric_limits<uint64_t>::digits - value_bits;

    /* a row, aligned to its size (it lies in one cache line): idx_word stores idx(x) and above it the low
       word_part_bits bits of the word, advance_word stores the advance and above it the high bits of the word */
    struct alignas(4 * sizeof(uint64_t)) row {
        uint64_t p;
        uint64_t q;
        uint64_t idx_word;
        uint64_t advance_word;
    };

    // the rows: [0..k] rows[x] = the row of input interval x; rows[k].p = J + 1, the number of ranks
    std::vector<row, cache_aligned_no_init_allocator<row>> rows;
    uint64_t k = 0; // the number of input intervals
    // marks the first ranks p(0), ..., p(k-1) of the input intervals
    run_start_bits starts;

public:
    psi_mds() = default;

    /**
     * @brief builds Psi from the intervals before the balancing
     * @param pairs [0..r-1] the first rank y_j of the j-th interval and Psi[y_j], by increasing y_j; freed
     * @param num_ranks J + 1, the number of ranks
     * @param words [0..r-1] entry j = the word of interval j (< 2^(2 word_part_bits))
     * @param advances [0..r-1] entry j = the advance of the phrase runs of interval j
     * @param threads number of threads
     */
    psi_mds(mds_pairs&& pairs, uint64_t num_ranks, const interleaved_bit_aligned_vectors& words,
        const interleaved_bit_aligned_vectors& advances, uint16_t threads);

    // returns the number of input intervals
    inline uint64_t num_intervals() const { return k; }

    // returns the input interval containing rank y
    inline uint64_t interval_of(uint64_t y) const { return starts.run_of(y); }

    // returns the first rank of input interval x
    inline uint64_t p(uint64_t x) const { return rows[x].p; }

    // returns Psi of the first rank of input interval x
    inline uint64_t q(uint64_t x) const { return rows[x].q; }

    // returns Psi[y] for a rank y of input interval x
    inline uint64_t at(uint64_t y, uint64_t x) const { return q(x) + (y - p(x)); }

    // returns Psi[y]
    inline uint64_t operator[](uint64_t y) const { return at(y, interval_of(y)); }

    // returns the word that the suffixes of the ranks of input interval x start with
    inline uint64_t word(uint64_t x) const
    {
        return (rows[x].idx_word >> value_bits) | ((rows[x].advance_word >> value_bits) << word_part_bits);
    }

    // returns the advance of the first phrase runs of the suffixes of the ranks of input interval x
    inline uint64_t advance(uint64_t x) const { return rows[x].advance_word & value_mask; }

    // changes the rank y of input interval x to Psi[y] and x to the input interval containing Psi[y]
    inline void move(uint64_t& y, uint64_t& x) const
    {
        y = at(y, x);
        x = rows[x].idx_word & value_mask;
        while (p(x + 1) <= y) x++;
    }

    // fetches the row idx(x), which move reads after row x
    inline void prefetch_move(uint64_t x) const { __builtin_prefetch(&rows[rows[x].idx_word & value_mask]); }

    // fetches the cache line that interval_of(y) reads
    inline void prefetch_interval_of(uint64_t y) const { starts.prefetch(y); }

    // fetches row x
    inline void prefetch_row(uint64_t x) const { __builtin_prefetch(&rows[x]); }
};

/* the ends of the phrase runs by rank, sampled (see sort_parse_recursively), with SA[y] = the phrase run starting the
   suffix of rank y of the parse of the phrase runs, F_i = the end of phrase run i, F_{-1} = w - 1 and F_J = F_{J-1}:
   F_{SA[y]} and F_{SA[y]-1} at the first and at the last position y of each run of the BWT of the parse (see
   parse_sa::bwt_runs), and F_{SA[y]} at the ranks y of the runs of the repeated words. At most run boundaries of the
   BWT of T$, the position has a rank y (see bwt_pos::rank) of a run of a repeated word, or y or psi[y] is the first or
   the last position of a run of the BWT of the parse; from the other ranks, end_at walks with Psi. */
struct sampled_ends {
    uint64_t num_runs = 0; // r, the number of runs of the BWT of the parse of the phrase runs (0: not sampled)
    /* the ends at the first and at the last position of each run: [0..4r-1] for the first position y of run k,
       values[4k] = F_{SA[y]} and values[4k+1] = F_{SA[y]-1}; for its last position y, values[4k+2] = F_{SA[y]} and
       values[4k+3] = F_{SA[y]-1} */
    uint_array values;
    // the ranges [beg, end) of the ranks of the runs of the repeated words, in increasing order
    std::vector<std::pair<uint64_t, uint64_t>> repeated_ranges;
    // [0..repeated_ranges.size()] repeated_offset[i] = the number of ranks in the ranges before range i
    std::vector<uint64_t> repeated_offset;
    /* the ends at the ranks of the ranges: repeated_ends[repeated_offset[i] + y - beg] = F_{SA[y]} for a rank y of the
       range i = [beg, end) */
    uint_array repeated_ends;

    // the result of the lookups below if F_{SA[y]} is not sampled
    static constexpr uint64_t no_end = std::numeric_limits<uint64_t>::max();

    // returns whether the ends are not sampled
    bool empty() const { return num_runs == 0; }

    /* returns F_{SA[y]} (offset 0) or F_{SA[y]-1} (offset 1) if y is the first or the last position of a run of the
       BWT of the parse (marked in bwt_runs), else no_end */
    inline uint64_t boundary_end(const run_start_bits& bwt_runs, uint64_t y, uint64_t offset) const
    {
        if (bwt_runs.is_start(y)) return values[4 * bwt_runs.run_of(y) + offset];
        if (bwt_runs.is_start(y + 1)) return values[4 * bwt_runs.run_of(y) + 2 + offset];
        return no_end;
    }

    // returns F_{SA[y]} if y is a rank of a run of a repeated word, else no_end
    inline uint64_t repeated_end(uint64_t y) const
    {
        const auto it = std::upper_bound(repeated_ranges.begin(), repeated_ranges.end(), std::make_pair(y, no_end));
        if (it == repeated_ranges.begin()) return no_end;
        const uint64_t i = it - repeated_ranges.begin() - 1;
        if (y >= repeated_ranges[i].second) return no_end;
        return repeated_ends[repeated_offset[i] + y - repeated_ranges[i].first];
    }

    /* returns F_{SA[y]} if y or psi_y = psi[y] is the first or the last position of a run of the BWT of the parse or y
       is a rank of a run of a repeated word, else no_end */
    inline uint64_t sampled_end(const run_start_bits& bwt_runs, uint64_t y, uint64_t psi_y) const
    {
        uint64_t end = boundary_end(bwt_runs, y, 0);
        if (end != no_end) return end;
        // F_{SA[y]} = F_{SA[psi[y]]-1}
        end = boundary_end(bwt_runs, psi_y, 1);
        if (end != no_end) return end;
        return repeated_end(y);
    }

    /* returns F_{SA[y]} (psi_y = psi[y]); if it is not sampled, F_{SA[z]} minus the advances (see psi_mds) of the
       phrase runs SA[y] + 1, ..., SA[z] for the first rank z = psi^i[y] (i >= 1) whose end is sampled */
    inline uint64_t end_at(const run_start_bits& bwt_runs, uint64_t y, uint64_t psi_y, const psi_mds& psi) const
    {
        uint64_t end = sampled_end(bwt_runs, y, psi_y);
        if (end != no_end) [[likely]] return end;
        /* z = psi^i[y] (i >= 1) lies in the input interval x; it is not the first or the last position of a run
           (sampled_end and the last step checked it), so z != 0 */
        uint64_t z = psi_y, x = psi.interval_of(z);
        uint64_t advances = 0; // the sum of the advances of the phrase runs SA[y] + 1, ..., SA[z] - 1

        while (true) {
            advances += psi.advance(x);
            const uint64_t prev = z;
            psi.move(z, x);
            end = boundary_end(bwt_runs, z, 1);
            if (end == no_end) end = repeated_end(prev);
            if (end != no_end) return end - advances;
        }
    }
};

/* the suffix array SA of the parse of the phrase runs (SA[x] = the phrase run starting the suffix of rank x, SA[0] =
   J), as Psi with the starts C of the ranges of ranks of the suffixes per first word, the runs of its BWT and, if they
   are sampled, the ends of the phrase runs by rank; for the words with a phrase run of two or more phrases, also the
   order of the suffixes of the parse of the phrases that start with them. The suffixes of the parse of the phrases
   that start at the first phrase (the head) of a phrase run are ordered as the suffixes of P, the parse of the phrase
   runs (see phrase_runs.hpp). */
template <typename idx_t>
struct parse_sa {
    // Psi: [0..J] psi[x] = ISA[SA[x]+1], the rank of the suffix one phrase run later (psi[0] = ISA[0], see psi_mds)
    psi_mds psi;
    /* the starts of the ranges of ranks per first word: [0..D+1] C[0] = 0, C[1] = 1 (rank 0 is the suffix of the end)
       and C[s+2] = C[s+1] + the number of phrase runs of word s, so the suffixes of the ranks [C[s+1], C[s+2]) start
       with word s; C[D+1] = J+1 */
    uint_array C;
    // [0..D+1] C_intervals[c] = the input interval of Psi that starts at rank C[c] (see psi_mds)
    uint_array C_intervals;
    uint64_t num_words = 0; // D, the number of words
    // the words with a phrase run of two phrases or more (repeated words), in increasing order
    std::vector<repeated_word_ranks> repeated;
    ranked_bits repeated_set; // the set of the repeated words
    // the set of the ranks of the first runs of the repeated words (see repeated_word_ranks::first_run)
    ranked_bits first_runs;
    /* [0..repeated.size()] extra_before[y] = number of phrases minus number of phrase runs of the words
       repeated[0..y-1] */
    std::vector<uint64_t> extra_before;
    // the ends of the phrase runs by rank, if they are sampled (else phrase_arrays::ends_by_rank stores them)
    sampled_ends ends;
    /* marks the first positions of the r runs of the BWT of the parse of the phrase runs, whose entries are values of
       the parse (0 for the end, v >= 1 for the letter v - 1, see phrase_arrays::P) */
    run_start_bits bwt_runs;
    // the values of the runs: [0..r-1] bwt_run_values[k] = the value of run k (see bwt_merger::symbol_at_rank)
    uint_array bwt_run_values;

    // returns the ranks of the suffixes of the parse of the phrases that start with word s (nullptr: s is not repeated)
    inline const repeated_word_ranks* find_repeated(uint64_t s) const
    {
        if (repeated.empty() || !repeated_set.contains(s)) return nullptr;
        return &repeated[repeated_set.rank(s)];
    }

    // returns the rank among the suffixes of the parse of the phrases of the suffix starting at the head of the phrase
    // run of rank x
    inline uint64_t rank_among_phrases(uint64_t x) const
    {
        // the number of repeated words whose runs start at a rank <= x
        const uint64_t y = first_runs.rank(x + 1);

        if (y > 0) {
            const repeated_word_ranks& rw = repeated[y - 1];
            if (x < rw.first_run + rw.num_runs) return rw.first_rank + rw.head_offset[x - rw.first_run];
        }

        return x + extra_before[y];
    }

    /**
     * @brief returns the rank among the suffixes of the parse of the phrases of the suffix one phrase after the suffix
     *        with offset q among those that start with word s
     * @param s the word
     * @param rw find_repeated(s)
     * @param q [0..num_occ[s]-1] the offset
     */
    uint64_t next_rank_among_phrases(uint64_t s, const repeated_word_ranks* rw, uint64_t q) const
    {
        if (rw == nullptr) return rank_among_phrases(psi[C[s + 1] + q]);
        // the suffixes in layer 1 are followed by the head of the run after their run, the others by non-heads of s
        const uint64_t up_layer_1 = rw->occ - rw->up.num_runs;
        if (q < rw->down.num_runs) return rank_among_phrases(psi[rw->first_run + rw->down.runs[q]]);
        if (q >= up_layer_1) return rank_among_phrases(psi[rw->first_run + rw->up.runs[q - up_layer_1]]);
        return rw->first_rank + rw->non_head_offset(q - rw->down.num_runs);
    }
};

/**
 * @brief the generalized suffix array of the dictionary (see sa_file), the positions where a suffix equals its
 *        predecessor and the word containing each position
 */
template <typename idx_t>
struct dictionary_sa {
    /* the generalized suffix array of dict.text (in a file or in RAM, see temporary_array): [0..size-1] SA[i] =
       starting position of the i-th smallest suffix (sizeof(idx_t) bytes per value) */
    temporary_array sa_file;
    /* bits [0..size-1] in size/64+1 words, bit i is bit i % 64 of equals_previous_bits[i / 64]: bit i (i >= 1) = 1
       <=> dict.text[SA[i]..] and dict.text[SA[i-1]..] are equal up to their separators and longer than w; bit 0 = 0 */
    malloc_array_t<uint64_t> equals_previous_bits;
    // bits [0..size-1] in size/64+1 words, stored as equals_previous_bits: bit x = 1 <=> a word starts at position x
    malloc_array_t<uint64_t> word_starts;
    // [0..size/64] word_starts_before[b] = number of word starts in dict.text[0..64b-1]
    uint_array word_starts_before;

    // returns whether the suffix of rank i equals the suffix of rank i-1 (see equals_previous_bits)
    inline bool equals_previous(uint64_t i) const { return (equals_previous_bits[i / 64] >> (i % 64)) & 1; }

    // returns the word containing position x of dict.text
    inline uint64_t word_of(uint64_t x) const
    {
        return word_starts_before[x / 64] + std::popcount(word_starts[x / 64] & (~uint64_t { 0 } >> (63 - x % 64))) - 1;
    }

    // fetches the cache lines that word_of(x) reads
    inline void prefetch_word_of(uint64_t x) const
    {
        __builtin_prefetch(&word_starts[x / 64]);
        word_starts_before.prefetch(x / 64);
    }

    /**
     * @brief reads SA[beg..end-1] from sa_file
     * @param file a reader of sa_file (with sizeof(idx_t) bytes per value)
     * @param beg first rank
     * @param end rank after the last one
     * @param values receives SA[beg..end-1]
     */
    static void read_sa(array_reader& file, uint64_t beg, uint64_t end, no_init_vector<idx_t>& values)
    {
        resize_no_init(values, end - beg);
        file.read(beg, end, values.data());
    }
};

/**
 * @brief returns the number of distinct values the suffix sorting of the dictionary uses (one per separator and
 *        one per symbol, see build_dictionary_sa)
 */
template <typename sym_t>
inline uint64_t dictionary_alphabet_size(const dictionary<sym_t>& dict, uint64_t max_symbol)
{
    if constexpr (sizeof(sym_t) == 1) return 256;
    else return dict.num_words + (max_symbol < dict.size ? max_symbol + 1 : dict.size);
}

// entries of extra space (libsais: fs) per alphabet value that hold all bucket arrays of libsais
constexpr uint64_t free_space_per_value = 6;

// extra space (libsais: fs) at the end of the suffix array that libsais uses for a byte alphabet
constexpr uint64_t byte_alphabet_free_space = free_space_per_value * 256;

/**
 * @brief computes the suffix array of the integer array T[0..n-1] with values in [0..k-1]
 * @param fs number of entries of extra space after SA[n-1] that libsais may use
 */
template <typename idx_t>
void build_integer_sa(idx_t* T, idx_t* SA, uint64_t n, uint64_t k, uint64_t fs, uint16_t threads)
{
    int64_t ret;
    if constexpr (std::is_same_v<idx_t, int32_t>)
        ret = libsais_int_omp(T, SA, n, k, std::min<uint64_t>(fs, std::numeric_limits<int32_t>::max() - n), threads);
    else if constexpr (std::is_same_v<idx_t, sa_int40_t>)
        ret = libsais40_impl::libsais40_long_omp(T, SA, n, k, fs, threads);
    else ret = libsais48_impl::libsais48_long_omp(T, SA, n, k, fs, threads);
    if (ret != 0) throw std::runtime_error("pfp_rlbwt: libsais failed on the parsing");
}

// computes the generalized suffix array of the string set T[0..n-1] (separated by 0)
template <typename idx_t>
void build_generalized_sa(const uint8_t* T, idx_t* SA, uint64_t n, uint16_t threads)
{
    int64_t ret;
    if constexpr (std::is_same_v<idx_t, int32_t>)
        ret = libsais_gsa_omp(T, SA, n, byte_alphabet_free_space, nullptr, threads);
    else if constexpr (std::is_same_v<idx_t, sa_int40_t>)
        ret = libsais40_impl::libsais40_gsa_omp(T, SA, n, byte_alphabet_free_space, nullptr, threads);
    else ret = libsais48_impl::libsais48_gsa_omp(T, SA, n, byte_alphabet_free_space, nullptr, threads);
    if (ret != 0) throw std::runtime_error("pfp_rlbwt: libsais failed on the dictionary");
}

/**
 * @brief computes the order of the suffixes of the parse of the phrases that start with a repeated word (see
 *        repeated_word_ranks)
 * @param psa the suffix array of the parse of the phrase runs, Psi and C
 * @param arr the phrase runs
 * @param J number of phrase runs
 * @param w length of a trigger string
 * @param dict the dictionary
 * @param threads number of threads
 */
template <typename idx_t, typename sym_t>
void build_repeated_word_ranks(parse_sa<idx_t>& psa, phrase_arrays<idx_t, sym_t>& arr, uint64_t J, uint64_t w,
    const dictionary<sym_t>& dict, uint16_t threads)
{
    const uint64_t R = arr.repeated_words.size();
    psa.repeated.resize(R);
    psa.repeated_set = std::move(arr.repeated_set);
    psa.first_runs = ranked_bits([&](uint64_t y) { return uint64_t(psa.C[arr.repeated_words[y].word + 1]); }, R, J + 2);
    psa.extra_before.assign(R + 1, 0);
    uint64_t repeated_runs = 0; // the number of phrase runs of the repeated words

    for (uint64_t y = 0; y < R; y++) {
        const uint64_t s = arr.repeated_words[y].word;
        psa.extra_before[y + 1] = psa.extra_before[y] + arr.repeated_words[y].extra_phrases;
        repeated_runs += psa.C[s + 2] - psa.C[s + 1];
    }

    // computes the ranks of the suffixes starting with word repeated[y] with p threads
    auto build_word = [&](uint64_t y, uint64_t p) {
        const uint64_t s = arr.repeated_words[y].word;
        repeated_word_ranks& rw = psa.repeated[y];
        rw.word = s;
        rw.occ = dict.num_occ[s];
        rw.advance = dict.length(s) - w;
        const uint64_t first = psa.C[s + 1], last = psa.C[s + 2];
        rw.first_run = first;
        rw.first_rank = first + psa.extra_before[y];
        rw.num_runs = last - first;

        // the number of phrases of the run starting the suffix of rank x >= 1
        auto copies = [&](uint64_t x) { return psa.psi.advance(psa.psi.interval_of(x)) / rw.advance; };

        // the runs followed by a smaller suffix (the down runs) come first: [first..mid-1]
        uint64_t mid = first;

        for (uint64_t hi = last; mid < hi;) {
            const uint64_t x = (mid + hi) / 2;
            if (psa.psi[x] < first) mid = x + 1;
            else hi = x;
        }

        for (bool is_down : { true, false }) {
            const uint64_t beg = is_down ? first : mid, end = is_down ? mid : last;
            const uint64_t num_runs = end - beg;
            // (rank of the suffix after the run, x) of the runs, sorted by the rank after them
            std::vector<std::pair<uint64_t, uint64_t>> order(num_runs);

            #pragma omp parallel for num_threads(p) if (p > 1)
            for (uint64_t x = beg; x < end; x++) order[x - beg] = { psa.psi[x], x };

            radix_sort(order.begin(), order.end(), key_of_first {}, p);
            uint_array runs(num_runs, value_bytes(last - first));

            #pragma omp parallel for num_threads(p) if (p > 1)
            for (uint64_t z = 0; z < num_runs; z++) runs.set(z, order[z].second - first);

            // returns the position in runs of the run starting the suffix of rank x
            auto index = [&](uint64_t x) {
                return std::lower_bound(order.begin(), order.end(), std::make_pair(uint64_t(psa.psi[x]), x))
                    - order.begin();
            };

            // the groups of runs with equal numbers of phrases, in the order of SA: (copies, (begin, end))
            std::vector<std::pair<uint64_t, std::pair<uint64_t, uint64_t>>> groups;

            for (uint64_t x = beg; x < end;) {
                const uint64_t m = copies(x);
                uint64_t e = x + 1;
                while (e < end && copies(e) == m) e++;
                groups.push_back({ m, { x, e } });
                x = e;
            }

            // the down groups have increasing, the up groups decreasing numbers of phrases
            if (!is_down) std::reverse(groups.begin(), groups.end());
            std::vector<std::pair<uint64_t, uint64_t>> copies_count;
            std::vector<uint64_t> long_idx;

            for (const auto& [m, range] : groups) copies_count.emplace_back(m, range.second - range.first);

            for (auto it = groups.rbegin(); it != groups.rend() && it->first >= 2; it++) {
                const uint64_t size = long_idx.size();
                long_idx.resize(size + it->second.second - it->second.first);

                #pragma omp parallel for num_threads(p) if (p > 1)
                for (uint64_t x = it->second.first; x < it->second.second; x++) {
                    long_idx[size + x - it->second.first] = index(x);
                }
            }

            order = std::vector<std::pair<uint64_t, uint64_t>>();

            (is_down ? rw.down : rw.up) = run_class(std::move(runs), num_runs, std::move(copies_count),
                std::move(long_idx));
        }

        rw.head_offset = uint_array(rw.num_runs, value_bytes(rw.occ));

        for (bool is_down : { true, false }) {
            const run_class& cls = is_down ? rw.down : rw.up;

            #pragma omp parallel for num_threads(p) if (p > 1)
            for (uint64_t z = 0; z < cls.num_runs; z++) {
                const uint64_t x = cls.runs[z];
                rw.head_offset.set(x, rw.offset_of(is_down, z, copies(first + x)));
            }
        }
    };

    // a word with at least a share of 1 / threads of the runs is processed with all threads, the others with one each
    std::vector<uint64_t> small;

    for (uint64_t y = 0; y < R; y++) {
        const uint64_t s = arr.repeated_words[y].word;
        if ((psa.C[s + 2] - psa.C[s + 1]) * threads >= repeated_runs) build_word(y, threads);
        else small.push_back(y);
    }

    #pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
    for (uint64_t z = 0; z < small.size(); z++) build_word(small[z], 1);

    arr.repeated_words = std::vector<repeated_word>();
}

/**
 * @brief computes the starts C of the ranges of ranks of the suffixes of the parse of the phrase runs per first word
 *        (see parse_sa::C) and the number of words
 * @param psa receives C and the number of words
 * @param arr the phrase runs (the repeated words)
 * @param J number of phrase runs
 * @param dict the dictionary
 */
template <typename idx_t, typename sym_t>
void build_word_ranges(parse_sa<idx_t>& psa, const phrase_arrays<idx_t, sym_t>& arr, uint64_t J,
    const dictionary<sym_t>& dict)
{
    const uint64_t D = dict.num_words;
    psa.num_words = D;
    psa.C = uint_array(D + 2, value_bytes(J + 1));
    psa.C.set(0, 0);
    psa.C.set(1, 1);

    if (arr.repeated_words.empty()) {
        for (uint64_t s = 0; s < D; s++) psa.C.set(s + 2, psa.C[s + 1] + dict.num_occ[s]);
    } else {
        // a repeated word has fewer phrase runs than phrases
        auto rw_it = arr.repeated_words.begin();

        for (uint64_t s = 0; s < D; s++) {
            uint64_t runs = dict.num_occ[s];

            if (rw_it != arr.repeated_words.end() && rw_it->word == s) {
                runs -= rw_it->extra_phrases;
                rw_it++;
            }

            psa.C.set(s + 2, psa.C[s + 1] + runs);
        }
    }
}

/**
 * @brief calls f(v, s, m, repeated) for each value v >= 1 of the parse of the phrase runs (the letter v - 1, see
 *        phrase_arrays::P), whose phrase runs have the word s and m phrases each; repeated <=> s is a repeated word
 *        (it has one value per key of its phrase runs, see run_key)
 * @param arr the phrase runs (the keys of the repeated words)
 * @param D the number of words
 * @param f the function
 */
template <typename idx_t, typename sym_t, typename value_fnc_t>
void for_each_value(const phrase_arrays<idx_t, sym_t>& arr, uint64_t D, value_fnc_t f)
{
    uint64_t v = 1;
    uint64_t y = 0; // index in arr.repeated_words of the next repeated word

    // the letters are ordered by word
    for (uint64_t s = 0; s < D; s++) {
        if (y < arr.repeated_words.size() && arr.repeated_words[y].word == s) {
            for (uint64_t key : arr.repeated_words[y].keys) f(v++, s, run_copies(key), true);
            y++;
        } else {
            f(v++, s, 1, false);
        }
    }
}

/**
 * @brief returns the BWT symbols of the values of the parse of the phrase runs: [0..num_letters] symbols[v] = the
 *        symbol before the last w symbols of the word of value v, which precedes the phrase after a phrase run of
 *        value v; symbols[0] = arr.last_symbol
 * @param arr the phrase runs (the keys of the repeated words)
 * @param w length of a trigger string
 * @param dict the dictionary
 */
template <typename idx_t, typename sym_t>
std::vector<sym_t> symbols_of_values(const phrase_arrays<idx_t, sym_t>& arr, uint64_t w,
    const dictionary<sym_t>& dict)
{
    std::vector<sym_t> symbols(arr.num_letters + 1);
    symbols[0] = arr.last_symbol;

    for_each_value(arr, dict.num_words, [&](uint64_t v, uint64_t s, uint64_t, bool) {
        symbols[v] = dict.text[dict.word_start[s] + dict.length(s) - w - 1];
    });

    return symbols;
}

// the rows word and advance of Psi (see psi_mds) per value of the parse of the phrase runs
struct value_rows {
    // the row of the words and the row of the advances of the intervals of Psi
    using rows_t = std::pair<interleaved_bit_aligned_vectors, interleaved_bit_aligned_vectors>;

    // [0..num_letters] words[v] = the word of the phrase runs of value v (0 for the end)
    std::vector<uint32_t> words;
    // [0..num_letters] advances[v] = the advance of the phrase runs of value v (0 for the end)
    std::vector<uint64_t> advances;
    uint8_t word_width = 0, advance_width = 0; // bits per word and per advance

    template <typename idx_t, typename sym_t>
    value_rows(const phrase_arrays<idx_t, sym_t>& arr, uint64_t w, const dictionary<sym_t>& dict)
        : words(arr.num_letters + 1, 0), advances(arr.num_letters + 1, 0)
    {
        uint64_t max_advance = 0;

        for_each_value(arr, dict.num_words, [&](uint64_t v, uint64_t s, uint64_t m, bool) {
            words[v] = uint32_t(s);
            advances[v] = m * (dict.length(s) - w);
            max_advance = std::max(max_advance, advances[v]);
        });

        word_width = uint8_t(std::bit_width(dict.num_words));
        advance_width = uint8_t(std::bit_width(max_advance));
    }

    // returns the rows of k intervals of Psi, not initialized (see set)
    rows_t new_rows(uint64_t k, uint16_t threads) const
    {
        interleaved_bit_aligned_vectors word_row({ word_width }), advance_row({ advance_width });
        word_row.resize_no_init(k, threads);
        advance_row.resize_no_init(k, threads);
        return { std::move(word_row), std::move(advance_row) };
    }

    // sets the word and the advance of interval j, whose value is v (several threads can set them at once)
    void set(rows_t& rows, uint64_t j, uint64_t v) const
    {
        rows.first.set_parallel<0, uint64_t>(j, words[v]);
        rows.second.set_parallel<0, uint64_t>(j, advances[v]);
    }
};

/**
 * @brief computes the input intervals of Psi that start at the ranks C (see parse_sa::C_intervals)
 * @param psa the suffix array of the parse of the phrase runs (C and Psi); receives C_intervals
 * @param threads number of threads
 */
template <typename idx_t>
void build_word_intervals(parse_sa<idx_t>& psa, uint16_t threads)
{
    const uint64_t num = psa.num_words + 2;
    psa.C_intervals = uint_array(num, value_bytes(psa.psi.num_intervals()));

    #pragma omp parallel for num_threads(threads)
    for (uint64_t c = 0; c < num; c++) psa.C_intervals.set(c, psa.psi.interval_of(psa.C[c]));
}

/**
 * @brief builds Psi as a move data structure (see psi_mds), the runs of the BWT of the parse of the phrase runs with
 *        their values (see parse_sa::bwt_runs) and the symbols of the values (see phrase_arrays::value_symbols) from
 *        Psi as an array: the intervals of Psi are the maximal ranges of the ranks of one value that Psi maps to
 *        consecutive ranks, i.e., to a run of the BWT
 * @param psa the suffix array of the parse of the phrase runs (C); receives Psi and the runs
 * @param arr the phrase runs (the letters of the repeated words); receives the symbols of the values
 * @param psi [0..J] Psi
 * @param J number of phrase runs
 * @param w length of a trigger string
 * @param dict the dictionary
 * @param threads number of threads
 */
template <typename idx_t, typename sym_t>
void build_psi_from_array(parse_sa<idx_t>& psa, phrase_arrays<idx_t, sym_t>& arr, const idx_t* psi, uint64_t J,
    uint64_t w, const dictionary<sym_t>& dict, uint16_t threads)
{
    const uint64_t num_values = arr.num_letters + 1;
    /* the first ranks of the values: [0..num_values] value_start[v] = the first rank of the suffixes that start with
       value v, value_start[num_values] = J + 1; the suffixes of the values of a repeated word are counted with the
       letters of its phrase runs */
    std::vector<uint64_t> value_start(num_values + 1, 0);
    value_start[1] = 1;

    for_each_value(arr, dict.num_words, [&](uint64_t v, uint64_t s, uint64_t, bool repeated) {
        if (!repeated) value_start[v + 1] = psa.C[s + 2] - psa.C[s + 1];
    });

    for (const auto& [i, letter] : arr.repeated_run_letters) value_start[letter + 2]++;
    for (uint64_t v = 0; v < num_values; v++) value_start[v + 1] += value_start[v];

    // part t of p parts processes the ranks [t * (J+1) / p, (t+1) * (J+1) / p)
    const uint64_t p = std::max<uint64_t>(1, std::min<uint64_t>(threads, J + 1));
    auto part_start = [&](uint64_t t) { return (t * (J + 1)) / p; };

    // calls f(y, v) for each rank y in part t that starts an interval, whose value is v
    auto for_interval_starts = [&](uint64_t t, auto f) {
        const uint64_t end = part_start(t + 1);
        uint64_t y = part_start(t);
        uint64_t v = std::upper_bound(value_start.begin(), value_start.end(), y) - value_start.begin() - 1;

        for (; y < end; y++) {
            while (value_start[v + 1] <= y) v++;
            if (y == value_start[v] || uint64_t(psi[y]) != uint64_t(psi[y - 1]) + 1) f(y, v);
        }
    };

    // [0..p] part_first[t] = the number of intervals of the parts before part t
    std::vector<uint64_t> part_first(p + 1, 0);

    #pragma omp parallel for num_threads(p)
    for (uint64_t t = 0; t < p; t++) for_interval_starts(t, [&](uint64_t, uint64_t) { part_first[t + 1]++; });

    for (uint64_t t = 0; t < p; t++) part_first[t + 1] += part_first[t];
    const uint64_t k = part_first[p];
    const value_rows by_value(arr, w, dict);
    mds_pairs pairs(k, J + 1);
    value_rows::rows_t rows = by_value.new_rows(k, threads);
    psa.bwt_runs = run_start_bits(J + 1, threads);

    #pragma omp parallel for num_threads(p)
    for (uint64_t t = 0; t < p; t++) {
        uint64_t j = part_first[t];

        for_interval_starts(t, [&](uint64_t y, uint64_t v) {
            pairs.set<0>(j, y);
            pairs.set<1>(j, psi[y]);
            by_value.set(rows, j, v);
            psa.bwt_runs.mark(psi[y]);
            j++;
        });
    }

    psa.bwt_runs.finish();
    psa.bwt_run_values = uint_array(k, value_bytes(num_values));

    #pragma omp parallel for num_threads(p)
    for (uint64_t t = 0; t < p; t++) {
        for_interval_starts(t, [&](uint64_t y, uint64_t v) { psa.bwt_run_values.set(psa.bwt_runs.run_of(psi[y]), v); });
    }

    psa.psi = psi_mds(std::move(pairs), J + 1, rows.first, rows.second, threads);
    build_word_intervals(psa, threads);
    arr.value_symbols = symbols_of_values(arr, w, dict);
}

/**
 * @brief computes Psi of the suffix array of the parse of the phrase runs and the runs of its BWT (see
 *        build_psi_from_array) and frees the parse; orders the ends of the phrase runs by rank (see
 *        phrase_arrays::ends_by_rank), reading them and the inverse suffix array, which is stored in a temporary array
 *        meanwhile (see temporary_array), and removes both; for the repeated words also computes the order of the
 *        suffixes of the parse of the phrases that start with them
 * @param arr the phrase runs
 * @param J number of phrase runs
 * @param w length of a trigger string
 * @param dict the dictionary
 * @param file path of the file that stores the inverse suffix array temporarily
 * @param prm parameters (threads and logging)
 * @param time phase timer (see log_phase_start)
 */
template <typename idx_t, typename sym_t>
parse_sa<idx_t> build_parse_sa(phrase_arrays<idx_t, sym_t>& arr, uint64_t J, uint64_t w,
    const dictionary<sym_t>& dict, const std::string& file, const params& prm,
    std::chrono::steady_clock::time_point& time)
{
    const uint16_t threads = prm.num_threads;
    parse_sa<idx_t> psa;
    log_phase_start(prm.log, time, "building SA of " + level_name("P"));
    const uint64_t free_space = free_space_per_value * (arr.num_letters + 1);
    malloc_array_t<idx_t> SA = new_malloc_array<idx_t>(J + 1 + free_space);
    build_integer_sa(arr.P.get(), SA.get(), J + 1, arr.num_letters + 1, free_space, threads);
    resize_malloc_array(SA, J + 1);

    build_word_ranges(psa, arr, J, dict);

    /* the inverse suffix array ISA, which is stored in the memory of arr.P (the parse is not needed any more), and Psi
       from it, which is stored in the memory of SA */
    malloc_array_t<idx_t> ISA = std::move(arr.P);

    #pragma omp parallel for num_threads(threads)
    for (uint64_t x = 0; x <= J; x++) ISA[uint64_t(SA[x])] = idx_t(x);

    #pragma omp parallel for num_threads(threads)
    for (uint64_t x = 0; x <= J; x++) SA[x] = ISA[uint64_t(SA[x]) == J ? 0 : uint64_t(SA[x]) + 1];

    log_phase_end(prm.log, time);
    temporary_array isa_file(file, temporary_arrays_in_memory((J + 1) * sizeof(idx_t)));

    if (isa_file.is_in_memory()) {
        isa_file.store(std::move(ISA), J + 1);
    } else {
        log_phase_start(prm.log, time, "storing data structures to disk");
        isa_file.store(std::move(ISA), J + 1);
        log_phase_end(prm.log, time);
    }

    log_phase_start(prm.log, time, "building M_Psi of " + level_name("P"));
    build_psi_from_array(psa, arr, SA.get(), J, w, dict, threads);
    SA.reset();
    log_phase_end(prm.log, time);

    /* ends_by_rank[ISA[v]] = F_{v-1} for v in [0..J]; thread t of p threads reads ISA[v] and F_{v-1} from their
       temporary arrays for v in [t * (J+1) / p, (t+1) * (J+1) / p), with readers that together read
       temporary_read_bytes bytes at once */
    log_phase_start(prm.log, time, "sorting the phrase ends by SA of " + level_name("P"));
    arr.ends_by_rank = uint_array(J + 1, arr.position_bytes);

    run_on_threads(threads, [&](uint64_t t, uint64_t p) {
        const uint64_t beg = (t * (J + 1)) / p, end = ((t + 1) * (J + 1)) / p;
        if (beg >= end) return;
        constexpr uint64_t num_files = 2;
        constexpr uint64_t read_bytes = temporary_read_bytes / num_files;
        array_reader isa(isa_file, sizeof(idx_t), read_bytes);
        array_reader ends(arr.ends_file, arr.position_bytes, read_bytes);
        for (uint64_t v = beg; v < end; v++) arr.ends_by_rank.set(isa[v], v == 0 ? w - 1 : ends[v - 1]);
    });

    isa_file.remove();
    arr.ends_file.remove();
    if (!arr.repeated_words.empty()) build_repeated_word_ranks(psa, arr, J, w, dict, threads);
    log_phase_end(prm.log, time);
    return psa;
}

/**
 * @brief computes the generalized suffix array of the dictionary, which is stored in a temporary array (see
 *        dictionary_sa::sa_file), and the positions where a suffix of length > w equals its predecessor; the equal
 *        suffixes are found with Phi as in the computation of the PLCP array (Phi is computed in the place of the
 *        suffix array if the suffix array is stored in a file)
 * @param dict the dictionary
 * @param max_symbol the largest symbol of the dictionary
 * @param w length of a trigger string
 * @param file path of the file that stores the suffix array (see dictionary_sa::sa_file)
 * @param prm parameters (threads and logging)
 * @param time phase timer (see log_phase_start)
 */
template <typename idx_t, typename sym_t>
dictionary_sa<idx_t> build_dictionary_sa(const dictionary<sym_t>& dict, uint64_t max_symbol, uint64_t w,
    const std::string& file, const params& prm, std::chrono::steady_clock::time_point& time)
{
    const uint16_t threads = prm.num_threads;
    const uint64_t size = dict.size;
    const uint64_t D = dict.num_words;
    const uint64_t num_bit_blocks = size / 64 + 1; // number of 64-bit words of a bit vector of size bits
    dictionary_sa<idx_t> dsa;
    const sym_t* text = dict.text.get();
    malloc_array_t<idx_t> SA;

    log_phase_start(prm.log, time, "building SA of " + level_name("D"));

    if constexpr (sizeof(sym_t) == 1) {
        SA = new_malloc_array<idx_t>(size + byte_alphabet_free_space);
        build_generalized_sa(text, SA.get(), size, threads);
    } else {
        /* the separator of word s is s (no two separators are equal) and symbol c is D + c, or D + the rank of c among
           the symbols of the dictionary if the symbols are larger than the dictionary */
        // the extra space of libsais: free_space_per_value entries per alphabet value, if they fit into size entries
        const uint64_t alphabet = dictionary_alphabet_size(dict, max_symbol);
        const uint64_t free_space = free_space_per_value * alphabet <= size ? free_space_per_value * alphabet : 0;
        SA = new_malloc_array<idx_t>(size + free_space);
        malloc_array_t<idx_t> T = new_malloc_array<idx_t>(size);
        // the distinct symbols of the dictionary in increasing order (only if max_symbol >= size)
        std::vector<sym_t> symbols;

        if (max_symbol >= size) {
            symbols.assign(text, text + size);
            radix_sort(symbols.begin(), symbols.end(), key_of_value {}, threads);
            symbols.erase(std::unique(symbols.begin(), symbols.end()), symbols.end());
        }

        #pragma omp parallel for num_threads(threads)
        for (uint64_t s = 0; s < D; s++) {
            for (uint64_t x = dict.word_start[s]; x + 1 < dict.word_start[s + 1]; x++) {
                const uint64_t c = symbols.empty() ? uint64_t(text[x])
                    : uint64_t(std::lower_bound(symbols.begin(), symbols.end(), text[x]) - symbols.begin());
                T[x] = idx_t(D + c);
            }

            T[dict.word_start[s + 1] - 1] = idx_t(s);
        }

        build_integer_sa(T.get(), SA.get(), size, alphabet, free_space, threads);
        T.reset();
        resize_malloc_array(SA, size);
    }

    log_phase_end(prm.log, time);

    // Phi is computed in the place of SA if SA is stored in a file
    dsa.sa_file = temporary_array(file, temporary_arrays_in_memory(size * sizeof(idx_t)));
    malloc_array_t<idx_t> phi_array;

    if (dsa.sa_file.is_in_memory()) {
        phi_array = new_malloc_array<idx_t>(size);
        dsa.sa_file.store(std::move(SA), size);
    } else {
        log_phase_start(prm.log, time, "storing data structures to disk");
        dsa.sa_file.write(SA.get(), size * sizeof(idx_t));
        log_phase_end(prm.log, time);
        phi_array = std::move(SA);
    }

    /* thread t of p threads processes the positions (or ranks) [64 * first_bit_block(t, p),
       64 * first_bit_block(t + 1, p)) */
    auto first_bit_block = [&](uint64_t t, uint64_t p) { return (t * num_bit_blocks) / p; };
    auto first_pos = [&](uint64_t t, uint64_t p) { return std::min(size, 64 * first_bit_block(t, p)); };

    // Phi[SA[i]] = SA[i-1] for i in [1, size)
    log_phase_start(prm.log, time, "building Phi of " + level_name("D"));
    idx_t* phi = phi_array.get();
    phase_progress phi_prog(prm.log, size);

    run_on_threads(threads, [&](uint64_t t, uint64_t p) {
        const uint64_t beg = std::max<uint64_t>(1, first_pos(t, p));
        const uint64_t end = first_pos(t + 1, p);
        if (beg >= end) return;
        array_reader in(dsa.sa_file, sizeof(idx_t));
        uint64_t prev = in[beg - 1];

        for (uint64_t base = beg; base < end; base += phi_prog.grain()) {
            const uint64_t stop = std::min(end, base + phi_prog.grain());

            for (uint64_t i = base; i < stop; i++) {
                const uint64_t x = in[i];
                phi[x] = idx_t(prev);
                prev = x;
            }

            phi_prog.advance(stop - base);
        }
    });

    phi_prog.finish();
    log_phase_end(prm.log, time);

    /* bit x of equals_phi_bits = 1 <=> the suffix at position x is longer than w and equals the suffix Phi[x] up to
       their separators; l is the length of their common prefix, which decreases by at most one from one position to
       the next (as in the computation of the PLCP array), and is 0 at a separator */
    log_phase_start(prm.log, time, "computing the equal suffixes of " + level_name("D"));
    malloc_array_t<uint64_t> equals_phi_bits = new_malloc_array<uint64_t>(num_bit_blocks);
    auto sym = [&](uint64_t x) { return uint64_t(text[x]); };
    phase_progress equal_prog(prm.log, num_bit_blocks);

    #pragma omp parallel num_threads(threads)
    {
        const uint64_t p = omp_get_num_threads();
        const uint64_t t = omp_get_thread_num();
        uint64_t l = 0;

        for (uint64_t b = first_bit_block(t, p); b < first_bit_block(t + 1, p); b++) {
            uint64_t bits = 0;

            for (uint64_t x = 64 * b; x < std::min(size, 64 * (b + 1)); x++) {
                if (sym(x) == separator) { l = 0; continue; }
                const uint64_t y = uint64_t(phi[x]);
                while (sym(x + l) != separator && sym(x + l) == sym(y + l)) l++;
                if (l > w && sym(x + l) == separator && sym(y + l) == separator) bits |= uint64_t { 1 } << (x % 64);
                if (l > 0) l--;
            }

            equals_phi_bits[b] = bits;
            equal_prog.at(b);
        }
    }

    equal_prog.finish();
    phi_array.reset();
    log_phase_end(prm.log, time);

    // bit i of dsa.equals_previous_bits is bit SA[i] of equals_phi_bits
    log_phase_start(prm.log, time, "reading data structures from disk");
    dsa.equals_previous_bits = new_malloc_array<uint64_t>(num_bit_blocks);
    phase_progress read_prog(prm.log, size);

    run_on_threads(threads, [&](uint64_t t, uint64_t p) {
        for (uint64_t b = first_bit_block(t, p); b < first_bit_block(t + 1, p); b++) dsa.equals_previous_bits[b] = 0;
        const uint64_t beg = first_pos(t, p);
        const uint64_t end = first_pos(t + 1, p);
        if (beg >= end) return;
        array_reader in(dsa.sa_file, sizeof(idx_t));

        for (uint64_t base = beg; base < end; base += read_prog.grain()) {
            const uint64_t stop = std::min(end, base + read_prog.grain());

            for (uint64_t i = base; i < stop; i++) {
                const uint64_t x = in[i];
                dsa.equals_previous_bits[i / 64] |= ((equals_phi_bits[x / 64] >> (x % 64)) & 1) << (i % 64);
            }

            read_prog.advance(stop - base);
        }
    });

    read_prog.finish();
    log_phase_end(prm.log, time);
    equals_phi_bits.reset();

    // mark the starts of the words
    dsa.word_starts = new_malloc_array<uint64_t>(num_bit_blocks);
    dsa.word_starts_before = uint_array(num_bit_blocks, value_bytes(D));

    #pragma omp parallel for num_threads(threads)
    for (uint64_t i = 0; i < num_bit_blocks; i++) dsa.word_starts[i] = 0;

    for (uint64_t s = 0; s < dict.num_words; s++) {
        dsa.word_starts[dict.word_start[s] / 64] |= uint64_t { 1 } << (dict.word_start[s] % 64);
    }

    uint64_t starts = 0;

    for (uint64_t i = 0; i < num_bit_blocks; i++) {
        dsa.word_starts_before.set(i, starts);
        starts += std::popcount(dsa.word_starts[i]);
    }

    return dsa;
}

} // namespace pfp_rlbwt

#endif // PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_SUFFIX_ARRAYS_HPP
