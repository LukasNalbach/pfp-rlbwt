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
 * @file suffix_arrays.hpp
 * @brief the suffix arrays of the dictionary and of the parse, and Psi as a move data structure
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
    void finish();

    // returns whether position y is the first position of a run
    inline bool is_start(uint64_t y) const { return bits.bit(y); }

    // returns the run containing position y
    inline uint64_t run_of(uint64_t y) const { return bits.rank(y + 1) - 1; }

    // fetches the cache line that run_of(y) reads
    inline void prefetch(uint64_t y) const { __builtin_prefetch(bits.line(y + 1)); }

    // returns the first position after position y < u that starts a run, or u if there is none
    inline uint64_t next_start(uint64_t y) const;
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
    inline uint64_t repeated_end(uint64_t y) const;

    /* returns F_{SA[y]} if y or psi_y = psi[y] is the first or the last position of a run of the BWT of the parse or y
       is a rank of a run of a repeated word, else no_end */
    inline uint64_t sampled_end(const run_start_bits& bwt_runs, uint64_t y, uint64_t psi_y) const;

    /* returns F_{SA[y]} (psi_y = psi[y]); if it is not sampled, F_{SA[z]} minus the advances (see psi_mds) of the
       phrase runs SA[y] + 1, ..., SA[z] for the first rank z = psi^i[y] (i >= 1) whose end is sampled */
    inline uint64_t end_at(const run_start_bits& bwt_runs, uint64_t y, uint64_t psi_y, const psi_mds& psi) const;
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
    inline uint64_t rank_among_phrases(uint64_t x) const;

    /**
     * @brief returns the rank among the suffixes of the parse of the phrases of the suffix one phrase after the suffix
     *        with offset q among those that start with word s
     * @param s the word
     * @param rw find_repeated(s)
     * @param q [0..num_occ[s]-1] the offset
     */
    uint64_t next_rank_among_phrases(uint64_t s, const repeated_word_ranks* rw, uint64_t q) const;
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
void build_integer_sa(idx_t* T, idx_t* SA, uint64_t n, uint64_t k, uint64_t fs, uint16_t threads);

// computes the generalized suffix array of the string set T[0..n-1] (separated by 0)
template <typename idx_t>
void build_generalized_sa(const uint8_t* T, idx_t* SA, uint64_t n, uint16_t threads);

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
    const dictionary<sym_t>& dict, uint16_t threads);

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
    const dictionary<sym_t>& dict);

/**
 * @brief calls f(v, s, m, repeated) for each value v >= 1 of the parse of the phrase runs (the letter v - 1, see
 *        phrase_arrays::P), whose phrase runs have the word s and m phrases each; repeated <=> s is a repeated word
 *        (it has one value per key of its phrase runs, see run_key)
 * @param arr the phrase runs (the keys of the repeated words)
 * @param D the number of words
 * @param f the function
 */
template <typename idx_t, typename sym_t, typename value_fnc_t>
void for_each_value(const phrase_arrays<idx_t, sym_t>& arr, uint64_t D, value_fnc_t f);

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
    const dictionary<sym_t>& dict);

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
    value_rows(const phrase_arrays<idx_t, sym_t>& arr, uint64_t w, const dictionary<sym_t>& dict);

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
void build_word_intervals(parse_sa<idx_t>& psa, uint16_t threads);

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
    uint64_t w, const dictionary<sym_t>& dict, uint16_t threads);

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
    std::chrono::steady_clock::time_point& time);

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
    const std::string& file, const params& prm, std::chrono::steady_clock::time_point& time);

} // namespace pfp_rlbwt

#include <pfp_rlbwt/construction/suffix_arrays.tpp>

#endif // PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_SUFFIX_ARRAYS_HPP
