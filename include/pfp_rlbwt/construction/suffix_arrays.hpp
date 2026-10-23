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
#include <cstring>
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
#include <omp.h>

#include <pfp_rlbwt/construction/parse.hpp>
#include <pfp_rlbwt/construction/phrase_runs.hpp>
#include <pfp_rlbwt/types.hpp>

namespace pfp_rlbwt {

/* Psi of the suffix array SA of the parse of the phrase runs (SA[x] = the phrase run starting the suffix of rank x,
   SA[0] = J), and the starts C of the ranges of ranks of the suffixes per first word; for the words with a phrase run
   of two or more phrases, also the order of the suffixes of P (the parse of the phrases) that start with them. The
   suffixes of P that start at the first phrase (the head) of a phrase run are ordered as the suffixes of the parse of
   the phrase runs (see phrase_runs.hpp). */
template <typename idx_t>
struct parse_sa {
    // Psi: [0..J] psi[x] = ISA[SA[x]+1], the rank of the suffix one phrase run later (psi[0] = ISA[0])
    index_array<idx_t> psi;
    /* [0..D+1] C[c] = number of suffixes of the parse of the phrase runs that start with a value < c (value 0 is the
       end, value s+1 is word s); C[0] = 0, C[1] = 1, C[D+1] = J+1 */
    uint_array C;
    uint64_t num_words = 0; // D, the number of words
    /* the ranks in blocks of 2^rank_shift ranks, about one block per word: [0..(J >> rank_shift) + 1]
       block_word[t] = the word that the suffix of rank min(max(1, t * 2^rank_shift), J) starts with */
    std::vector<uint32_t> block_word;
    uint8_t rank_shift = 0;
    // the words with a phrase run of two phrases or more (repeated words), in increasing order
    std::vector<repeated_word_ranks> repeated;
    ranked_bits repeated_set; // the set of the repeated words
    // the set of the ranks of the first runs of the repeated words (see repeated_word_ranks::first_run)
    ranked_bits first_runs;
    /* [0..repeated.size()] extra_before[y] = number of phrases minus number of phrase runs of the words
       repeated[0..y-1] */
    std::vector<uint64_t> extra_before;

    // returns the word that the suffix of rank x >= 1 starts with
    inline uint64_t word_of_rank(uint64_t x) const
    {
        const uint64_t t = x >> rank_shift;
        uint64_t lo = block_word[t] + 1, hi = block_word[t + 1] + 2; // C[lo] <= x < C[hi]

        while (hi - lo > 1) {
            const uint64_t mid = (lo + hi) / 2;
            if (C[mid] <= x) lo = mid;
            else hi = mid;
        }

        return lo - 1;
    }

    // returns the ranks of the suffixes of P that start with word s, or nullptr if s is not repeated
    inline const repeated_word_ranks* find_repeated(uint64_t s) const
    {
        if (repeated.empty() || !repeated_set.contains(s)) return nullptr;
        return &repeated[repeated_set.rank(s)];
    }

    // returns the rank among the suffixes of P of the suffix starting at the head of the phrase run of rank x
    inline uint64_t rank_in_P(uint64_t x) const
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
     * @brief returns the rank among the suffixes of P of the suffix one phrase after the suffix with offset q among
     *        the suffixes of P that start with word s
     * @param s the word
     * @param rw find_repeated(s)
     * @param q [0..num_occ[s]-1] the offset
     */
    uint64_t next_rank_in_P(uint64_t s, const repeated_word_ranks* rw, uint64_t q) const
    {
        if (rw == nullptr) return rank_in_P(psi[C[s + 1] + q]);
        // the suffixes in layer 1 are followed by the head of the run after their run, the others by non-heads of s
        const uint64_t up_layer_1 = rw->occ - rw->up.num_runs;
        if (q < rw->down.num_runs) return rank_in_P(psi[rw->first_run + rw->down.runs[q]]);
        if (q >= up_layer_1) return rank_in_P(psi[rw->first_run + rw->up.runs[q - up_layer_1]]);
        return rw->first_rank + rw->non_head_offset(q - rw->down.num_runs);
    }
};

/**
 * @brief the generalized suffix array of the dictionary in a file, the positions where a suffix equals its
 *        predecessor and the word containing each position
 */
template <typename idx_t>
struct dictionary_sa {
    /* the file of the generalized suffix array of dict.text: [0..size-1] SA[i] = starting position of the i-th
       smallest suffix (sizeof(idx_t) bytes per value) */
    temporary_file sa_file;
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

    /**
     * @brief reads SA[beg..end-1] from sa_file
     * @param file a reader of sa_file (with sizeof(idx_t) bytes per value)
     * @param beg first rank
     * @param end rank after the last one
     * @param values receives SA[beg..end-1]
     */
    static void read_sa(array_file_reader& file, uint64_t beg, uint64_t end, no_init_vector<idx_t>& values)
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
 * @brief computes the order of the suffixes of P that start with a repeated word (see repeated_word_ranks)
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
        auto copies = [&](uint64_t x) {
            return (arr.ends_by_rank[psa.psi[x]] - arr.ends_by_rank[x]) / rw.advance;
        };

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
 * @brief computes Psi of the suffix array of the parse of the phrase runs and frees the parse; orders the ends and the
 *        symbols of the phrase runs by rank (see phrase_arrays::ends_by_rank and phrase_arrays::symbols_by_rank),
 *        reading them and the inverse suffix array, which is stored in a file meanwhile, from disk, and removes their
 *        files; for the repeated words also computes the order of the suffixes of P that start with them
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
    const uint64_t D = dict.num_words;
    parse_sa<idx_t> psa;
    psa.num_words = D;
    log_phase_start(prm.log, time, "building SA of P");
    const uint64_t free_space = free_space_per_value * (arr.num_letters + 1);
    malloc_array_t<idx_t> SA = new_malloc_array<idx_t>(J + 1 + free_space);
    build_integer_sa(arr.P.get(), SA.get(), J + 1, arr.num_letters + 1, free_space, threads);
    idx_t* shrunk = (idx_t*) std::realloc(SA.release(), (J + 1) * sizeof(idx_t));
    if (shrunk == nullptr) throw std::bad_alloc();
    SA.reset(shrunk);

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

    /* the inverse suffix array ISA, which is stored in the memory of P (P is not needed any more), and Psi from it,
       which is stored in the memory of SA */
    malloc_array_t<idx_t> ISA = std::move(arr.P);

    #pragma omp parallel for num_threads(threads)
    for (uint64_t x = 0; x <= J; x++) ISA[uint64_t(SA[x])] = idx_t(x);

    #pragma omp parallel for num_threads(threads)
    for (uint64_t x = 0; x <= J; x++) SA[x] = ISA[uint64_t(SA[x]) == J ? 0 : uint64_t(SA[x]) + 1];

    psa.psi = index_array<idx_t>(std::move(SA), J + 1, J);
    psa.rank_shift = J + 1 > D ? std::bit_width((J + 1) / D) - 1 : 0;
    psa.block_word.resize((J >> psa.rank_shift) + 2);

    #pragma omp parallel for num_threads(threads)
    for (uint64_t t = 0; t < psa.block_word.size(); t++) {
        const uint64_t x = std::min(std::max<uint64_t>(1, t << psa.rank_shift), J);
        uint64_t lo = 1, hi = D + 1; // C[lo] <= x < C[hi]

        while (hi - lo > 1) {
            const uint64_t mid = (lo + hi) / 2;
            if (psa.C[mid] <= x) lo = mid;
            else hi = mid;
        }

        psa.block_word[t] = uint32_t(lo - 1);
    }

    log_phase_end(prm.log, time);

    log_phase_start(prm.log, time, "spilling ISA of P");
    temporary_file isa_file = write_temporary_file(file, ISA.get(), (J + 1) * sizeof(idx_t));
    ISA.reset();
    log_spilled_size(prm.log, file);
    log_phase_end(prm.log, time);

    /* ends_by_rank[ISA[v]] = F_{v-1} and symbols_by_rank[ISA[v]] = value v-2 of symbols_file (see phrase_arrays) for
       v in [0..J]; thread t of p threads reads ISA[v], F_{v-1} and value v-2 from their files for
       v in [t * (J+1) / p, (t+1) * (J+1) / p), with readers that together read direct_reader::block bytes at once */
    log_phase_start(prm.log, time, "reading ISA of P, permuting the phrase runs");
    arr.ends_by_rank = uint_array(J + 1, arr.position_bytes);
    arr.symbols_by_rank = new_malloc_array<sym_t>(J + 1);

    run_on_threads(threads, [&](uint64_t t, uint64_t p) {
        const uint64_t beg = (t * (J + 1)) / p, end = ((t + 1) * (J + 1)) / p;
        if (beg >= end) return;
        constexpr uint64_t num_files = 3;
        constexpr uint64_t read_bytes = direct_reader::block / num_files;
        array_file_reader isa(isa_file.path(), sizeof(idx_t), read_bytes);
        array_file_reader ends(arr.ends_file.path(), arr.position_bytes, read_bytes);
        array_file_reader symbols(arr.symbols_file.path(), sizeof(sym_t), read_bytes);

        for (uint64_t v = beg; v < end; v++) {
            const uint64_t x = isa[v];
            arr.ends_by_rank.set(x, v == 0 ? w - 1 : ends[v - 1]);
            if (v <= 1) arr.symbols_by_rank[x] = arr.last_symbol;
            else std::memcpy(&arr.symbols_by_rank[x], symbols.bytes(v - 2), sizeof(sym_t));
        }
    });

    log_spilled_size(prm.log, file);
    isa_file.remove();
    arr.ends_file.remove();
    arr.symbols_file.remove();
    if (!arr.repeated_words.empty()) build_repeated_word_ranks(psa, arr, J, w, dict, threads);
    log_phase_end(prm.log, time);
    return psa;
}

/**
 * @brief computes the generalized suffix array of the dictionary, which is stored in a file, and the positions where
 *        a suffix of length > w equals its predecessor; Phi is computed in the place of the suffix array, and the
 *        equal suffixes are found with Phi as in the computation of the PLCP array
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

    log_phase_start(prm.log, time, "building SA of D");

    if constexpr (sizeof(sym_t) == 1) {
        SA = new_malloc_array<idx_t>(size + byte_alphabet_free_space);
        build_generalized_sa(text, SA.get(), size, threads);
    } else {
        /* the separator of word s is s (no two separators are equal) and symbol c is D + c, or D + the rank of c among
           the symbols of the dictionary if the symbols are larger than the dictionary */
        SA = new_malloc_array<idx_t>(size);
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

        build_integer_sa(T.get(), SA.get(), size, dictionary_alphabet_size(dict, max_symbol), 0, threads);
    }

    log_phase_end(prm.log, time);

    log_phase_start(prm.log, time, "spilling SA of D");
    dsa.sa_file = write_temporary_file(file, SA.get(), size * sizeof(idx_t));
    log_spilled_size(prm.log, file);
    log_phase_end(prm.log, time);

    /* thread t of p threads processes the positions (or ranks) [64 * first_bit_block(t, p),
       64 * first_bit_block(t + 1, p)) */
    auto first_bit_block = [&](uint64_t t, uint64_t p) { return (t * num_bit_blocks) / p; };
    auto first_pos = [&](uint64_t t, uint64_t p) { return std::min(size, 64 * first_bit_block(t, p)); };

    // Phi[SA[i]] = SA[i-1] for i in [1, size), in the place of SA
    log_phase_start(prm.log, time, "building Phi of D");
    idx_t* phi = SA.get();

    run_on_threads(threads, [&](uint64_t t, uint64_t p) {
        const uint64_t beg = std::max<uint64_t>(1, first_pos(t, p));
        const uint64_t end = first_pos(t + 1, p);
        if (beg >= end) return;
        array_file_reader in(file, sizeof(idx_t));
        uint64_t prev = in[beg - 1];

        for (uint64_t i = beg; i < end; i++) {
            const uint64_t x = in[i];
            phi[x] = idx_t(prev);
            prev = x;
        }
    });

    log_phase_end(prm.log, time);

    /* bit x of equals_phi_bits = 1 <=> the suffix at position x is longer than w and equals the suffix Phi[x] up to
       their separators; l is the length of their common prefix, which decreases by at most one from one position to
       the next (as in the computation of the PLCP array), and is 0 at a separator */
    log_phase_start(prm.log, time, "finding the equal suffixes of D");
    malloc_array_t<uint64_t> equals_phi_bits = new_malloc_array<uint64_t>(num_bit_blocks);
    auto sym = [&](uint64_t x) { return uint64_t(text[x]); };

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
        }
    }

    SA.reset();
    log_phase_end(prm.log, time);

    // bit i of dsa.equals_previous_bits is bit SA[i] of equals_phi_bits
    log_phase_start(prm.log, time, "reading SA of D");
    dsa.equals_previous_bits = new_malloc_array<uint64_t>(num_bit_blocks);

    run_on_threads(threads, [&](uint64_t t, uint64_t p) {
        for (uint64_t b = first_bit_block(t, p); b < first_bit_block(t + 1, p); b++) dsa.equals_previous_bits[b] = 0;
        const uint64_t beg = first_pos(t, p);
        const uint64_t end = first_pos(t + 1, p);
        if (beg >= end) return;
        array_file_reader in(file, sizeof(idx_t));

        for (uint64_t i = beg; i < end; i++) {
            const uint64_t x = in[i];
            dsa.equals_previous_bits[i / 64] |= ((equals_phi_bits[x / 64] >> (x % 64)) & 1) << (i % 64);
        }
    });

    log_spilled_size(prm.log, file);
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
