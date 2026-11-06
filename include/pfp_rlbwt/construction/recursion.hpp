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

#ifndef PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_RECURSION_HPP
#define PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_RECURSION_HPP

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include <misc/log.hpp>
#include <pfp_rlbwt/construction/merge.hpp>

namespace pfp_rlbwt {

/* The recursive prefix-free parsing: the parse P of the phrase runs (see phrase_arrays::P; P[J] = 0, the end) is
   suffix-sorted with the prefix-free parsing of the text T' = (P[0] + 1) (P[1] + 1) ... (P[J-1] + 1), whose symbols are
   at least min_symbol (see letter_source). SA of T'$ is SA of P, and the BWT of T'$ is the BWT of P with each value
   v >= 1 replaced by v + 1 and the end replaced by the terminator. For each value v, Psi maps the ranks of the suffixes
   that start with v, in increasing order, to the positions of v in the BWT of P, in increasing order; hence Psi follows
   from the runs of this BWT (see for_each_run_by_value). */

/**
 * @brief computes the runs of the BWT of T$ in the plain format, with the SA samples if prm.write_sa_samples is
 *        set, and the suffix array of T$ (or its values in ranges, see level_output) if prm.write_sa is set, and keeps
 *        them in RAM (a level of the recursive prefix-free parsing; it is recursive itself if prm.recursive is set)
 * @param text the text T (all symbols >= min_symbol)
 * @param prefix path prefix of the temporary files
 * @param prm parameters (the format is plain; no LCP samples, full BWT, full LCP array or extra sections)
 * @param out receives the runs, the suffix array and their number; out.parsed can stop the level after its parsing
 */
template <typename sym_t>
void build_level(const text_source<sym_t>& text, const std::string& prefix, const params& prm, level_output& out);

// the levels after the first are compiled with the instantiations of their symbol types (see src/instantiations.cpp)
extern template void build_level<uint24_t>(const text_source<uint24_t>&, const std::string&, const params&,
    level_output&);
extern template void build_level<uint32_t>(const text_source<uint32_t>&, const std::string&, const params&,
    level_output&);
extern template void build_level<uint64_t>(const text_source<uint64_t>&, const std::string&, const params&,
    level_output&);

// returns the position field at offset of the record of run k of the BWT of T'$ (offset 0: the start of run k)
inline uint64_t record_field(const level_output& out, uint64_t k, uint8_t offset)
{
    return read_bytes(out.records.data() + k * out.layout.bytes + offset, out.layout.position_bytes);
}

// returns the value in P of run k of the BWT of T'$: its symbol minus one, or 0 (the end) for the terminator
inline uint64_t record_value(const level_output& out, uint64_t k)
{
    const uint64_t c = read_bytes(out.records.data() + k * out.layout.bytes + out.layout.symbol,
        out.layout.symbol_bytes);
    return c == terminator ? 0 : c - 1;
}

/**
 * @brief marks the first positions of the runs of the BWT of T'$ (see parse_sa::bwt_runs)
 * @param psa receives the marks
 * @param out the runs (see level_output)
 * @param J number of phrase runs
 * @param threads number of threads
 */
template <typename idx_t>
void build_bwt_runs(parse_sa<idx_t>& psa, const level_output& out, uint64_t J, uint16_t threads);

/**
 * @brief calls run(t, j, k, a, len, y) for the runs of the BWT of T'$ ordered by their values in P and then by their
 *        starts; in this order, they are the images of consecutive ranges of ranks under Psi (the suffixes that start
 *        with a value v follow those that start with a smaller value): the j-th run in this order is run k of the BWT,
 *        starts at position a, has length len and is the image of the ranks [y, y + len); thread t processes part t
 *        of the runs in this order
 * @param out the runs (see level_output)
 * @param J number of phrase runs
 * @param num_values number of values of P (num_letters + 1)
 * @param threads number of threads
 * @param run function called for each run
 */
template <typename run_fnc_t>
void for_each_run_by_value(const level_output& out, uint64_t J, uint64_t num_values, uint16_t threads,
    run_fnc_t run);

/**
 * @brief builds Psi (see psi_mds) and the values of the runs of the BWT of P (see parse_sa::bwt_run_values) from the
 *        runs of the BWT of T'$
 * @param psa receives Psi and the values
 * @param arr the phrase runs (the letters of the repeated words)
 * @param out the runs (see level_output)
 * @param J number of phrase runs
 * @param w length of a trigger string
 * @param dict the dictionary
 * @param threads number of threads
 */
template <typename idx_t, typename sym_t>
void build_psi_and_run_values(parse_sa<idx_t>& psa, const phrase_arrays<idx_t, sym_t>& arr, const level_output& out,
    uint64_t J, uint64_t w, const dictionary<sym_t>& dict, uint16_t threads);

/**
 * @brief computes the sampled ends of the phrase runs by rank (see sampled_ends) from the runs of the BWT of T'$ with
 *        their SA samples and from SA of T'$ at the ranks of the runs of the repeated words; reads the ends of the
 *        phrase runs in text order once and removes them
 * @param psa receives the sampled ends (psa.ends)
 * @param arr the phrase runs (the ends in text order)
 * @param out the runs with SA samples (see level_output)
 * @param ranges the ranges [beg, end) of the ranks of the runs of the repeated words, in increasing order
 * @param repeated_sa SA[y] for the ranks y of the ranges, one range after another
 * @param J number of phrase runs
 * @param w length of a trigger string
 * @param threads number of threads
 */
template <typename idx_t, typename sym_t>
void build_sampled_ends(parse_sa<idx_t>& psa, phrase_arrays<idx_t, sym_t>& arr, const level_output& out,
    std::vector<std::pair<uint64_t, uint64_t>> ranges, const uint_array& repeated_sa, uint64_t J, uint64_t w,
    uint16_t threads);

/**
 * @brief the first part of the recursive counterpart of build_parse_sa: suffix-sorts P with the prefix-free parsing of
 *        T' (see build_level), whose symbols are read from the word ids of the phrase runs (see letter_source), and
 *        computes the ends of the phrase runs by rank, sampled (see sampled_ends) or, with the full SA or LCP array, in
 *        phrase_arrays::ends_by_rank; frees the word ids and removes the ends of the phrase runs in text order
 * @param arr the phrase runs
 * @param J number of phrase runs
 * @param w length of a trigger string
 * @param dict the dictionary
 * @param prefix path prefix of the temporary files of the next level
 * @param prm parameters (threads, logging and the recursion)
 * @param time phase timer (see log_phase_start)
 * @param out receives the runs of the BWT of T'$ (see level_output), from which build_psi_from_runs computes Psi; if
 *        out.completed is false, the PFP of T' is larger than half of T', P is not sorted and the word ids are kept
 * @return the suffix array of the parse of the phrase runs without Psi
 */
template <typename idx_t, typename sym_t>
parse_sa<idx_t> sort_parse_recursively(phrase_arrays<idx_t, sym_t>& arr, uint64_t J, uint64_t w,
    const dictionary<sym_t>& dict, const std::string& prefix, const params& prm,
    std::chrono::steady_clock::time_point& time, level_output& out);

/**
 * @brief the second part of the recursive counterpart of build_parse_sa (see sort_parse_recursively): computes Psi of
 *        the suffix array of the parse of the phrase runs (see psi_mds), the values of the runs of its BWT and the
 *        symbols of the values from the runs of the BWT of T'$ and, for the repeated words, the order of the suffixes
 *        of the parse of the phrases that start with them (see build_repeated_word_ranks)
 * @param psa the suffix array of the parse of the phrase runs (see sort_parse_recursively)
 * @param arr the phrase runs
 * @param out the runs of the BWT of T'$; freed
 * @param J number of phrase runs
 * @param w length of a trigger string
 * @param dict the dictionary
 * @param prm parameters (threads and logging)
 * @param time phase timer (see log_phase_start)
 */
template <typename idx_t, typename sym_t>
void build_psi_from_runs(parse_sa<idx_t>& psa, phrase_arrays<idx_t, sym_t>& arr, level_output& out, uint64_t J,
    uint64_t w, const dictionary<sym_t>& dict, const params& prm, std::chrono::steady_clock::time_point& time);

} // namespace pfp_rlbwt

#include <pfp_rlbwt/construction/recursion.tpp>

#endif // PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_RECURSION_HPP
