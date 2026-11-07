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
 * @file recursion.hpp
 * @brief the recursive prefix-free parsing of the parse: the BWT runs of the levels
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
void build_bwt_runs(parse_sa<idx_t>& psa, const level_output& out, uint64_t J, uint16_t threads)
{
    psa.bwt_runs = run_start_bits(J + 1, threads);

    #pragma omp parallel for num_threads(threads)
    for (uint64_t k = 0; k < out.r; k++) psa.bwt_runs.mark(record_field(out, k, 0));

    psa.bwt_runs.finish();
}

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
    run_fnc_t run)
{
    const uint64_t r = out.r;
    // returns the start of run k (J + 1 for k = r)
    auto start = [&](uint64_t k) { return k == r ? J + 1 : record_field(out, k, 0); };

    /* the runs ordered by value and then by start: [0..r-1] order[j] = the j-th run in this order; [0..num_values]
       first_of_value[v] = number of runs with a value < v */
    std::vector<uint64_t> first_of_value(num_values + 1, 0);
    for (uint64_t k = 0; k < r; k++) first_of_value[record_value(out, k) + 1]++;
    for (uint64_t v = 0; v < num_values; v++) first_of_value[v + 1] += first_of_value[v];
    uint_array order(r, value_bytes(r));

    {
        std::vector<uint64_t> next = first_of_value;
        for (uint64_t k = 0; k < r; k++) order.set(next[record_value(out, k)]++, k);
    }

    first_of_value = std::vector<uint64_t>();

    // part t of p parts processes the runs order[part_first[t]..part_first[t+1]-1], from the rank part_rank[t] on
    const uint64_t p = std::max<uint64_t>(1, std::min<uint64_t>(threads, r));
    std::vector<uint64_t> part_first(p + 1), part_rank(p + 1, 0);
    for (uint64_t t = 0; t <= p; t++) part_first[t] = (t * r) / p;

    #pragma omp parallel for num_threads(p)
    for (uint64_t t = 0; t < p; t++) {
        uint64_t sum = 0;
        for (uint64_t j = part_first[t]; j < part_first[t + 1]; j++) sum += start(order[j] + 1) - start(order[j]);
        part_rank[t + 1] = sum;
    }

    for (uint64_t t = 0; t < p; t++) part_rank[t + 1] += part_rank[t];

    #pragma omp parallel for num_threads(p)
    for (uint64_t t = 0; t < p; t++) {
        uint64_t y = part_rank[t];

        for (uint64_t j = part_first[t]; j < part_first[t + 1]; j++) {
            const uint64_t k = order[j];
            const uint64_t a = start(k), len = start(k + 1) - a;
            run(t, j, k, a, len, y);
            y += len;
        }
    }
}

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
    uint64_t J, uint64_t w, const dictionary<sym_t>& dict, uint16_t threads)
{
    const uint64_t num_values = arr.num_letters + 1;
    // the run of position a with the value v is the image of the ranks [y, y + len), whose rows are the ones of v
    const value_rows by_value(arr, w, dict);
    mds_pairs pairs(out.r, J + 1);
    value_rows::rows_t rows = by_value.new_rows(out.r, threads);

    for_each_run_by_value(out, J, num_values, threads,
        [&](uint64_t, uint64_t j, uint64_t k, uint64_t a, uint64_t, uint64_t y) {
            pairs.set<0>(j, y);
            pairs.set<1>(j, a);
            by_value.set(rows, j, record_value(out, k));
        });

    psa.psi = psi_mds(std::move(pairs), J + 1, rows.first, rows.second, threads);
    build_word_intervals(psa, threads);
    psa.bwt_run_values = uint_array(out.r, value_bytes(num_values));

    #pragma omp parallel for num_threads(threads)
    for (uint64_t k = 0; k < out.r; k++) psa.bwt_run_values.set(k, record_value(out, k));
}

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
    uint16_t threads)
{
    sampled_ends& samples = psa.ends;
    const uint64_t r = out.r;
    const uint8_t bytes = arr.position_bytes;
    samples.num_runs = r;
    samples.repeated_offset.assign(ranges.size() + 1, 0);

    for (uint64_t i = 0; i < ranges.size(); i++) {
        samples.repeated_offset[i + 1] = samples.repeated_offset[i] + ranges[i].second - ranges[i].first;
    }

    const uint64_t num_repeated_ranks = samples.repeated_offset.back();
    /* the requests (v, slot), ordered by v: F_{v-1} and F_v for the first (slot 2k) and the last (slot 2k + 1) position
       of run k, whose SA values are v, and F_v for the ranks of the ranges (slots 2r, ...) */
    std::vector<std::pair<uint64_t, uint64_t>> requests(2 * r + num_repeated_ranks);

    #pragma omp parallel for num_threads(threads)
    for (uint64_t k = 0; k < r; k++) {
        requests[2 * k] = { record_field(out, k, out.layout.sa_start), 2 * k };
        requests[2 * k + 1] = { record_field(out, k, out.layout.sa_end), 2 * k + 1 };
    }

    for (uint64_t i = 0; i < num_repeated_ranks; i++) requests[2 * r + i] = { repeated_sa[i], 2 * r + i };
    radix_sort(requests.begin(), requests.end(), key_of_first {}, threads);
    samples.values = uint_array(4 * r, bytes);
    samples.repeated_ends = uint_array(num_repeated_ranks, bytes);

    // thread t of p threads answers the requests [t * size / p, (t+1) * size / p), reading the ends in increasing order
    run_on_threads(threads, [&](uint64_t t, uint64_t p) {
        const uint64_t beg = (t * requests.size()) / p, end = ((t + 1) * requests.size()) / p;
        if (beg >= end) return;
        array_reader ends(arr.ends_file, bytes);

        for (uint64_t i = beg; i < end; i++) {
            const auto [v, slot] = requests[i];
            // F_{v-1} and F_v, with F_{-1} = w - 1 and F_J = F_{J-1}
            const uint64_t prev_end = v == 0 ? w - 1 : ends[v - 1];
            const uint64_t cur_end = v == J ? prev_end : ends[v];

            if (slot >= 2 * r) {
                samples.repeated_ends.set(slot - 2 * r, cur_end);
            } else {
                samples.values.set(2 * slot, cur_end);
                samples.values.set(2 * slot + 1, prev_end);
            }
        }
    });

    arr.ends_file.remove();
    samples.repeated_ranges = std::move(ranges);
}

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
    std::chrono::steady_clock::time_point& time, level_output& out)
{
    const uint16_t threads = prm.num_threads;
    const uint64_t num_values = arr.num_letters + 1;
    parse_sa<idx_t> psa;
    build_word_ranges(psa, arr, J, dict);

    /* without the full SA and LCP array, the ends of the phrase runs by rank are sampled (see sampled_ends): the next
       level computes the SA samples and SA at the ranks of the runs of the repeated words */
    const bool sampled = !prm.write_sa && !prm.write_lcp;
    // the ranges [beg, end) of the ranks of the runs of the repeated words, in increasing order (if sampled)
    std::vector<std::pair<uint64_t, uint64_t>> repeated_ranges;
    uint64_t num_repeated_ranks = 0; // the number of ranks in repeated_ranges

    if (sampled) {
        for (const repeated_word& rw : arr.repeated_words) {
            repeated_ranges.emplace_back(psa.C[rw.word + 1], psa.C[rw.word + 2]);
            num_repeated_ranks += psa.C[rw.word + 2] - psa.C[rw.word + 1];
        }
    }

    // SA of P at the ranks of repeated_ranges, one range after another
    uint_array repeated_sa(num_repeated_ranks, value_bytes(J));
    params level_prm = prm;
    level_prm.window = prm.recursion_window;
    level_prm.modulus = prm.recursion_modulus;
    level_prm.position_bytes = 0;
    level_prm.format = file_format::plain;
    level_prm.write_sa_samples = sampled;
    level_prm.write_lcp_samples = false;
    level_prm.write_bwt = level_prm.write_lcp = false;
    level_prm.write_sa = !sampled || num_repeated_ranks != 0;
    level_prm.extra_sections.clear();
    const std::string level_prefix = prefix + "_" + std::to_string(current_level + 1);

    /* the next level writes SA of P to ends_by_rank (allocated once the word ids are freed), whose values become the
       ends below, or its values at the ranks of the runs of the repeated words to repeated_sa */
    out.sa = sampled ? &repeated_sa : &arr.ends_by_rank;
    out.sa_ranges = repeated_ranges;

    /* computes the next level with T' read with symbols of type text_t from the word ids of the phrase runs; it stops
       after its parsing if the PFP of T' is larger than half of T' (then the word ids are kept) */
    auto next_level = [&]<typename text_t>() {
        out.parsed = [&](const result& sizes) {
            if (2 * (sizes.num_phrases + sizes.dictionary_size) > J) return false;
            arr.free_word_ids();
            if (!sampled) arr.ends_by_rank = uint_array(J + 1, arr.position_bytes);
            return true;
        };

        build_level<text_t>(letter_source<idx_t, sym_t, text_t>(arr, J), level_prefix, level_prm, out);
    };

    // the type of the symbols of T' (they are at most num_values, see letter_source)
    if (num_values + 1 <= uint24_t::max) next_level.template operator()<uint24_t>();
    else if (num_values + 1 <= std::numeric_limits<uint32_t>::max()) next_level.template operator()<uint32_t>();
    else next_level.template operator()<uint64_t>();

    out.parsed = nullptr;
    out.sa = nullptr;
    if (!out.completed) return psa;

    if (prm.log) {
        log_out() << "r(BWT(" << level_name("P") << ")) = " << format_number(out.r) << ", |" << level_name("P") << "|/r = "
                  << format_ratio(double(J) / double(out.r)) << std::endl;
    }

    log_phase_start(prm.log, time, "sorting the phrase ends by SA of " + level_name("P"));

    build_bwt_runs(psa, out, J, threads);

    if (sampled) {
        build_sampled_ends(psa, arr, out, std::move(repeated_ranges), repeated_sa, J, w, threads);
        repeated_sa = uint_array();
        // the records keep only the starts and the symbols of the runs
        const record_layout with_samples = out.layout;
        out.layout = record_layout(with_samples.position_bytes, with_samples.symbol_bytes, false, false);

        for (uint64_t k = 0; k < out.r; k++) {
            std::memmove(out.records.data() + k * out.layout.bytes, out.records.data() + k * with_samples.bytes,
                out.layout.bytes);
        }

        out.records.resize(out.r * out.layout.bytes);
        out.records.shrink_to_fit();
    } else {
        // the ends in text order: [0..J-1] ends[i] = F_i
        uint_array ends;

        if (arr.ends_file.is_in_memory()) {
            ends = uint_array(arr.ends_file.release(), arr.position_bytes);
        } else {
            ends = uint_array(J, arr.position_bytes);

            run_on_threads(threads, [&](uint64_t t, uint64_t p) {
                const uint64_t beg = (t * J) / p, end = ((t + 1) * J) / p;
                if (beg >= end) return;
                array_reader in(arr.ends_file, arr.position_bytes);
                in.read(beg, end, ends.raw() + beg * arr.position_bytes);
            });

            arr.ends_file.remove();
        }

        // ends_by_rank[x] = SA[x] (written by the next level) becomes F_{SA[x]-1}, with F_{-1} = w - 1
        #pragma omp parallel for num_threads(threads)
        for (uint64_t x = 0; x <= J; x++) {
            const uint64_t v = arr.ends_by_rank[x];
            arr.ends_by_rank.set(x, v == 0 ? w - 1 : ends[v - 1]);
        }
    }

    log_phase_end(prm.log, time);
    return psa;
}

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
    uint64_t w, const dictionary<sym_t>& dict, const params& prm, std::chrono::steady_clock::time_point& time)
{
    const uint16_t threads = prm.num_threads;
    log_phase_start(prm.log, time, "building M_Psi of " + level_name("P") + " from its BWT");
    build_psi_and_run_values(psa, arr, out, J, w, dict, threads);
    arr.value_symbols = symbols_of_values(arr, w, dict);
    log_phase_end(prm.log, time);

    if (!arr.repeated_words.empty()) {
        log_phase_start(prm.log, time, "sorting the repeated phrases of " + level_name("P"));
        build_repeated_word_ranks(psa, arr, J, w, dict, threads);
        log_phase_end(prm.log, time);
    }

    out.records = no_init_vector<uint8_t>();
}

} // namespace pfp_rlbwt

#endif // PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_RECURSION_HPP
