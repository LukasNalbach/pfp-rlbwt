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

#ifndef PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_RECURSION_TPP
#define PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_RECURSION_TPP

#include <pfp_rlbwt/construction/recursion.hpp>

namespace pfp_rlbwt {

template <typename idx_t>
inline void build_bwt_runs(parse_sa<idx_t>& psa, const level_output& out, uint64_t J, uint16_t threads)
{
    psa.bwt_runs = run_start_bits(J + 1, threads);

    #pragma omp parallel for num_threads(threads)
    for (uint64_t k = 0; k < out.r; k++) psa.bwt_runs.mark(record_field(out, k, 0));

    psa.bwt_runs.finish();
}

template <typename run_fnc_t>
inline void for_each_run_by_value(const level_output& out, uint64_t J, uint64_t num_values, uint16_t threads, run_fnc_t run)
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

template <typename idx_t, typename sym_t>
inline void build_psi_and_run_values(parse_sa<idx_t>& psa, const phrase_arrays<idx_t, sym_t>& arr, const level_output& out,
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

template <typename idx_t, typename sym_t>
inline void build_sampled_ends(parse_sa<idx_t>& psa, phrase_arrays<idx_t, sym_t>& arr, const level_output& out,
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

template <typename idx_t, typename sym_t>
inline parse_sa<idx_t> sort_parse_recursively(phrase_arrays<idx_t, sym_t>& arr, uint64_t J, uint64_t w,
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

template <typename idx_t, typename sym_t>
inline void build_psi_from_runs(parse_sa<idx_t>& psa, phrase_arrays<idx_t, sym_t>& arr, level_output& out, uint64_t J,
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

#endif // PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_RECURSION_TPP
