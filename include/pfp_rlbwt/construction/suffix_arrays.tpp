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

#ifndef PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_SUFFIX_ARRAYS_TPP
#define PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_SUFFIX_ARRAYS_TPP

#include <pfp_rlbwt/construction/suffix_arrays.hpp>

namespace pfp_rlbwt {

inline void run_start_bits::finish()
{
    constexpr uint64_t line_bits = cache_line_rank_bit_vector::bits_per_line;
    uint64_t marks = 0;

    for (uint64_t l = 0; l <= bits.size() / line_bits; l++) {
        bits.set_line_rank(l, marks);
        const uint64_t* line = bits.line(l * line_bits);
        for (uint64_t j = 1; j < cache_line_rank_bit_vector::words_per_line; j++) marks += std::popcount(line[j]);
    }
}

inline uint64_t run_start_bits::next_start(uint64_t y) const
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

inline uint64_t sampled_ends::repeated_end(uint64_t y) const
{
    const auto it = std::upper_bound(repeated_ranges.begin(), repeated_ranges.end(), std::make_pair(y, no_end));
    if (it == repeated_ranges.begin()) return no_end;
    const uint64_t i = it - repeated_ranges.begin() - 1;
    if (y >= repeated_ranges[i].second) return no_end;
    return repeated_ends[repeated_offset[i] + y - repeated_ranges[i].first];
}

inline uint64_t sampled_ends::sampled_end(const run_start_bits& bwt_runs, uint64_t y, uint64_t psi_y) const
{
    uint64_t end = boundary_end(bwt_runs, y, 0);
    if (end != no_end) return end;
    // F_{SA[y]} = F_{SA[psi[y]]-1}
    end = boundary_end(bwt_runs, psi_y, 1);
    if (end != no_end) return end;
    return repeated_end(y);
}

inline uint64_t sampled_ends::end_at(const run_start_bits& bwt_runs, uint64_t y, uint64_t psi_y, const psi_mds& psi) const
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

template <typename idx_t>
inline uint64_t parse_sa<idx_t>::rank_among_phrases(uint64_t x) const
{
    // the number of repeated words whose runs start at a rank <= x
    const uint64_t y = first_runs.rank(x + 1);

    if (y > 0) {
        const repeated_word_ranks& rw = repeated[y - 1];
        if (x < rw.first_run + rw.num_runs) return rw.first_rank + rw.head_offset[x - rw.first_run];
    }

    return x + extra_before[y];
}

template <typename idx_t>
inline uint64_t parse_sa<idx_t>::next_rank_among_phrases(uint64_t s, const repeated_word_ranks* rw, uint64_t q) const
{
    if (rw == nullptr) return rank_among_phrases(psi[C[s + 1] + q]);
    // the suffixes in layer 1 are followed by the head of the run after their run, the others by non-heads of s
    const uint64_t up_layer_1 = rw->occ - rw->up.num_runs;
    if (q < rw->down.num_runs) return rank_among_phrases(psi[rw->first_run + rw->down.runs[q]]);
    if (q >= up_layer_1) return rank_among_phrases(psi[rw->first_run + rw->up.runs[q - up_layer_1]]);
    return rw->first_rank + rw->non_head_offset(q - rw->down.num_runs);
}

template <typename idx_t>
inline void build_integer_sa(idx_t* T, idx_t* SA, uint64_t n, uint64_t k, uint64_t fs, uint16_t threads)
{
    int64_t ret;
    if constexpr (std::is_same_v<idx_t, int32_t>)
        ret = libsais_int_omp(T, SA, n, k, std::min<uint64_t>(fs, std::numeric_limits<int32_t>::max() - n), threads);
    else if constexpr (std::is_same_v<idx_t, sa_int40_t>)
        ret = libsais40_impl::libsais40_long_omp(T, SA, n, k, fs, threads);
    else ret = libsais48_impl::libsais48_long_omp(T, SA, n, k, fs, threads);
    if (ret != 0) throw std::runtime_error("pfp_rlbwt: libsais failed on the parsing");
}

template <typename idx_t>
inline void build_generalized_sa(const uint8_t* T, idx_t* SA, uint64_t n, uint16_t threads)
{
    int64_t ret;
    if constexpr (std::is_same_v<idx_t, int32_t>)
        ret = libsais_gsa_omp(T, SA, n, byte_alphabet_free_space, nullptr, threads);
    else if constexpr (std::is_same_v<idx_t, sa_int40_t>)
        ret = libsais40_impl::libsais40_gsa_omp(T, SA, n, byte_alphabet_free_space, nullptr, threads);
    else ret = libsais48_impl::libsais48_gsa_omp(T, SA, n, byte_alphabet_free_space, nullptr, threads);
    if (ret != 0) throw std::runtime_error("pfp_rlbwt: libsais failed on the dictionary");
}

template <typename idx_t, typename sym_t>
inline void build_repeated_word_ranks(parse_sa<idx_t>& psa, phrase_arrays<idx_t, sym_t>& arr, uint64_t J, uint64_t w,
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

template <typename idx_t, typename sym_t>
inline void build_word_ranges(parse_sa<idx_t>& psa, const phrase_arrays<idx_t, sym_t>& arr, uint64_t J, const dictionary<sym_t>& dict)
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

template <typename idx_t, typename sym_t, typename value_fnc_t>
inline void for_each_value(const phrase_arrays<idx_t, sym_t>& arr, uint64_t D, value_fnc_t f)
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

template <typename idx_t, typename sym_t>
inline std::vector<sym_t> symbols_of_values(const phrase_arrays<idx_t, sym_t>& arr, uint64_t w, const dictionary<sym_t>& dict)
{
    std::vector<sym_t> symbols(arr.num_letters + 1);
    symbols[0] = arr.last_symbol;

    for_each_value(arr, dict.num_words, [&](uint64_t v, uint64_t s, uint64_t, bool) {
        symbols[v] = dict.text[dict.word_start[s] + dict.length(s) - w - 1];
    });

    return symbols;
}

template <typename idx_t, typename sym_t>
inline value_rows::value_rows(const phrase_arrays<idx_t, sym_t>& arr, uint64_t w, const dictionary<sym_t>& dict)
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

template <typename idx_t>
inline void build_word_intervals(parse_sa<idx_t>& psa, uint16_t threads)
{
    const uint64_t num = psa.num_words + 2;
    psa.C_intervals = uint_array(num, value_bytes(psa.psi.num_intervals()));

    #pragma omp parallel for num_threads(threads)
    for (uint64_t c = 0; c < num; c++) psa.C_intervals.set(c, psa.psi.interval_of(psa.C[c]));
}

template <typename idx_t, typename sym_t>
inline void build_psi_from_array(parse_sa<idx_t>& psa, phrase_arrays<idx_t, sym_t>& arr, const idx_t* psi, uint64_t J,
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

template <typename idx_t, typename sym_t>
inline parse_sa<idx_t> build_parse_sa(phrase_arrays<idx_t, sym_t>& arr, uint64_t J, uint64_t w,
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

template <typename idx_t, typename sym_t>
inline dictionary_sa<idx_t> build_dictionary_sa(const dictionary<sym_t>& dict, uint64_t max_symbol, uint64_t w,
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

#endif // PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_SUFFIX_ARRAYS_TPP
