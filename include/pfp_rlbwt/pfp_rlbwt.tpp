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

#ifndef PFP_RLBWT_PFP_RLBWT_PFP_RLBWT_TPP
#define PFP_RLBWT_PFP_RLBWT_PFP_RLBWT_TPP

#include <pfp_rlbwt/pfp_rlbwt.hpp>

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <misc/files.hpp>
#include <misc/log.hpp>

#include <pfp_rlbwt/construction/recursion.hpp>

namespace pfp_rlbwt {

// returns the name of the text of the current level in the log: T, P, P_2, ... (see current_level)
inline std::string level_text_name()
{
    return current_level == 1 ? "T" : current_level == 2 ? "P" : "P_" + std::to_string(current_level - 1);
}

// returns the name of the phase that computes the outputs selected in prm (see merge)
inline std::string merge_phase_name(const params& prm, bool writes_run_data)
{
    std::vector<std::string> parts;

    if (writes_run_data) {
        parts.push_back(prm.format == file_format::compressed ? "compressed RLBWT" : "RLBWT");
        if (prm.write_sa_samples && prm.write_lcp_samples) parts.push_back("SA- and LCP-samples");
        else if (prm.write_sa_samples) parts.push_back("SA-samples");
        else if (prm.write_lcp_samples) parts.push_back("LCP-samples");
    }

    if (prm.write_bwt) parts.push_back("BWT");
    if (prm.write_sa) parts.push_back("SA");
    if (prm.write_lcp) parts.push_back("LCP");
    std::string name = "computing";

    for (uint64_t i = 0; i < parts.size(); i++)
        name += (i == 0 ? " " : (i + 1 == parts.size() ? " and " : ", ")) + parts[i];

    return current_level == 1 ? name : name + " of " + level_text_name();
}

/**
 * @brief builds the output from the parsing
 * @tparam idx_t type of the values of the suffix arrays of the parsing and of the dictionary (int32_t, sa_int40_t
 *         or sa_int48_t)
 * @param prs the parsing; its phrases are freed
 * @param prm parameters
 * @param prefix path prefix of the output files
 * @param level_out if not nullptr, receives the runs and the suffix array (see build_level)
 * @return the number r of runs
 */
template <typename sym_t, typename idx_t>
uint64_t build_from_parsing(parsing<sym_t>& prs, const params& prm, const std::string& prefix,
    level_output* level_out)
{
    std::chrono::steady_clock::time_point time;
    const uint16_t p = prm.num_threads;
    const uint64_t J = prs.num_phrase_runs;

    // P is suffix-sorted with its prefix-free parsing if J >= recursion_min_phrase_runs and J <= n / 2
    bool recursive = prm.recursive && J >= prm.recursion_min_phrase_runs && 2 * J <= prs.n;

    log_phase_start(prm.log, time, "building " + level_name("P"));
    phrase_arrays<idx_t, sym_t> arr = build_phrase_arrays<idx_t, sym_t>(prs, prefix, p, recursive);
    log_phase_end(prm.log, time);

    /* the suffix array of the dictionary is built before the one of P and stays on disk; if P is suffix-sorted
       recursively, the next level runs before (it frees the word ids of the phrase runs) and Psi is built after */
    parse_sa<idx_t> psa;
    level_output next_level; // the runs of the BWT of P, if P is suffix-sorted recursively
    if (recursive) {
        psa = sort_parse_recursively<idx_t, sym_t>(arr, J, prs.w, prs.dict, prefix, prm, time, next_level);

        // the PFP of P does not halve P: P is suffix-sorted directly
        if (!next_level.completed) {
            recursive = false;
            psa = parse_sa<idx_t>();
            log_phase_start(prm.log, time, "storing data structures to disk");
            store_parse(arr, J, p);
            log_phase_end(prm.log, time);
        }
    }

    dictionary_sa<idx_t> dsa = build_dictionary_sa<idx_t, sym_t>(prs.dict, prs.max_symbol, prs.w,
        prefix + ext_dictionary_sa, prm, time);

    if (recursive) build_psi_from_runs<idx_t, sym_t>(psa, arr, next_level, J, prs.w, prs.dict, prm, time);
    else psa = build_parse_sa<idx_t, sym_t>(arr, J, prs.w, prs.dict, prefix + ext_parse_isa, prm, time);

    uint64_t r;

    if (prm.write_lcp) {
        /* the full LCP array is computed from PLCP sampled at the SA values of the run starts: the first pass writes
           the outputs of the runs and collects these samples, the second pass writes the full arrays */
        params run_prm = prm;
        run_prm.write_bwt = run_prm.write_sa = run_prm.write_lcp = false;
        std::vector<run_start_values> run_starts;
        log_phase_start(prm.log, time, merge_phase_name(run_prm, true));
        r = merge<sym_t, idx_t>(prs, dsa, psa, arr, run_prm, prefix, true, &run_starts);
        log_phase_end(prm.log, time);

        log_phase_start(prm.log, time, "sampling PLCP");
        const sampled_plcp plcp = build_sampled_plcp(run_starts, prs.n + 1, p);
        log_phase_end(prm.log, time);

        log_phase_start(prm.log, time, merge_phase_name(prm, false));
        merge<sym_t, idx_t>(prs, dsa, psa, arr, prm, prefix, false, nullptr, &plcp);
        log_phase_end(prm.log, time);
    } else {
        log_phase_start(prm.log, time, merge_phase_name(prm, true));
        r = merge<sym_t, idx_t>(prs, dsa, psa, arr, prm, prefix, true, nullptr, nullptr, level_out);
        log_phase_end(prm.log, time);
    }

    return r;
}

// checks the parameters that do not depend on the length of the text (see build)
inline void check_params(const params& prm)
{
    if (prm.window == 0 || prm.modulus == 0 || prm.num_threads == 0)
        throw std::invalid_argument("pfp_rlbwt: invalid params");
    if (prm.recursive && (prm.recursion_window == 0 || prm.recursion_modulus == 0))
        throw std::invalid_argument("pfp_rlbwt: invalid recursion params");
    if (prm.runs_per_block == 0 || prm.runs_per_block % byte_bits != 0)
        throw std::invalid_argument("pfp_rlbwt: runs_per_block is not a positive multiple of 8");
    if (prm.runs_per_symbol_sample == 0 || prm.runs_per_symbol_sample > prm.runs_per_block)
        throw std::invalid_argument("pfp_rlbwt: runs_per_symbol_sample is not in [1..runs_per_block]");
    if (prm.lcp_values_per_block == 0) throw std::invalid_argument("pfp_rlbwt: lcp_values_per_block is not positive");

    for (uint64_t i = 0; i < prm.extra_sections.size(); i++) {
        const uint64_t id = prm.extra_sections[i].id;
        if (id < first_extra_section) throw std::invalid_argument("pfp_rlbwt: an extra section has an identifier "
            "smaller than first_extra_section");
        for (uint64_t j = 0; j < i; j++) if (prm.extra_sections[j].id == id)
            throw std::invalid_argument("pfp_rlbwt: two extra sections have the same identifier");
    }
}

/**
 * @brief checks that a text of length n fits the parameters and chooses position_bytes if it is 0
 * @param n the length of the text
 * @param prm the parameters (position_bytes is set)
 */
inline void fit_text_length(uint64_t n, params& prm)
{
    // the text positions and the window fit into max_value_bytes bytes: n + window < 2^(8 * max_value_bytes)
    constexpr uint64_t position_limit = uint64_t { 1 } << (8 * max_value_bytes);
    if (prm.window >= position_limit || n >= position_limit - prm.window)
        throw std::invalid_argument("pfp_rlbwt: the text is too long (at most 2^48 - window - 1 symbols)");
    const uint8_t min_bytes = value_bytes(n + prm.window);
    if (prm.position_bytes == 0) prm.position_bytes = min_bytes;
    if (prm.position_bytes < min_bytes || prm.position_bytes > max_value_bytes)
        throw std::invalid_argument("pfp_rlbwt: invalid position_bytes");
}

/**
 * @brief computes the prefix-free parsing of the text, logs its sizes and builds the output from it (see build and
 *        build_level)
 * @param text the text T
 * @param prefix path prefix of the output files
 * @param prm parameters (checked; position_bytes is set if the length of T is known, else after the parse)
 * @param level_out if not nullptr, receives the runs and the suffix array (see build_level)
 * @return sizes of the parsing and of the output
 */
template <typename sym_t>
result build_parsed(text_ranges<sym_t>& text, const std::string& prefix, params prm, level_output* level_out)
{
    std::chrono::steady_clock::time_point time;
    log_phase_start(prm.log, time, current_level == 1 ? "building PFP" : "building PFP of " + level_text_name());
    text.split(prm.num_threads);
    const uint64_t seed = choose_trigger_seed<sym_t>(text, prm);
    /* on a level after the first, level_out->parsed decides after the parse whether the level continues (see parse),
       and a degenerate parse is not repeated (its prefix-free parsing does not halve the text) */
    bool stopped = false;
    std::function<bool(const result&)> parsed = [&](const result& sizes) {
        stopped = !level_out->parsed(sizes);
        return !stopped;
    };

    const bool decides = level_out != nullptr && level_out->parsed;
    parsing<sym_t> prs = parse<sym_t>(text, prm, seed, decides ? &parsed : nullptr);
    log_phase_end(prm.log, time);
    const uint64_t n = prs.n;

    if (!text.size_known()) {
        fit_text_length(n, prm);
        prs.position_bytes = prm.position_bytes;
    }

    // a degenerate parse: parse again with the other seed and keep the parse with more phrases
    if (level_out == nullptr && prs.num_phrases < n / (degenerate_parse_tolerance * prm.modulus)) {
        log_phase_start(prm.log, time, "building PFP again (only " + format_count(prs.num_phrases, "phrase") + ")");
        parsing<sym_t> retry = parse<sym_t>(text, prm, seed == 0 ? retry_trigger_seed : 0);
        log_phase_end(prm.log, time);
        if (retry.num_phrases > prs.num_phrases) prs = std::move(retry);
    }

    result res;
    res.n_with_terminator = n + 1;
    res.num_phrases = prs.num_phrases;
    res.num_phrase_runs = prs.num_phrase_runs;
    res.num_words = prs.dict.num_words;
    res.dictionary_size = prs.dict.size;

    if (prm.log) {
        /* |P| is the number of phrase runs; the share of the phrases that equal the phrase before them, i.e., that are
           merged into a run, follows if there are any (their number, if the share rounds to 0.00 %) */
        const uint64_t pfp_size = res.num_phrases + res.dictionary_size;
        const uint64_t merged = res.num_phrases - res.num_phrase_runs;
        const std::string merged_share = merged == 0 ? std::string()
            : 20000 * merged < res.num_phrases ? " (" + format_count(merged, "phrase") + " in runs)"
            : " (" + format_percentage(merged, res.num_phrases) + " runs)";
        log_message("|" + level_name("P") + "| = " + format_number(res.num_phrase_runs) + merged_share
            + ", |" + level_name("D") + "| = " + format_number(res.dictionary_size) + ", " + level_name("d") + " = "
            + format_number(res.num_words) + ", n/|PFP| = " + format_ratio(double(n) / double(pfp_size)) + "\n");
    }

    if (stopped) return res;

    // the largest length and alphabet size of the inputs of libsais choose the type of the suffix arrays
    const uint64_t max_index = std::max({ res.num_phrase_runs + 1, res.dictionary_size,
        dictionary_alphabet_size(prs.dict, prs.max_symbol) });

    if (max_index <= PFP_RLBWT_MAX_INDEX32) {
        res.r = build_from_parsing<sym_t, int32_t>(prs, prm, prefix, level_out);
    } else if (max_index <= PFP_RLBWT_MAX_INDEX40) {
        res.r = build_from_parsing<sym_t, sa_int40_t>(prs, prm, prefix, level_out);
    } else if (max_index <= uint64_t(SA_INT48_MAX_N)) {
        res.r = build_from_parsing<sym_t, sa_int48_t>(prs, prm, prefix, level_out);
    } else {
        throw std::length_error("pfp_rlbwt: the parsing is too large for 48-bit suffix arrays");
    }

    return res;
}

template <typename sym_t>
void build_level(const text_source<sym_t>& text, const std::string& prefix, const params& prm_in, level_output& out)
{
    params prm = prm_in;
    prm.position_bytes = value_bytes(text.size() + prm.window);
    source_ranges<sym_t> ranges(text);

    // current_level is the number of this level while it is built
    struct level_scope {
        level_scope() { current_level++; }
        ~level_scope() { current_level--; }
    } scope;

    out.r = build_parsed<sym_t>(ranges, prefix, prm, &out).r;
    out.completed = out.r != 0;
}

template <typename sym_t>
result build(text_ranges<sym_t>& text, const std::string& prefix, const params& prm_in)
{
    params prm = prm_in;
    check_params(prm);
    if (text.size_known()) fit_text_length(text.max_size(), prm);
    return build_parsed<sym_t>(text, prefix, prm, nullptr);
}

template <typename sym_t>
result build(const text_source<sym_t>& text, const std::string& prefix, const params& prm)
{
    source_ranges<sym_t> ranges(text);
    return build<sym_t>(ranges, prefix, prm);
}

} // namespace pfp_rlbwt

#endif // PFP_RLBWT_PFP_RLBWT_PFP_RLBWT_TPP
