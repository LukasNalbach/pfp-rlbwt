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
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <misc/files.hpp>
#include <misc/log.hpp>

#include <pfp_rlbwt/construction/merge.hpp>

namespace pfp_rlbwt {

// a reader of a file_source
template <typename sym_t>
class file_reader : public text_source<sym_t>::reader {
private:
    direct_reader file;
    uint64_t n; // n, the number of symbols in the file
    bool reverse; // true <=> the text is the reversed file
    const std::array<uint8_t, 256>* map; // the map of the bytes of the file_source (nullptr: not mapped)

public:
    file_reader(const std::string& path, uint64_t n, bool reverse, const std::array<uint8_t, 256>* map)
        : file(path), n(n), reverse(reverse), map(map)
    { }

    void read(uint64_t pos, uint64_t len, sym_t* dst) override
    {
        const uint64_t beg = reverse ? n - pos - len : pos;
        std::memcpy(dst, file.read(beg * sizeof(sym_t), (beg + len) * sizeof(sym_t)), len * sizeof(sym_t));
        if (reverse) std::reverse(dst, dst + len);
        if (map != nullptr) for (uint64_t i = 0; i < len; i++) dst[i] = (*map)[dst[i]];
    }
};

// a reader of a memory_source
template <typename sym_t>
class memory_reader : public text_source<sym_t>::reader {
private:
    const sym_t* text; // text[0..n-1]
    uint64_t n; // n, the number of symbols
    bool reverse; // true <=> the text is text[0..n-1] reversed

public:
    memory_reader(const sym_t* text, uint64_t n, bool reverse) : text(text), n(n), reverse(reverse) { }

    void read(uint64_t pos, uint64_t len, sym_t* dst) override
    {
        if (reverse) std::reverse_copy(text + n - pos - len, text + n - pos, dst);
        else std::memcpy(dst, text + pos, len * sizeof(sym_t));
    }
};

template <typename sym_t>
file_source<sym_t>::file_source(const std::string& path, bool reverse, const std::array<uint8_t, 256>* map)
    : path(path), reverse(reverse), mapped(map != nullptr)
{
    std::error_code ec;
    const uint64_t bytes = std::filesystem::file_size(path, ec);
    if (ec) throw std::runtime_error("pfp_rlbwt: cannot read " + path);
    if (bytes % sizeof(sym_t) != 0)
        throw std::invalid_argument("pfp_rlbwt: the file size is not a multiple of the symbol width");
    if (mapped && sizeof(sym_t) != 1)
        throw std::invalid_argument("pfp_rlbwt: a map is only applied to one-byte symbols");
    n = bytes / sizeof(sym_t);
    if (mapped) this->map = *map;
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

    return name;
}

/**
 * @brief builds the output from the parsing
 * @tparam idx_t type of the values of the suffix arrays of the parsing and of the dictionary (int32_t, sa_int40_t
 *         or int64_t)
 * @param prs the parsing; its phrases are freed
 * @param prm parameters
 * @param prefix path prefix of the output files
 * @return the number r of runs
 */
template <typename sym_t, typename idx_t>
uint64_t build_from_parsing(parsing<sym_t>& prs, const params& prm, const std::string& prefix)
{
    std::chrono::steady_clock::time_point time;
    const uint16_t p = prm.num_threads;

    log_phase_start(prm.log, time, "storing PFP (the ends of the phrase runs on disk)");
    phrase_arrays<idx_t, sym_t> arr = build_phrase_arrays<idx_t, sym_t>(prs, prefix, p);
    log_phase_end(prm.log, time);

    // the suffix array of the dictionary is built before the one of the parse and stays on disk
    dictionary_sa<idx_t> dsa = build_dictionary_sa<idx_t, sym_t>(prs.dict, prs.max_symbol, prs.w,
        prefix + ext_dictionary_sa, prm, time);

    parse_sa<idx_t> psa = build_parse_sa<idx_t, sym_t>(arr, prs.num_phrase_runs, prs.w, prs.dict,
        prefix + ext_parse_isa, prm, time);

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

        log_phase_start(prm.log, time, "sampling the PLCP array");
        const sampled_plcp plcp = build_sampled_plcp(run_starts, prs.n + 1, p);
        log_phase_end(prm.log, time);

        log_phase_start(prm.log, time, merge_phase_name(prm, false));
        merge<sym_t, idx_t>(prs, dsa, psa, arr, prm, prefix, false, nullptr, &plcp);
        log_phase_end(prm.log, time);
    } else {
        log_phase_start(prm.log, time, merge_phase_name(prm, true));
        r = merge<sym_t, idx_t>(prs, dsa, psa, arr, prm, prefix, true);
        log_phase_end(prm.log, time);
    }

    return r;
}

template <typename sym_t>
result build(const text_source<sym_t>& text, const std::string& prefix, const params& prm_in)
{
    // check the parameters and choose the defaults
    params prm = prm_in;
    const uint64_t n = text.size();
    if (prm.window == 0 || prm.modulus == 0 || prm.num_threads == 0)
        throw std::invalid_argument("pfp_rlbwt: invalid params");
    // the text positions and the window fit into max_value_bytes bytes: n + window < 2^(8 * max_value_bytes)
    constexpr uint64_t position_limit = uint64_t { 1 } << (8 * max_value_bytes);
    if (prm.window >= position_limit || n >= position_limit - prm.window)
        throw std::invalid_argument("pfp_rlbwt: the text is too long (at most 2^48 - window - 1 symbols)");
    const uint8_t min_bytes = value_bytes(n + prm.window);
    if (prm.position_bytes == 0) prm.position_bytes = min_bytes;
    if (prm.position_bytes < min_bytes || prm.position_bytes > max_value_bytes)
        throw std::invalid_argument("pfp_rlbwt: invalid position_bytes");
    if (prm.runs_per_block == 0 || prm.runs_per_block % byte_bits != 0)
        throw std::invalid_argument("pfp_rlbwt: runs_per_block is not a positive multiple of 8");
    if (prm.runs_per_symbol_sample == 0 || prm.runs_per_symbol_sample > prm.runs_per_block)
        throw std::invalid_argument("pfp_rlbwt: runs_per_symbol_sample is not in [1..runs_per_block]");

    for (uint64_t i = 0; i < prm.extra_sections.size(); i++) {
        const uint64_t id = prm.extra_sections[i].id;
        if (id < first_extra_section) throw std::invalid_argument("pfp_rlbwt: an extra section has an identifier "
            "smaller than first_extra_section");
        for (uint64_t j = 0; j < i; j++) if (prm.extra_sections[j].id == id)
            throw std::invalid_argument("pfp_rlbwt: two extra sections have the same identifier");
    }

    std::chrono::steady_clock::time_point time;
    log_phase_start(prm.log, time, "building PFP");
    const uint64_t seed = choose_trigger_seed<sym_t>(text, prm);
    parsing<sym_t> prs = parse<sym_t>(text, prm, seed);
    log_phase_end(prm.log, time);

    // a degenerate parse: parse again with the other seed and keep the parse with more phrases
    if (prs.num_phrases < n / (degenerate_parse_tolerance * prm.modulus)) {
        log_phase_start(prm.log, time, "only " + format_number(prs.num_phrases) + " phrases, parsing again with other "
            "trigger strings");
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
        const uint64_t pfp_size = res.num_phrases + res.dictionary_size;
        std::ostringstream ratio;
        ratio << std::fixed << std::setprecision(2) << double(n) / double(pfp_size);
        std::cout << "|P| = " << format_number(res.num_phrases) << " in " << format_number(res.num_phrase_runs)
                  << " phrase runs; sigma(P) = " << format_number(res.num_words)
                  << "; |D| = " << format_number(res.dictionary_size) << std::endl
                  << "|PFP| = |P| + |D| = " << format_number(pfp_size) << "; n / |PFP| = " << ratio.str() << std::endl;
    }

    // the largest length and alphabet size of the inputs of libsais choose the type of the suffix arrays
    const uint64_t max_index = std::max({ res.num_phrase_runs + 1, res.dictionary_size,
        dictionary_alphabet_size(prs.dict, prs.max_symbol) });

    if (max_index <= PFP_RLBWT_MAX_INDEX32) res.r = build_from_parsing<sym_t, int32_t>(prs, prm, prefix);
    else if (max_index <= PFP_RLBWT_MAX_INDEX40) res.r = build_from_parsing<sym_t, sa_int40_t>(prs, prm, prefix);
    else res.r = build_from_parsing<sym_t, int64_t>(prs, prm, prefix);

    return res;
}

} // namespace pfp_rlbwt

#endif // PFP_RLBWT_PFP_RLBWT_PFP_RLBWT_TPP
