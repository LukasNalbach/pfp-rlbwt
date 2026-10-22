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
 * @file compress.tpp
 * @brief compress: the BWT construction with pfp-rlbwt, the choice of the start pairs, the coding of the runs and the
 *        FASTA mode
 */

#ifndef PFP_RLBWT_PFP_RLBWT_COMP_COMPRESS_TPP
#define PFP_RLBWT_PFP_RLBWT_COMP_COMPRESS_TPP

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <filesystem>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <coding/byte_codes.hpp>
#include <fasta/split.hpp>
#include <misc/bits.hpp>
#include <misc/files.hpp>
#include <misc/log.hpp>
#include <misc/packed_integers.hpp>
#include <misc/temp_files.hpp>
#include <misc/threads.hpp>
#include <pfp_rlbwt/fasta.hpp>
#include <pfp_rlbwt/output/readers.hpp>
#include <pfp_rlbwt/pfp_rlbwt.hpp>
#include <pfp_rlbwt/comp/pfp_rlbwt_comp.hpp>
#include <pfp_rlbwt/comp/run_model.hpp>

namespace pfp_rlbwt::comp {

/**
 * @brief the files pfp-rlbwt writes with a prefix; they are removed when the object is destroyed
 */
class pfp_output_files {
private:
    std::string prefix; // the prefix

public:
    explicit pfp_output_files(const std::string& prefix) : prefix(prefix) { }
    pfp_output_files(const pfp_output_files&) = delete;

    ~pfp_output_files()
    {
        std::error_code ec;
        std::filesystem::remove(prefix + pfp_rlbwt::ext_runs, ec);
    }
};

/**
 * @brief writes bytes to a file
 * @param path the file
 * @param bytes the bytes
 */
inline void write_file(const std::string& path, const std::vector<uint8_t>& bytes)
{
    direct_ofstream out(path);
    if (!out.is_open()) throw std::runtime_error("cannot create " + path);
    write_in_chunks(out, reinterpret_cast<const char*>(bytes.data()), bytes.size());
    out.close();
    if (!out.good()) throw std::runtime_error("cannot write " + path);
}

/**
 * @brief appends the run-length encoded BWT of a text T and the start pairs to the compressed file (see the file
 *        format in the README) and adds the sizes of its parts to the result
 * @param n the length of T
 * @param occurs the bytes that occur in T
 * @param make_source called as make_source(map), where map is a std::array of uint8_t or of uint16_t with [0..255]
 *        map[c] = the symbol of pfp-rlbwt of the byte c; returns T with its bytes mapped as a pfp_rlbwt::text_source
 * @param prefix the prefix of the temporary files of pfp-rlbwt
 * @param params parameters
 * @param out the bytes of the compressed file
 * @param res the result
 */
template <typename make_source_t>
inline void append_rlbwt(uint64_t n, const byte_set& occurs, make_source_t make_source, const std::string& prefix,
    const compress_params& params, std::vector<uint8_t>& out, compress_result& res)
{
    const uint64_t p = std::max<uint64_t>(1, params.num_threads);
    const uint64_t begin = out.size(); // the size of out before the run-length encoded BWT
    write_varint(out, n);

    if (n == 0) {
        res.header_bytes += out.size() - begin;
        return;
    }

    // [0..255] symbol[c] = the symbol of pfp-rlbwt of the byte c: its rank among the bytes of the text plus min_symbol
    std::array<uint16_t, num_byte_values> symbol {};
    uint64_t sigma = 0;

    for (uint64_t c = 0; c < num_byte_values; c++) {
        if (occurs[c]) symbol[c] = uint16_t(pfp_rlbwt::min_symbol + sigma++);
    }

    // the alphabet: one bit per byte value, set iff the byte occurs in the text
    for (uint64_t c = 0; c < num_byte_values; c += std::numeric_limits<uint8_t>::digits) {
        uint8_t bits = 0;
        for (uint64_t i = 0; i < std::numeric_limits<uint8_t>::digits; i++) bits |= uint8_t(occurs[c + i]) << i;
        out.push_back(bits);
    }

    // the BWT of T$, the SA values at the starts of its runs and the temporary files of pfp-rlbwt
    pfp_output_files files(prefix);
    pfp_rlbwt::params pfp_params;
    pfp_params.num_threads = uint16_t(std::min<uint64_t>(p, std::numeric_limits<uint16_t>::max()));
    pfp_params.write_sa_samples = true;
    pfp_params.log = params.log;
    pfp_rlbwt::result bwt;

    if (pfp_rlbwt::min_symbol + sigma <= num_byte_values) {
        std::array<uint8_t, num_byte_values> symbol8 {};
        for (uint64_t c = 0; c < num_byte_values; c++) symbol8[c] = uint8_t(symbol[c]);
        bwt = pfp_rlbwt::build(make_source(symbol8), prefix, pfp_params);
    } else {
        bwt = pfp_rlbwt::build(make_source(symbol), prefix, pfp_params);
    }

    auto time = now();
    log_phase_start(params.log, time, "coding the runs");

    /* [0, n) is split into k ranges of equal length; range g contains the text positions s with s * k / n = g; a
       start pair can have the text position s if 1 <= s <= n - 1 */
    const uint64_t k = std::min<uint64_t>(bwt.r, std::max<uint64_t>(1, params.max_start_pairs));
    auto range_of = [&](uint64_t s) { return uint64_t((uint128_t(s) * k) / n); };

    const uint64_t r = bwt.r;
    // the runs are coded in blocks of at most runs_per_block runs; block b contains the runs [b * r / num_blocks,
    // (b + 1) * r / num_blocks)
    const uint64_t num_blocks = div_ceil(r, std::max<uint64_t>(1, params.runs_per_block));
    auto block_start = [&](uint64_t b) { return uint64_t((uint128_t(b) * r) / num_blocks); };
    const std::string run_data_path = prefix + pfp_rlbwt::ext_runs;
    pfp_rlbwt::run_data_file run_data;

    if (!run_data.open(run_data_path) || !run_data.has_sa_samples()) {
        throw std::runtime_error("cannot read " + run_data_path);
    }

    const uint64_t threads = std::min(p, num_blocks);
    /* [t * k + g] smallest = the smallest SA value in range g at a run start of a block coded by thread t (n, if
       there is none), run_of = its run */
    std::vector<uint64_t> smallest(threads * k, n);
    no_init_vector<uint64_t> run_of;
    resize_no_init(run_of, threads * k);
    std::atomic<uint64_t> terminator_run = 0; // the run of $
    // [0..num_blocks-1] the coded run symbols and run lengths of each block
    std::vector<std::vector<uint8_t>> symbol_bytes(num_blocks), length_bytes(num_blocks);
    std::atomic<uint64_t> next_block = 0;

    run_on_threads(threads, [&](uint64_t t, uint64_t) {
        uint64_t* thread_smallest = smallest.data() + t * k;
        uint64_t* thread_run_of = run_of.data() + t * k;

        for (uint64_t b = next_block++; b < num_blocks; b = next_block++) {
            const uint64_t first = block_start(b), end = block_start(b + 1);
            pfp_rlbwt::run_stream runs(run_data, first, end);
            pfp_rlbwt::run_record run;
            run_model model(sigma, n);
            range_encoder symbols, lengths;

            for (uint64_t x = first; x < end; x++) {
                runs.next(run);
                const uint64_t c = run.symbol, len = run.length, s = run.sa_start;

                if (c == pfp_rlbwt::terminator) {
                    terminator_run = x;
                    model.skip_terminator();
                } else {
                    model.encode(symbols, lengths, c - pfp_rlbwt::min_symbol, len);
                }

                if (s != 0 && s < n && s < thread_smallest[range_of(s)]) {
                    thread_smallest[range_of(s)] = s;
                    thread_run_of[range_of(s)] = x;
                }
            }

            symbol_bytes[b] = symbols.finish();
            length_bytes[b] = lengths.finish();
        }
    });

    // smallest[g] and run_of[g] become the smallest SA value in range g at a run start and its run
    for (uint64_t t = 1; t < threads; t++) {
        for (uint64_t g = 0; g < k; g++) {
            if (smallest[t * k + g] < smallest[g]) {
                smallest[g] = smallest[t * k + g];
                run_of[g] = run_of[t * k + g];
            }
        }
    }

    // the start pairs (run, SA value at its start), sorted by the SA value
    std::vector<std::pair<uint64_t, uint64_t>> pairs;

    for (uint64_t g = 0; g < k; g++) {
        if (smallest[g] < n) pairs.emplace_back(run_of[g], smallest[g]);
    }

    write_varint(out, r);
    write_varint(out, terminator_run);
    write_varint(out, pairs.size());
    const uint64_t pairs_begin = out.size(); // the size of out before the start pairs

    if (!pairs.empty()) {
        // the text positions are stored as differences to the previous one (the first one to 0)
        uint64_t max_difference = 0;

        for (uint64_t j = 0; j < pairs.size(); j++) {
            max_difference = std::max(max_difference, pairs[j].second - (j == 0 ? 0 : pairs[j - 1].second));
        }

        const uint8_t run_width = std::bit_width(r - 1);
        const uint8_t difference_width = std::bit_width(max_difference);
        out.push_back(difference_width);
        bit_packer packer(out);

        for (uint64_t j = 0; j < pairs.size(); j++) {
            packer.push(pairs[j].first, run_width);
            packer.push(pairs[j].second - (j == 0 ? 0 : pairs[j - 1].second), difference_width);
        }

        packer.finish();
    }

    res.start_pair_bytes = out.size() - pairs_begin;
    write_varint(out, num_blocks);

    for (uint64_t b = 0; b < num_blocks; b++) {
        write_varint(out, symbol_bytes[b].size());
        write_varint(out, length_bytes[b].size());
    }

    res.header_bytes += out.size() - begin - res.start_pair_bytes;

    for (uint64_t b = 0; b < num_blocks; b++) {
        out.insert(out.end(), symbol_bytes[b].begin(), symbol_bytes[b].end());
        out.insert(out.end(), length_bytes[b].begin(), length_bytes[b].end());
        res.symbol_bytes += symbol_bytes[b].size();
        res.length_bytes += length_bytes[b].size();
    }

    log_phase_end(params.log, time);

    res.sigma = sigma;
    res.r = r;
    res.num_run_blocks = num_blocks;
    res.num_start_pairs = pairs.size();
}

inline compress_result compress(const byte_text& text, const std::string& output_file, const compress_params& params)
{
    compress_result res;
    res.n = res.sequence_length = text.size();
    const uint64_t p = std::max<uint64_t>(1, params.num_threads);
    std::vector<uint8_t> out;
    out.push_back(params.fasta ? text_mode_fasta : text_mode_plain);
    res.header_bytes = out.size();
    const std::string prefix = temporary_file_path("pfp-rlbwt-comp-", output_file);

    if (params.fasta) {
        // the FASTA file is scanned once; pfp-rlbwt reads the sequence text from it
        const fasta_parts parts = text.bytes() != nullptr
            ? split_fasta([&] { return pfp_rlbwt::memory_file(text.bytes(), text.size()); }, p, params.log)
            : split_fasta([&] { return direct_reader(text.file()); }, p, params.log);

        parts.section.write(out);
        res.num_records = parts.section.num_records;
        res.sequence_length = parts.text.text_size;
        res.fasta_header_bytes = parts.section.headers.size();
        res.fasta_layout_bytes = parts.section.layout.size();
        res.header_bytes = out.size() - res.fasta_header_bytes - res.fasta_layout_bytes;
        byte_set occurs {};
        for (uint64_t c = 0; c < num_byte_values; c++) occurs[c] = parts.text.freq[c] > 0;
        const fasta_sample* samples = parts.text.samples.data();
        const uint64_t num_samples = parts.text.samples.size();

        append_rlbwt(parts.text.text_size, occurs, [&](const auto& map) {
            using sym_t = typename std::decay_t<decltype(map)>::value_type;

            if (text.bytes() != nullptr) {
                return pfp_rlbwt::fasta_source<sym_t>(text.bytes(), text.size(), parts.text.text_size, parts.transform,
                    samples, num_samples, false, &map);
            }

            return pfp_rlbwt::fasta_source<sym_t>(text.file(), parts.text.text_size, parts.transform, samples,
                num_samples, false, &map);
        }, prefix, params, out, res);
    } else {
        auto time = now();
        log_phase_start(params.log, time, "reading the alphabet");
        const byte_set occurs = text.occurring_bytes(p);
        log_phase_end(params.log, time);

        append_rlbwt(text.size(), occurs, [&](const auto& map) {
            using sym_t = typename std::decay_t<decltype(map)>::value_type;
            return mapped_byte_source<sym_t>(text, map);
        }, prefix, params, out, res);
    }

    write_file(output_file, out);
    res.file_bytes = out.size();
    return res;
}

} // namespace pfp_rlbwt::comp

#endif // PFP_RLBWT_PFP_RLBWT_COMP_COMPRESS_TPP
