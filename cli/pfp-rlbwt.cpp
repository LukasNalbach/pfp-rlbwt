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
 * @file pfp-rlbwt.cpp
 * @brief computes the run-length encoded BWT of a file and optionally its samples with the prefix-free parsing
 *        (pfp-rlbwt)
 */

#include <omp.h>

#include <array>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <misc/cli_help.hpp>
#include <misc/files.hpp>
#include <misc/log.hpp>
#include <misc/measurement.hpp>
#include <pfp_rlbwt/pfp_rlbwt.hpp>

namespace {

// the options of pfp-rlbwt
struct pfp_rlbwt_cli_t {
    std::string path_input_file;
    std::string prefix; // -o: prefix of the output files; empty <=> <input_file>
    std::string text_name; // <name> of -m
    uint64_t p = omp_get_max_threads(); // -p
    uint64_t symbol_bytes = 1; // -symbol-bytes: bytes per symbol of the input file
    bool reverse = false; // -rev
    bool map_bytes = false; // -map
    pfp_rlbwt::params params; // the parameters of the construction
    direct_ofstream measurement_file; // -m
} cli;

/**
 * @brief prints the help of pfp-rlbwt
 * @param out the stream the help is written to
 */
void print_help(std::ostream& out)
{
    using namespace pfp_rlbwt;
    const params defaults;
    const std::string min_symbol_str = std::to_string(min_symbol);
    help_printer h(out);
    h.intro(std::string("pfp-rlbwt: computes the run-length encoded BWT of a file with the prefix-free parsing and "
        "writes it to <prefix>") + ext_runs + ", optionally with the suffix array values at the starts and at the "
        "ends of its runs and the LCP values at the starts of its runs; optionally also the BWT, the suffix array and "
        "the LCP array with one value per position.",
        "pfp-rlbwt [options] <input_file>");
    h.option("<input_file>", "a path", "the text");
    h.group("input");
    h.option("-symbol-bytes <b>", "integer in [1, " + std::to_string(sizeof(sym48_t)) + "] (default: 1)",
        "the number of bytes per symbol of the input file (little endian); every symbol has to be at least "
        + min_symbol_str);
    h.option("-map", "", "map the bytes of the input file to " + min_symbol_str + ", " + std::to_string(min_symbol + 1)
        + ", ... in increasing order (for a file that contains smaller bytes; only with -symbol-bytes 1); the BWT "
        "consists of the mapped bytes");
    h.option("-rev", "", "compute the output for the reversed file");
    h.group("parsing");
    h.option("-w <w>", "integer >= 1 (default: " + std::to_string(defaults.window) + ")",
        "the length of a trigger string in symbols");
    h.option("-mod <m>", "integer >= 1 (default: " + std::to_string(defaults.modulus) + ")",
        "a window is a trigger string with probability 1/m");
    h.option("-pos-bytes <b>", "integer in [" + std::to_string(min_value_bytes) + ", " + std::to_string(max_value_bytes)
        + "] (default: the fewest that hold n + w)", "the number of bytes per text position in the arrays of the "
        "construction");
    h.group("output");
    h.option("-o <prefix>", "a path (default: <input_file>)", "the prefix of the output files");
    h.option("-format <f>", "plain or compressed (default: plain)", std::string("the format of <prefix>") + ext_runs
        + ": plain stores one record per run (its start, its symbol and its samples, each with a fixed number of "
        "bytes); compressed stores the starts Elias-Fano coded and the symbols Huffman coded in blocks, the SA "
        "samples bit-packed and the LCP samples with RLE-GEF, each readable with O(1) disk reads");
    h.option("-sa", "", std::string("also write the SA values at the starts and at the ends of the runs to <prefix>")
        + ext_runs);
    h.option("-lcp", "", std::string("also write the LCP values at the starts of the runs to <prefix>") + ext_runs);
    h.option("-full-bwt", "", std::string("write <prefix>") + ext_bwt + " (the BWT, one symbol of the input width "
        "per position)");
    h.option("-full-sa", "", std::string("write <prefix>") + ext_sa + " (the suffix array, the fewest bytes per "
        "value that hold n + 1)");
    h.option("-full-lcp", "", std::string("write <prefix>") + ext_lcp + " (the LCP array, the fewest bytes per value "
        "that hold n + 1)");
    h.group("execution");
    h.option_threads();
    h.group("measurement");
    h.option_measurement();
}

/**
 * @brief reads the arguments of pfp-rlbwt into cli
 * @param argc the number of command-line arguments
 * @param argv the command-line arguments
 */
void read_arguments(int argc, char** argv)
{
    pfp_rlbwt::params& params = cli.params;
    cli_arguments args(argc, argv);

    while (args.at_option()) {
        const std::string s = args.next();

        if (s == "-o") {
            cli.prefix = args.value(s);
        } else if (s == "-p") {
            cli.p = args.integer(s, 1, omp_get_max_threads());
        } else if (s == "-symbol-bytes") {
            cli.symbol_bytes = args.integer(s, 1, sizeof(pfp_rlbwt::sym48_t));
        } else if (s == "-map") {
            cli.map_bytes = true;
        } else if (s == "-rev") {
            cli.reverse = true;
        } else if (s == "-w") {
            params.window = args.integer(s, 1);
        } else if (s == "-mod") {
            params.modulus = args.integer(s, 1);
        } else if (s == "-format") {
            const std::string format = args.value(s);
            if (format == "plain") params.format = pfp_rlbwt::file_format::plain;
            else if (format == "compressed") params.format = pfp_rlbwt::file_format::compressed;
            else throw usage_error("-format must be plain or compressed");
        } else if (s == "-pos-bytes") {
            params.position_bytes = args.integer(s, pfp_rlbwt::min_value_bytes, pfp_rlbwt::max_value_bytes);
        } else if (s == "-sa") {
            params.write_sa_samples = true;
        } else if (s == "-lcp") {
            params.write_lcp_samples = true;
        } else if (s == "-full-bwt") {
            params.write_bwt = true;
        } else if (s == "-full-sa") {
            params.write_sa = true;
        } else if (s == "-full-lcp") {
            params.write_lcp = true;
        } else if (s == "-m") {
            read_measurement_option(args, cli.measurement_file, cli.text_name);
        } else {
            throw unknown_option(s);
        }
    }

    if (args.empty()) throw usage_error("missing <input_file>");
    cli.path_input_file = args.next();
    if (!args.empty()) throw usage_error("expected exactly one <input_file>");
    if (cli.prefix.empty()) cli.prefix = cli.path_input_file;
    if (cli.map_bytes && cli.symbol_bytes != 1) throw usage_error("-map requires -symbol-bytes 1");
}

// returns the map of the bytes occurring in the input file to min_symbol, min_symbol + 1, ... in increasing order
std::array<uint8_t, 256> byte_map()
{
    const pfp_rlbwt::file_source<uint8_t> text(cli.path_input_file);
    const uint64_t n = text.size();
    const uint64_t p = cli.p;
    const uint64_t blocks = (n + pfp_rlbwt::read_block_size - 1) / pfp_rlbwt::read_block_size;
    // the bytes each thread has read: [0..p-1][0..255] seen[i_p][c] = 1 <=> thread i_p has read the byte c
    std::vector<std::array<uint8_t, 256>> seen(p);

    #pragma omp parallel num_threads(p)
    {
        std::array<uint8_t, 256>& seen_t = seen[omp_get_thread_num()];
        const std::unique_ptr<pfp_rlbwt::text_source<uint8_t>::reader> reader = text.new_reader();
        std::vector<uint8_t> buffer(pfp_rlbwt::read_block_size);

        #pragma omp for schedule(dynamic)
        for (uint64_t b = 0; b < blocks; b++) {
            const uint64_t pos = b * pfp_rlbwt::read_block_size;
            const uint64_t len = std::min(pfp_rlbwt::read_block_size, n - pos);
            reader->read(pos, len, buffer.data());
            for (uint64_t i = 0; i < len; i++) seen_t[buffer[i]] = 1;
        }
    }

    /* the symbols the bytes are mapped to: [0..255] map[c] is the symbol the byte c is mapped to (min_symbol,
       min_symbol + 1, ... in increasing order of the occurring bytes; 0 for bytes that do not occur) */
    std::array<uint8_t, 256> map {};
    uint64_t next = pfp_rlbwt::min_symbol;

    for (uint64_t c = 0; c < 256; c++) {
        bool occurs = false;
        for (const std::array<uint8_t, 256>& seen_t : seen) occurs |= seen_t[c] != 0;
        if (!occurs) continue;

        if (next > UINT8_MAX) {
            throw std::runtime_error("the input file contains more than " + std::to_string(256 - pfp_rlbwt::min_symbol)
                + " distinct bytes, so they cannot be mapped");
        }

        map[c] = next++;
    }

    return map;
}

// runs pfp_rlbwt::build on the input file read as symbols of type sym_t
template <typename sym_t>
pfp_rlbwt::result build(const std::array<uint8_t, 256>* map)
{
    const pfp_rlbwt::file_source<sym_t> text(cli.path_input_file, cli.reverse, map);
    return pfp_rlbwt::build(text, cli.prefix, cli.params);
}

// computes the output files of the input file
void run()
{
    pfp_rlbwt::params& params = cli.params;
    const uint64_t p = cli.p;
    std::error_code ec;
    const uint64_t bytes = std::filesystem::file_size(cli.path_input_file, ec);
    if (ec) throw std::runtime_error("cannot read " + cli.path_input_file);
    if (bytes % cli.symbol_bytes != 0) {
        throw std::runtime_error("the size of the input file is not a multiple of " + std::to_string(cli.symbol_bytes)
            + " bytes");
    }
    params.num_threads = p;
    params.log = true;
    omp_set_num_threads(p);

    // the output files, each with whether it is written
    const std::vector<std::pair<bool, const char*>> outputs = {
        { true, pfp_rlbwt::ext_runs }, { params.write_bwt, pfp_rlbwt::ext_bwt },
        { params.write_sa, pfp_rlbwt::ext_sa }, { params.write_lcp, pfp_rlbwt::ext_lcp }
    };
    // the format of the runs file and its samples
    const std::string format = params.format == pfp_rlbwt::file_format::plain ? "plain" : "compressed";
    const std::string samples = " (" + format + (params.write_sa_samples && params.write_lcp_samples
        ? ", with SA- and LCP-samples)" : params.write_sa_samples ? ", with SA-samples)" : params.write_lcp_samples
        ? ", with LCP-samples)" : ")");

    std::string output_names; // the written output files, separated by ", "

    for (const auto& [written, ext] : outputs) {
        if (!written) continue;
        output_names += (output_names.empty() ? "" : ", ") + cli.prefix + ext;
        if (ext == pfp_rlbwt::ext_runs) output_names += samples;
    }

    std::cout << "building the BWT of " << (cli.reverse ? "the reversed " : "") << cli.path_input_file << " ("
              << format_size(bytes) << ") with the prefix-free parsing using " << format_threads(p)
              << ", window length w = " << params.window << ", modulus " << params.modulus
              << (cli.symbol_bytes == 1 ? std::string() : ", " + std::to_string(cli.symbol_bytes) + " bytes per symbol")
              << std::endl;
    std::cout << "output: " << output_names << std::endl << std::endl;

    const auto time_start = now();
    const uint64_t baseline_mem_usage = malloc_count_current();
    malloc_count_reset_peak();
    std::array<uint8_t, 256> map {};

    if (cli.map_bytes) {
        auto time = now();
        log_phase_start(true, time, "mapping the bytes");
        map = byte_map();
        log_phase_end(true, time);
    }

    pfp_rlbwt::result result;

    switch (cli.symbol_bytes) {
        case 1: result = build<uint8_t>(cli.map_bytes ? &map : nullptr); break;
        case 2: result = build<uint16_t>(nullptr); break;
        case 3: result = build<pfp_rlbwt::sym24_t>(nullptr); break;
        case 4: result = build<uint32_t>(nullptr); break;
        case 5: result = build<pfp_rlbwt::sym40_t>(nullptr); break;
        default: result = build<pfp_rlbwt::sym48_t>(nullptr); break;
    }

    const uint64_t time_ns = time_diff_ns(time_start, now());
    // peak memory usage over all phases (see log_phase_start)
    const uint64_t peak = overall_peak() - baseline_mem_usage;
    uint64_t output_size = 0; // the total size of the output files in bytes
    std::vector<std::pair<std::string, uint64_t>> output_sizes; // the written output files with their sizes in bytes

    for (const auto& [written, ext] : outputs) {
        if (!written) continue;
        const uint64_t size = std::filesystem::file_size(cli.prefix + ext, ec);
        if (ec) throw std::runtime_error("cannot read " + cli.prefix + ext);
        output_sizes.emplace_back(cli.prefix + ext, size);
        output_size += size;
    }

    std::cout << std::endl << "n = " << format_number(result.n_with_terminator) << " (with the terminator), r = "
              << format_number(result.r) << ", n/r = "
              << format_ratio(double(result.n_with_terminator) / double(result.r)) << std::endl << std::endl;

    if (begin_listing(output_size)) {
        for (const auto& [path, size] : output_sizes) log_component(path, size);
        end_listing();
    }

    std::cout << std::endl;
    log_construction_summary(time_ns, bytes, peak);
    log_output_size("output", output_size, bytes);

    if (cli.measurement_file.is_open()) {
        result_line(cli.measurement_file)
            ("type", "pfp")
            ("text", cli.text_name)
            ("n", result.n_with_terminator)
            ("r", result.r)
            ("symbol_bytes", cli.symbol_bytes)
            ("reverse", cli.reverse)
            ("window", params.window)
            ("modulus", params.modulus)
            ("num_threads", p)
            ("phrases", result.num_phrases)
            ("phrase_runs", result.num_phrase_runs)
            ("words", result.num_words)
            ("dictionary_size", result.dictionary_size)
            ("format", format)
            ("sa_samples", params.write_sa_samples)
            ("lcp_samples", params.write_lcp_samples)
            ("full_bwt", params.write_bwt)
            ("full_sa", params.write_sa)
            ("full_lcp", params.write_lcp)
            ("size_output", output_size)
            ("time_build", time_ns)
            ("throughput_mb_per_s", throughput_in_mb_per_s(bytes, time_ns))
            ("peak_memory_usage", peak);
        cli.measurement_file.close();
    }
}

} // namespace

int main(int argc, char** argv)
{
    return run_cli(argc, argv, print_help, [&] {
        read_arguments(argc, argv);
        run();
    });
}
