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

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <misc/cli_help.hpp>
#include <misc/files.hpp>
#include <misc/log.hpp>
#include <misc/measurement.hpp>
#include <pfp_rlbwt/fasta.hpp>
#include <pfp_rlbwt/pfp_rlbwt.hpp>
#include <rlbwt/byte_file.hpp>

namespace {

// the options of pfp-rlbwt
struct pfp_rlbwt_cli_t {
    std::string path_input_file;
    std::string base_name; // -o: the base name of the output files; empty <=> <input_file>
    std::string text_name; // <name> of -m
    uint64_t p = omp_get_max_threads(); // -p
    uint64_t symbol_bytes = 1; // -symbol-bytes: bytes per symbol of the input file
    fasta_option fasta = fasta_option::detect; // -fasta
    bool fasta_mode = false; // whether the input file is read in the FASTA mode (see run)
    std::string fasta_alphabet = std::string(fasta_default_alphabet); // -A
    bool reverse = false; // -rev
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
        "writes it to <base_name>") + ext_runs + ", optionally with the suffix array values at the starts and at the "
        "ends of its runs and the LCP values at the starts of its runs; optionally also the BWT, the suffix array and "
        "the LCP array with one value per position.",
        "pfp-rlbwt [options] <input_file>");
    h.group("arguments");
    h.option("<input_file>", "", "the text; in the FASTA mode a FASTA file");
    h.group("input");
    h.option("-symbol-bytes <b>", "integer in [1, " + std::to_string(sizeof(uint64_t)) + "] (default: 1)",
        "the number of bytes per symbol of the input file (little endian); with more than 1, every symbol has to be "
        "at least " + min_symbol_str + "; with 1, if a byte below " + min_symbol_str + " occurs, the distinct bytes "
        "are mapped to " + min_symbol_str + ", " + std::to_string(min_symbol + 1) + ", ... in increasing order (to "
        "16-bit symbols if more than " + std::to_string(num_byte_values - min_symbol) + " distinct bytes occur), the "
        "BWT consists of the mapped bytes and the map is stored in <base_name>" + ext_runs + " (see the README)");
    h.option_fasta("FASTA mode (only with -symbol-bytes 1): replace each header line and each run of symbols outside "
        "the alphabet of -A (after upper-casing) with one separator symbol, skip the lines before the first header and "
        "drop carriage returns");
    h.option_fasta_alphabet();
    h.option("-rev", "", "compute the output for the reversed file");
    h.group("parsing");
    h.option("-w <w>", "integer >= 1 (default: " + std::to_string(defaults.window) + ")",
        "the length of a trigger string in symbols");
    h.option("-mod <m>", "integer >= 1 (default: " + std::to_string(defaults.modulus) + ")",
        "a window is a trigger string with probability 1/m");
    h.option("-pos-bytes <b>", "integer in [" + std::to_string(min_value_bytes) + ", " + std::to_string(max_value_bytes)
        + "] (default: the smallest number in which n + w fits)", "the number of bytes per text position in the arrays "
        "of the construction");
    h.option("-recursive <r>", "on | off (default: on)", "how the parse of the phrase runs is suffix-sorted\n"
        "on: with its prefix-free parsing, recursively, as long as a parse has at least <j> phrase runs (-rec-min), at "
        "most half as many phrase runs as its text has symbols, and a prefix-free parsing at most half as large as the "
        "parse\n"
        "off: with a suffix array construction");
    h.option("-rec-w <w>", "integer >= 1 (default: " + std::to_string(defaults.recursion_window) + ")",
        "the length of a trigger string in the prefix-free parsings of the parses (-recursive on)");
    h.option("-rec-mod <m>", "integer >= 1 (default: " + std::to_string(defaults.recursion_modulus) + ")",
        "a window of a parse is a trigger string with probability 1/m (-recursive on)");
    h.option("-rec-min <j>", "integer >= 0 (default: " + std::to_string(defaults.recursion_min_phrase_runs) + ")",
        "the smallest number of phrase runs of a parse that is suffix-sorted with its prefix-free parsing "
        "(-recursive on)");
    h.group("output");
    h.option("-o <base_name>", "a path (default: <input_file>)", "the base name of the output files");
    h.option("-format <f>", "plain | compressed (default: plain)", std::string("the format of <base_name>") + ext_runs
        + "\nplain: one record per run (its start, its symbol and its samples, each with a fixed number of bytes)\n"
        "compressed: the starts Elias-Fano coded and the symbols Huffman coded in blocks, the SA samples bit-packed "
        "and the LCP samples with RLE-GEF, each readable with O(1) disk reads");
    h.option("-sa", "", std::string("also write the SA values at the starts and at the ends of the runs to <base_name>")
        + ext_runs);
    h.option("-lcp", "", std::string("also write the LCP values at the starts of the runs to <base_name>") + ext_runs);
    h.option("-full-bwt", "", std::string("write <base_name>") + ext_bwt + " (the BWT, one symbol per position, with "
        "the bytes per symbol of the input or of the mapped bytes)");
    h.option("-full-sa", "", std::string("write <base_name>") + ext_sa + " (the suffix array, with the smallest number "
        "of bytes per value in which n + 1 fits)");
    h.option("-full-lcp", "", std::string("write <base_name>") + ext_lcp + " (the LCP array, with the smallest number "
        "of bytes per value in which n + 1 fits)");
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
            cli.base_name = args.value(s);
        } else if (s == "-p") {
            cli.p = args.integer(s, 1, omp_get_max_threads());
        } else if (s == "-symbol-bytes") {
            cli.symbol_bytes = args.integer(s, 1, sizeof(uint64_t));
        } else if (s == "-fasta") {
            cli.fasta = read_fasta_option(args);
        } else if (s == "-A") {
            cli.fasta_alphabet = args.value(s);
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
            else throw usage_error("-format expects plain or compressed");
        } else if (s == "-recursive") {
            const std::string recursive = args.value(s);
            if (recursive == "on") params.recursive = true;
            else if (recursive == "off") params.recursive = false;
            else throw usage_error("-recursive expects on or off");
        } else if (s == "-rec-w") {
            params.recursion_window = args.integer(s, 1);
        } else if (s == "-rec-mod") {
            params.recursion_modulus = args.integer(s, 1);
        } else if (s == "-rec-min") {
            params.recursion_min_phrase_runs = args.integer(s, 0);
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

    if (args.empty()) throw usage_error("expected <input_file>");
    cli.path_input_file = args.next();
    if (!args.empty()) throw usage_error("expected exactly one <input_file>");
    if (cli.base_name.empty()) cli.base_name = cli.path_input_file;
    if (cli.fasta == fasta_option::on && cli.symbol_bytes != 1) throw usage_error("-fasta on needs -symbol-bytes 1");
}

/* runs pfp_rlbwt::build on the text of the input file with one-byte symbols: a FASTA file is read once if it is not
   read reversed and no symbol of its text can be below min_symbol; else the file is scanned first and the bytes of the
   text are mapped if one is below min_symbol (see rlbwt::build_byte_file) */
pfp_rlbwt::result build_bytes()
{
    pfp_rlbwt::params& params = cli.params;
    auto build_text = [&](auto&& text) { return pfp_rlbwt::build(text, cli.base_name, params); };
    const pfp_rlbwt::fasta_transform transform(cli.fasta_alphabet); // the transformation of the FASTA mode
    // whether a symbol of the text of the FASTA mode can be below min_symbol: the separator or a symbol of the alphabet
    const bool small_symbol = transform.separator < pfp_rlbwt::min_symbol || std::any_of(transform.in_alphabet.begin(),
        transform.in_alphabet.begin() + pfp_rlbwt::min_symbol, [](uint8_t in) { return in != 0; });

    if (cli.fasta_mode && !cli.reverse && !small_symbol) {
        pfp_rlbwt::fasta_ranges<> text(cli.path_input_file, transform);
        return build_text(text);
    }

    return rlbwt::build_byte_file(cli.path_input_file, cli.fasta_mode ? &transform : nullptr, cli.reverse,
        pfp_rlbwt::min_symbol, params, cli.p, true, build_text);
}

// runs pfp_rlbwt::build on the input file read as symbols of type sym_t
template <typename sym_t>
pfp_rlbwt::result build()
{
    if constexpr (sizeof(sym_t) == 1) {
        return build_bytes();
    } else {
        const pfp_rlbwt::file_source<sym_t> text(cli.path_input_file, cli.reverse);
        return pfp_rlbwt::build(text, cli.base_name, cli.params);
    }
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
    cli.fasta_mode = cli.symbol_bytes == 1 && uses_fasta_mode(cli.fasta, cli.path_input_file);
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
        output_names += (output_names.empty() ? "" : ", ") + cli.base_name + ext;
        if (ext == pfp_rlbwt::ext_runs) output_names += samples;
    }

    log_run_start(std::string("building the RLBWT of") + (cli.reverse ? " the reversed" : ""), cli.path_input_file,
        bytes, p, "w = " + format_number(params.window) + ", m = " + format_number(params.modulus)
        + (cli.symbol_bytes == 1 ? std::string() : ", symbol-bytes = " + std::to_string(cli.symbol_bytes))
        + (cli.fasta_mode ? ", fasta = on" : "")
        + (params.recursive ? ", recursive with w = " + format_number(params.recursion_window) + ", m = "
            + format_number(params.recursion_modulus) : std::string()),
        output_names);

    const auto time_start = now();
    const uint64_t baseline_mem_usage = malloc_count_current();
    malloc_count_reset_peak();
    pfp_rlbwt::result result;

    switch (cli.symbol_bytes) {
        #define PFP_RLBWT_SYMBOL_CASE(sym_t) case sizeof(sym_t): result = build<sym_t>(); break;
        STRING_UTILS_FOR_EACH_SYMBOL_TYPE(PFP_RLBWT_SYMBOL_CASE)
        #undef PFP_RLBWT_SYMBOL_CASE
    }

    const uint64_t time_ns = time_diff_ns(time_start, now());
    // peak memory usage over all phases (see log_phase_start)
    const uint64_t peak = overall_peak() - baseline_mem_usage;
    uint64_t output_size = 0; // the total size of the output files in bytes

    for (const auto& [written, ext] : outputs) {
        if (!written) continue;
        output_size += std::filesystem::file_size(cli.base_name + ext, ec);
        if (ec) throw std::runtime_error("cannot read " + cli.base_name + ext);
    }

    const uint64_t n = result.n_with_terminator - 1;
    log_message("n = " + format_number(n) + ", r = " + format_number(result.r) + ", n/r = "
        + format_ratio(double(n) / double(result.r)) + "\n");
    log_construction_summary("RLBWT", time_ns, bytes, peak);
    log_output_size("output", output_size, bytes);

    if (cli.measurement_file.is_open()) {
        result_line(cli.measurement_file)
            ("algo", "build_pfp_rlbwt")
            ("text", cli.text_name)
            ("n", n)
            ("r", result.r)
            ("symbol_bytes", cli.symbol_bytes)
            ("fasta", cli.fasta_mode)
            ("reversed", cli.reverse)
            ("window", params.window)
            ("modulus", params.modulus)
            ("recursive", params.recursive)
            ("recursion_window", params.recursion_window)
            ("recursion_modulus", params.recursion_modulus)
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
