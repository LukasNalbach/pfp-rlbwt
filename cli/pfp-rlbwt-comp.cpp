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
 * @file pfp-rlbwt-comp.cpp
 * @brief compresses a text into its run-length encoded BWT (see pfp_rlbwt::comp::compress) and restores it with -d
 *        (see pfp_rlbwt::comp::decompress)
 */

#include <omp.h>

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

#include <misc/compressor_cli.hpp>
#include <misc/log.hpp>
#include <misc/measurement.hpp>
#include <pfp_rlbwt/comp/pfp_rlbwt_comp.hpp>

namespace {

// number of significant digits of the printed bits per input symbol
constexpr int bits_per_symbol_digits = 4;

constexpr compressor_spec spec {
    .name = "pfp-rlbwt-comp",
    .ending = "rlbwt",
    .description = "compresses a text into its run-length encoded BWT; -d restores it",
    .fasta = true
};

// compresses the input file into the output file
void compress(compressor_cli_t& cli)
{
    const pfp_rlbwt::comp::byte_text text(cli.path_input_file);
    const uint64_t n = text.size();
    log_out() << "compressing " << cli.path_input_file << (cli.fasta ? " as a FASTA file" : "") << " ("
              << format_size(n) << ") using " << format_threads(cli.threads) << std::endl;
    log_out() << "output: " << cli.path_output_file << std::endl << std::endl;

    const uint64_t baseline = malloc_count_current();
    malloc_count_reset_peak();
    const auto time_start = now();
    const pfp_rlbwt::comp::compress_result res =
        pfp_rlbwt::comp::compress(text, cli.path_output_file,
        { .num_threads = cli.threads, .fasta = cli.fasta, .log = cli.log });
    const uint64_t time_ns = time_diff_ns(time_start, now());
    const uint64_t peak = overall_peak() - baseline;

    log_out() << std::endl;

    if (cli.fasta) {
        log_out() << format_number(res.num_records) << " records, sequence length "
                  << format_number(res.sequence_length) << std::endl;
    }

    log_out() << "n = " << format_number(res.sequence_length) << ", sigma = " << res.sigma << ", r = "
              << format_number(res.r) << ", n/r = "
              << format_ratio(double(res.sequence_length) / double(std::max<uint64_t>(1, res.r))) << ", "
              << format_number(res.num_start_pairs) << " start pairs" << std::endl;
    log_out() << "header: " << format_size(res.header_bytes) << ", start pairs: " << format_size(res.start_pair_bytes)
              << ", run symbols: " << format_size(res.symbol_bytes) << ", run lengths: "
              << format_size(res.length_bytes);
    if (cli.fasta) {
        log_out() << ", headers: " << format_size(res.fasta_header_bytes) << ", line layout: "
                  << format_size(res.fasta_layout_bytes);
    }
    log_out() << std::endl;
    log_out() << "compression time: " << format_time(time_ns) << " (" << format_throughput(n, time_ns) << ")"
              << std::endl;
    if constexpr (is_malloc_count_linked) log_peak_memory_usage("peak memory usage", peak, n);
    log_compression_ratio(res.file_bytes, n);
    log_out() << "bits per input symbol: " << std::setprecision(bits_per_symbol_digits)
              << double(std::numeric_limits<uint8_t>::digits * res.file_bytes) / double(std::max<uint64_t>(1, n))
              << std::endl;

    if (cli.measurement_file.is_open()) {
        result_line(cli.measurement_file)
            ("type", "pfp-rlbwt-compress")
            ("text", cli.text_name)
            ("n", res.n)
            ("fasta", cli.fasta)
            ("num_records", res.num_records)
            ("sequence_length", res.sequence_length)
            ("sigma", res.sigma)
            ("r", res.r)
            ("num_threads", cli.threads)
            ("num_start_pairs", res.num_start_pairs)
            ("num_run_blocks", res.num_run_blocks)
            ("size_header", res.header_bytes)
            ("size_start_pairs", res.start_pair_bytes)
            ("size_symbols", res.symbol_bytes)
            ("size_lengths", res.length_bytes)
            ("size_fasta_headers", res.fasta_header_bytes)
            ("size_fasta_layout", res.fasta_layout_bytes)
            ("size_file", res.file_bytes)
            ("time_compress", time_ns)
            ("peak_memory_usage", peak);
    }
}

// restores the text of the input file into the output file
void decompress(compressor_cli_t& cli)
{
    log_out() << "decompressing " << cli.path_input_file << " using " << format_threads(cli.threads) << std::endl;
    log_out() << "output: " << cli.path_output_file << std::endl << std::endl;

    const uint64_t baseline = malloc_count_current();
    malloc_count_reset_peak();
    const auto time_start = now();
    const uint64_t n = pfp_rlbwt::comp::decompress(cli.path_input_file, cli.path_output_file,
        { .num_threads = cli.threads, .log = cli.log });
    const uint64_t time_ns = time_diff_ns(time_start, now());
    const uint64_t peak = overall_peak() - baseline;

    log_out() << std::endl << "decompression time: " << format_time(time_ns) << " ("
              << format_throughput(n, time_ns) << ")" << std::endl;
    if constexpr (is_malloc_count_linked) log_peak_memory_usage("peak memory usage", peak, n);

    if (cli.measurement_file.is_open()) {
        result_line(cli.measurement_file)
            ("type", "pfp-rlbwt-decompress")
            ("text", cli.text_name)
            ("n", n)
            ("num_threads", cli.threads)
            ("time_decompress", time_ns)
            ("peak_memory_usage", peak);
    }
}

} // namespace

int main(int argc, char** argv)
{
    return run_compressor_cli(argc, argv, spec, [](help_printer&) {},
        [](const std::string&, cli_arguments&) { return false; }, compress, decompress);
}
