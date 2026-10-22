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

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <iterator>
#include <limits>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <fasta/split.hpp>
#include <misc/temp_files.hpp>
#include <pfp_rlbwt/comp/pfp_rlbwt_comp.hpp>
#include "test-fuzz.hpp"
#include "test-inputs.hpp"

namespace {

std::mt19937_64 gen(test_seed());

// maximum length of a random text; adjusted by run_fuzz
uint64_t max_input_size = 100000;
// smallest value run_fuzz reduces max_input_size to
constexpr uint64_t min_input_size = 64;
// the numbers of threads of the decompression that every round trip checks
const std::vector<uint64_t> decompression_threads = { 1, 2, 3, 8 };
// the largest number of start pairs drawn by verify_random_texts (besides the default)
constexpr uint64_t max_drawn_start_pairs = 16;
// the largest number of blocks of runs drawn by draw_runs_per_block (besides the default)
constexpr uint64_t max_drawn_run_blocks = 16;
// the numbers of threads of the compression and the decompression of the FASTA edge cases
const std::vector<uint64_t> fasta_threads = { 1, 3, 8 };
// parameters of the random FASTA files: the largest number of records, of sequences the records copy, of mutations
// of a copy, of the length of an N run or a lowercase run, and of the line width
constexpr uint64_t max_fasta_records = 64;
constexpr uint64_t max_fasta_sources = 4;
constexpr uint64_t max_mutations = 8;
constexpr uint64_t max_run = 100;
constexpr uint64_t max_line_width = 100;
// parameters of the random headers: the largest number of tokens, difference of a number to the number of the
// previous header, number of leading zeros of a number and length of a word
constexpr uint64_t max_header_tokens = 12;
constexpr uint64_t max_header_delta = 3;
constexpr uint64_t max_leading_zeros = 3;
constexpr uint64_t max_header_word = 5;

// returns a random integer in [min, max]
uint64_t random_in(uint64_t min, uint64_t max) { return std::uniform_int_distribution<uint64_t>(min, max)(gen); }

// a temporary file that is removed when the object is destroyed
struct temporary_file {
    std::string path = temporary_file_path("test-pfp-rlbwt-comp-");
    ~temporary_file() { std::error_code ec; std::filesystem::remove(path, ec); }
};

/**
 * @brief compresses a text from a file and from memory, and checks that both files are equal and that decompressing
 *        them into memory and into a file with each number of threads of decompression_threads restores the text
 * @param T the text
 * @param params the parameters of the compression
 */
void verify_round_trip(const std::string& T, const pfp_rlbwt::comp::compress_params& params)
{
    temporary_file input, compressed, compressed_from_memory, output;
    std::ofstream(input.path, std::ios::binary).write(T.data(), std::streamsize(T.size()));

    const pfp_rlbwt::comp::compress_result res = pfp_rlbwt::comp::compress(input.path, compressed.path, params);
    ASSERT_EQ(res.n, T.size());
    ASSERT_EQ(res.file_bytes, std::filesystem::file_size(compressed.path));
    ASSERT_EQ(res.file_bytes, res.header_bytes + res.start_pair_bytes + res.symbol_bytes + res.length_bytes +
        res.fasta_header_bytes + res.fasta_layout_bytes);
    ASSERT_LE(res.num_start_pairs, std::min(res.r, params.max_start_pairs));

    const pfp_rlbwt::comp::byte_text in_memory(reinterpret_cast<const uint8_t*>(T.data()), T.size());
    pfp_rlbwt::comp::compress(in_memory, compressed_from_memory.path, params);
    std::ifstream a(compressed.path, std::ios::binary), b(compressed_from_memory.path, std::ios::binary);
    ASSERT_EQ(std::string(std::istreambuf_iterator<char>(a), {}), std::string(std::istreambuf_iterator<char>(b), {}));

    for (uint64_t p : decompression_threads) {
        const uint64_t k = random_in(1, pfp_rlbwt::comp::default_num_strands);
        ASSERT_EQ(pfp_rlbwt::comp::decompress(compressed.path, { .num_threads = p, .num_strands = k }), T)
            << p << " threads, " << k << " strands";
        ASSERT_EQ(pfp_rlbwt::comp::decompress(compressed.path, output.path, { .num_threads = p, .num_strands = k }),
            T.size());
        std::ifstream restored(output.path, std::ios::binary);
        ASSERT_EQ(std::string(std::istreambuf_iterator<char>(restored), {}), T)
            << p << " threads, " << k << " strands, file";
    }
}

/* checks that the parts of T as a FASTA file do not depend on how its scan is split into chunks: chunks of a random
   minimum size, which can start inside header and sequence lines, against one chunk */
void verify_fasta_chunks(const std::string& T)
{
    auto open = [&] { return memory_file(reinterpret_cast<const uint8_t*>(T.data()), T.size()); };
    const fasta_parts one = split_fasta(open, 1, false, T.size() + 1);
    const fasta_parts split = split_fasta(open,
        random_in(1, std::max<uint64_t>(1, omp_get_max_threads())), false, random_in(1, T.size() + 1));
    ASSERT_EQ(split.section.flags, one.section.flags);
    ASSERT_EQ(split.section.num_records, one.section.num_records);
    ASSERT_EQ(split.section.headers, one.section.headers);
    ASSERT_EQ(split.section.layout, one.section.layout);
    ASSERT_EQ(split.text.text_size, one.text.text_size);
    ASSERT_EQ(split.text.freq, one.text.freq);
}

// returns parameters of the compression with a random number of threads and the default number of start pairs
pfp_rlbwt::comp::compress_params random_params()
{
    return { .num_threads = random_in(1, std::max<uint64_t>(1, omp_get_max_threads())) };
}

/* sets runs_per_block, for half of the calls, to a random value for which the at most n + 1 runs of the BWT of a text
   of length n form at most max_drawn_run_blocks + 1 blocks */
void draw_runs_per_block(pfp_rlbwt::comp::compress_params& params, uint64_t n)
{
    if (random_in(0, 1) == 0) {
        params.runs_per_block = random_in(std::max<uint64_t>(1, (n + 1) / max_drawn_run_blocks), n + 1);
    }
}

// returns parameters of the compression in the FASTA mode with a random number of threads
pfp_rlbwt::comp::compress_params random_fasta_params()
{
    pfp_rlbwt::comp::compress_params params = random_params();
    params.fasta = true;
    return params;
}

// returns a random header: tokens of digits (numbers that often differ from the previous header by a small value),
// letters and other bytes
std::string random_header(std::vector<uint64_t>& numbers)
{
    std::string header;
    const uint64_t num_tokens = random_in(0, max_header_tokens);
    numbers.resize(std::max<uint64_t>(numbers.size(), num_tokens));

    for (uint64_t j = 0; j < num_tokens; j++) {
        switch (random_in(0, 3)) {
            case 0: {
                if (random_in(0, 1) == 0) numbers[j] += random_in(0, max_header_delta);
                else numbers[j] = random_in(0, std::numeric_limits<uint64_t>::max());
                std::string number = std::to_string(numbers[j]);
                if (random_in(0, 3) == 0) number.insert(0, random_in(1, max_leading_zeros), '0');
                header += number;
                break;
            }
            case 1: header += std::string(random_in(1, max_header_word), char('A' + random_in(0, 'z' - 'A'))); break;
            case 2: header += " ._-|/:,"[random_in(0, 7)]; break;
            default: header += char(random_in(0, std::numeric_limits<uint8_t>::max())); break;
        }
    }

    // a header line ends at the first line feed
    std::replace(header.begin(), header.end(), '\n', ' ');
    return header;
}

// a random FASTA file: records with random headers and sequences over ACGT with N runs and lowercase runs, lines of
// a random width, sometimes irregular lines, empty lines, \r\n line ends and a missing final line end
void verify_random_fasta()
{
    const uint64_t max_length = std::max<uint64_t>(1, max_input_size / max_fasta_records);
    const uint64_t num_records = random_in(0, max_fasta_records);
    const bool crlf = random_in(0, 3) == 0;
    const std::string line_end = crlf ? "\r\n" : "\n";
    // a few sequences that the records copy with mutations, so that the text is repetitive
    std::vector<std::string> sources(random_in(1, max_fasta_sources));

    for (std::string& source : sources) {
        source.resize(random_in(0, max_length));
        for (char& c : source) c = "ACGT"[random_in(0, 3)];
    }

    std::string T;
    std::vector<uint64_t> numbers;
    uint64_t width = random_in(1, max_line_width);
    if (random_in(0, 3) == 0) T += "ACGT" + line_end; // a record without a header

    for (uint64_t r = 0; r < num_records; r++) {
        T += ">" + random_header(numbers) + line_end;
        std::string seq = sources[random_in(0, sources.size() - 1)];

        for (uint64_t m = random_in(0, max_mutations); m > 0 && !seq.empty(); m--) {
            const uint64_t i = random_in(0, seq.size() - 1);

            switch (random_in(0, 2)) {
                case 0: seq[i] = "ACGTRYKM"[random_in(0, 7)]; break;
                case 1: seq.replace(i, 0, std::string(random_in(1, max_run), 'N')); break;
                default: {
                    const uint64_t len = std::min<uint64_t>(random_in(1, max_run), seq.size() - i);
                    for (uint64_t k = i; k < i + len; k++) seq[k] = char(std::tolower(seq[k]));
                }
            }
        }

        if (random_in(0, 7) == 0) width = random_in(1, max_line_width);
        const bool irregular = random_in(0, 7) == 0;

        for (uint64_t i = 0; i < seq.size();) {
            const uint64_t len = std::min<uint64_t>(irregular ? random_in(0, max_line_width) : width, seq.size() - i);
            T += seq.substr(i, len) + line_end;
            i += len;
        }

        if (random_in(0, 7) == 0) T += line_end;
    }

    if (random_in(0, 3) == 0 && !T.empty()) T.resize(T.size() - line_end.size());
    verify_fasta_chunks(T);
    pfp_rlbwt::comp::compress_params params = random_fasta_params();
    draw_runs_per_block(params, T.size());
    verify_round_trip(T, params);
}

// random bytes compressed in the FASTA mode
void verify_random_bytes_as_fasta()
{
    std::string T(random_in(0, max_input_size), 0);
    const std::string alphabet = ">\r\nACGTacgtN 1";

    for (char& c : T) {
        c = random_in(0, 1) == 0 ? alphabet[random_in(0, alphabet.size() - 1)]
                                 : char(random_in(0, std::numeric_limits<uint8_t>::max()));
    }

    verify_fasta_chunks(T);
    verify_round_trip(T, random_fasta_params());
}

// a random text over a random alphabet of up to 256 bytes
void verify_random_texts()
{
    const uint64_t n = random_in(1, max_input_size);
    const uint64_t sigma = random_in(1, num_byte_values);
    const uint64_t first = random_in(0, num_byte_values - sigma);
    std::string T(n, 0);
    for (char& c : T) c = char(first + random_in(0, sigma - 1));
    pfp_rlbwt::comp::compress_params params = random_params();
    if (random_in(0, 1) == 0) params.max_start_pairs = random_in(1, max_drawn_start_pairs);
    draw_runs_per_block(params, n);
    verify_round_trip(T, params);
}

// a repetitive text
void verify_repetitive_texts()
{
    const uint8_t max_sym = uint8_t(random_in(1, std::numeric_limits<uint8_t>::max()));
    const std::vector<uint8_t> T = random_repetitive_input<std::vector<uint8_t>>(1, max_input_size, 0, max_sym, gen());
    pfp_rlbwt::comp::compress_params params = random_params();
    draw_runs_per_block(params, T.size());
    verify_round_trip(std::string(T.begin(), T.end()), params);
}

// a text that contains every byte value
void verify_all_byte_values()
{
    std::string T;
    for (uint64_t c = 0; c < num_byte_values; c++) T.push_back(char(c));
    const uint64_t n = random_in(0, max_input_size);
    for (uint64_t i = 0; i < n; i++) T.push_back(T[random_in(0, T.size() - 1)]);
    std::shuffle(T.begin(), T.end(), gen);
    verify_round_trip(T, random_params());
}

// a text of one byte value repeated (the BWT of T$ has two runs)
void verify_one_symbol()
{
    const char c = char(random_in(0, std::numeric_limits<uint8_t>::max()));
    verify_round_trip(std::string(random_in(1, max_input_size), c), random_params());
}

} // namespace

TEST(test_pfp_rlbwt_comp, edge_cases)
{
    for (uint64_t p : decompression_threads) {
        verify_round_trip("", { .num_threads = p });
        verify_round_trip("a", { .num_threads = p });
        verify_round_trip(std::string(1, '\0'), { .num_threads = p });
        verify_round_trip("ab", { .num_threads = p });
        verify_round_trip("aa", { .num_threads = p });
        verify_round_trip("abracadabra", { .num_threads = p, .max_start_pairs = 1 });
        verify_round_trip("abracadabra", { .num_threads = p, .max_start_pairs = 2 });
        verify_round_trip("abracadabra", { .num_threads = p });
        verify_round_trip("abracadabra", { .num_threads = p, .runs_per_block = 1 });
        verify_round_trip("abracadabra", { .num_threads = p, .runs_per_block = 2 });
    }
}

TEST(test_pfp_rlbwt_comp, fasta_edge_cases)
{
    const std::vector<std::string> files = {
        "", // no record
        ">", // an empty header without a line end
        ">\n", // an empty header
        ">\n>\n\n", // empty headers, empty sequences
        ">seq\nACGT\n", // a single record
        ">seq\nACGT", // a missing final line end
        ">seq\r\nACGT\r\nAC\r\n", // \r\n line ends
        ">seq\r\nACGT\r\nAC", // \r\n line ends, missing final line end
        ">seq\r\nAC\nGT\r\n", // a \r that is part of a line
        ">seq\nACGT\r", // a final \r
        "ACGT\nAC\n>seq\nGG\n", // a record without a header
        "\n\n", // empty lines without a header
        ">seq1\nACGTA\nCGTAC\nGT\n>seq2\nACG\nTACG\nT\n>seq3\n\n>seq4\nACGTACGT\n", // varying line widths
        ">a\nACGT\nACGT\n\n>b\nACGT\nAC\n\n\n", // empty lines after records
        ">a\nAC\n\nGT\n", // an empty line inside a record
        ">chr1 soft-masked\nACGTacgtacgtNNNNNNNNacgtACGT\nacgtACGTNNNNnnnnACGTRYKM\nac\n", // lowercase, N, IUPAC
        ">NC_000001.1 x\nACGT\n>NC_000002.1 x\nACGT\n>NC_000010.1 x\nACGT\n>NC_000009.12 y\nA\n", // numbers
        ">r1:54-99\nACGT\n>r1:132-199\nACGT\n>r2:0099-100\nACGT\n>r2:099-00100\nACGT\n", // leading zeros
        ">99999999999999999999\nA\n>100000000000000000000\nA\n>1\nA\n", // numbers longer than 18 digits
        ">\xff\x80 \t\x01\nAC\n>\xff\x81\nAC\n", // other bytes in headers
    };

    for (uint64_t p : fasta_threads) {
        for (const std::string& file : files) verify_round_trip(file, { .num_threads = p, .fasta = true });
    }

    for (const std::string& file : files) {
        for (uint64_t chunk_bytes = 1; chunk_bytes <= file.size(); chunk_bytes++) {
            auto open = [&] {
                return memory_file(reinterpret_cast<const uint8_t*>(file.data()), file.size());
            };
            const fasta_parts one = split_fasta(open, 1, false, file.size() + 1);

            for (uint64_t p : fasta_threads) {
                const fasta_parts split = split_fasta(open, p, false, chunk_bytes);
                ASSERT_EQ(split.section.headers, one.section.headers) << file << " " << chunk_bytes;
                ASSERT_EQ(split.section.layout, one.section.layout) << file << " " << chunk_bytes;
                ASSERT_EQ(split.section.flags, one.section.flags) << file << " " << chunk_bytes;
                ASSERT_EQ(split.section.num_records, one.section.num_records) << file << " " << chunk_bytes;
            }
        }
    }
}

TEST(test_pfp_rlbwt_comp, everything)
{
    // the tests use multiple threads themselves and run one after another
    run_fuzz("pfp-rlbwt-comp", {
        { "random-texts", [](uint64_t) { verify_random_texts(); }, false },
        { "repetitive-texts", [](uint64_t) { verify_repetitive_texts(); }, false },
        { "all-byte-values", [](uint64_t) { verify_all_byte_values(); }, false },
        { "one-symbol", [](uint64_t) { verify_one_symbol(); }, false },
        { "random-fasta", [](uint64_t) { verify_random_fasta(); }, false },
        { "random-bytes-as-fasta", [](uint64_t) { verify_random_bytes_as_fasta(); }, false }
    }, fuzz_iterations(), &max_input_size, min_input_size);
}
