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
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <omp.h>

#include <misc/repetitive_input.hpp>
#include <misc/temp_files.hpp>
#include <pfp_rlbwt/fasta.hpp>
#include <pfp_rlbwt/output/readers.hpp>
#include <pfp_rlbwt/pfp_rlbwt.hpp>
#include "test-fuzz.hpp"
#include "test-references.hpp"

namespace {

// initial maximum length of a fuzzed text; length of the text of the file_source test
constexpr uint64_t max_text_length = 20000;
// maximum length of a fuzzed text, adjusted by run_fuzz
uint64_t text_length = max_text_length;

// maximum alphabet size, window length, modulus and number of threads of a fuzzed construction
constexpr uint64_t max_sigma = 6;
constexpr uint64_t max_window = 8;
constexpr uint64_t max_modulus = 24;
constexpr uint64_t max_threads = 8;

// maximum number of extra sections of the runs file and maximum number of bytes of one
constexpr uint64_t max_extra_sections = 3;
constexpr uint64_t max_extra_section_bytes = 100;

// largest number of runs per block of a fuzzed construction (a multiple of 8) and number of random accesses per reader
constexpr uint64_t max_runs_per_block = 64;
constexpr uint64_t random_accesses = 16;

// the bytes of a random FASTA file, the separators its transformation uses and the largest sample period of its scan
const std::string fasta_bytes = ">>\n\n\n\r\rACGTNacgtx \t";
const std::string fasta_separators = "#\n";
constexpr uint64_t max_sample_period = 64;

// the number of random ranges read from the text of a FASTA file
constexpr uint64_t fasta_random_reads = 8;

thread_local std::mt19937_64 gen(test_seed() + omp_get_thread_num());

// returns a random integer in [lo..hi]
uint64_t random_in(uint64_t lo, uint64_t hi) { return std::uniform_int_distribution<uint64_t>(lo, hi)(gen); }

// the expected output for a text T
struct reference_t {
    // [0..n-1] the BWT, the suffix array and the LCP array of T$ (n = |T$|)
    std::vector<uint64_t> bwt, sa, lcp_full;
    /* [0..r-1] for run i: its start, its symbol, its length, the SA values at its start and at its end and the LCP
       value at its start (0 for run 0) */
    std::vector<uint64_t> start, sym, len, ssa, esa, lcp;
};

// removes files when it goes out of scope
struct file_remover {
    std::vector<std::string> paths;

    ~file_remover()
    {
        std::error_code ec;
        for (const std::string& path : paths) std::filesystem::remove(path, ec);
    }
};

// returns the paths of all output files of pfp_rlbwt::build with a prefix
std::vector<std::string> output_files(const std::string& prefix)
{
    std::vector<std::string> paths;

    for (const char* ext : { pfp_rlbwt::ext_runs, pfp_rlbwt::ext_bwt, pfp_rlbwt::ext_sa, pfp_rlbwt::ext_lcp }) {
        paths.push_back(prefix + ext);
    }

    return paths;
}

// computes the expected output for a text T from the suffix array of T$
template <typename sym_t>
reference_t compute_reference(const std::vector<sym_t>& text)
{
    const uint64_t n = text.size() + 1;
    // T as integers; the terminator $ is smaller than every symbol (pfp_rlbwt::min_symbol > pfp_rlbwt::terminator)
    const std::vector<uint64_t> values(text.begin(), text.end());
    const std::vector<int64_t> SA = compute_reference_sa<uint64_t>(values, n - 1, omp_get_max_threads());
    // returns T$[i]
    auto at = [&](uint64_t i) -> uint64_t { return i + 1 == n ? pfp_rlbwt::terminator : values[i]; };

    // [0..n-1] the LCP array of T$
    const std::vector<uint64_t> LCP = compute_reference_lcp(values, SA);

    reference_t expected;

    for (uint64_t i = 0; i < n; i++) {
        const uint64_t c = SA[i] == 0 ? pfp_rlbwt::terminator : at(SA[i] - 1);
        expected.bwt.push_back(c);
        expected.sa.push_back(SA[i]);
        expected.lcp_full.push_back(LCP[i]);

        if (i == 0 || c != expected.sym.back()) {
            if (i != 0) expected.esa.push_back(SA[i - 1]);
            expected.start.push_back(i);
            expected.sym.push_back(c);
            expected.len.push_back(0);
            expected.ssa.push_back(SA[i]);
            expected.lcp.push_back(i == 0 ? 0 : LCP[i]);
        }

        expected.len.back()++;
    }

    expected.esa.push_back(SA[n - 1]);
    return expected;
}

// reads the bytes of a file
std::vector<uint8_t> read_bytes(const std::string& file)
{
    std::ifstream in(file, std::ios::binary);
    return std::vector<uint8_t>(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

// returns the number of reads of the file of a reader on the disk (0 for a representation in RAM)
template <typename reader_t>
uint64_t reads_of(const reader_t& reader)
{
    if constexpr (requires { reader.num_reads(); }) return reader.num_reads();
    else return 0;
}

// checks a run read from a runs file against run i of the expected output
void check_run(const pfp_rlbwt::run_record& rec, const reference_t& expected, uint64_t i, bool sa, bool lcp)
{
    ASSERT_EQ(rec.start, expected.start[i]) << "run " << i;
    ASSERT_EQ(rec.length, expected.len[i]) << "run " << i;
    ASSERT_EQ(rec.symbol, expected.sym[i]) << "run " << i;
    ASSERT_EQ(rec.sa_start, sa ? expected.ssa[i] : 0) << "run " << i;
    ASSERT_EQ(rec.sa_end, sa ? expected.esa[i] : 0) << "run " << i;
    ASSERT_EQ(rec.lcp, lcp ? expected.lcp[i] : 0) << "run " << i;
}

/* checks the random access of a reader on the disk or in RAM against the expected output: record for every run (in RAM)
   or for random runs (on the disk), the other accesses for random runs and BWT positions; on the disk, record reads the
   file record_reads times and each other access once */
template <typename reader_t>
void check_random_access(reader_t& reader, const reference_t& expected, uint64_t record_reads, bool sa, bool lcp)
{
    const uint64_t r = expected.sym.size(), n = expected.bwt.size();
    ASSERT_EQ(reader.num_runs(), r);
    ASSERT_EQ(reader.text_length(), n);

    constexpr bool on_disk = requires { reader.num_reads(); };

    // returns the result of access, after checking that it has read the file reads times (on the disk)
    auto counted = [&](uint64_t reads, auto access) {
        const uint64_t before = reads_of(reader);
        const auto value = access();
        if constexpr (on_disk) EXPECT_EQ(reads_of(reader) - before, reads);
        return value;
    };

    for (uint64_t k = 0; k < (on_disk ? random_accesses : r); k++) {
        const uint64_t i = on_disk ? random_in(0, r - 1) : k;
        ASSERT_NO_FATAL_FAILURE(check_run(counted(record_reads, [&] { return reader.record(i); }), expected, i, sa,
            lcp));
    }

    for (uint64_t k = 0; k < random_accesses; k++) {
        const uint64_t i = random_in(0, r - 1);
        const pfp_rlbwt::run_head head = counted(1, [&] { return reader.run(i); });
        ASSERT_EQ(head.start, expected.start[i]) << "run " << i;
        ASSERT_EQ(head.symbol, expected.sym[i]) << "run " << i;
        ASSERT_EQ(counted(1, [&] { return reader.length(i); }), expected.len[i]) << "run " << i;

        if (sa) {
            ASSERT_EQ(counted(1, [&] { return reader.sa_start(i); }), expected.ssa[i]) << "run " << i;
            ASSERT_EQ(counted(1, [&] { return reader.sa_end(i); }), expected.esa[i]) << "run " << i;
        }

        if (lcp) ASSERT_EQ(counted(1, [&] { return reader.lcp(i); }), expected.lcp[i]) << "run " << i;
        const uint64_t j = random_in(0, n - 1);
        const uint64_t run = std::upper_bound(expected.start.begin(), expected.start.end(), j) - expected.start.begin() - 1;
        ASSERT_EQ(counted(1, [&] { return reader.run_of(j); }), run) << "position " << j;
    }
}

// builds the output for a text with the given parameters and compares it with the expected output
template <typename sym_t>
void verify_build(const std::vector<sym_t>& text, const pfp_rlbwt::params& params)
{
    const std::string prefix = temporary_file_path("pfp-rlbwt-test-");
    const file_remover remover { output_files(prefix) };
    pfp_rlbwt::memory_source<sym_t> source(text.data(), text.size());
    const pfp_rlbwt::result result = pfp_rlbwt::build(source, prefix, params);
    const reference_t expected = compute_reference(text);
    const uint64_t n = text.size() + 1;
    const uint64_t r = expected.sym.size();
    const bool sa = params.write_sa_samples, lcp = params.write_lcp_samples;
    const bool plain = params.format == pfp_rlbwt::file_format::plain;
    const uint8_t position_bytes = byte_width(n);
    EXPECT_EQ(result.n_with_terminator, n);
    ASSERT_EQ(result.r, r);

    // the runs file stores exactly the sections of its format, and the temporary files are removed
    const std::string runs = prefix + pfp_rlbwt::ext_runs;
    pfp_rlbwt::run_data_file file;
    ASSERT_TRUE(file.open(runs));
    ASSERT_EQ(file.format(), params.format);
    std::vector<uint64_t> ids = plain
        ? std::vector<uint64_t> { pfp_rlbwt::section_fields, pfp_rlbwt::section_records, pfp_rlbwt::section_record_samples }
        : std::vector<uint64_t> { pfp_rlbwt::section_fields, pfp_rlbwt::section_run_blocks,
            pfp_rlbwt::section_run_block_index, pfp_rlbwt::section_symbol_code };

    if (!plain && sa) ids.insert(ids.end(), { pfp_rlbwt::section_sa_run_starts, pfp_rlbwt::section_sa_run_ends });
    if (!plain && lcp) ids.insert(ids.end(), { pfp_rlbwt::section_lcp_blocks, pfp_rlbwt::section_lcp_block_index });
    for (uint64_t id : ids) EXPECT_TRUE(file.contains(id)) << id;
    EXPECT_EQ(file.sections().size(), ids.size() + params.extra_sections.size());
    EXPECT_FALSE(std::filesystem::exists(runs + pfp_rlbwt::ext_temporary_records));
    EXPECT_FALSE(std::filesystem::exists(runs + "." + std::to_string(pfp_rlbwt::section_lcp_blocks)));
    for (const char* ext : { pfp_rlbwt::ext_dictionary_sa, pfp_rlbwt::ext_phrase_run_ends,
             pfp_rlbwt::ext_phrase_run_symbols, pfp_rlbwt::ext_parse_isa }) {
        EXPECT_FALSE(std::filesystem::exists(prefix + ext)) << ext;
    }

    // the sections are stored one after another, followed by the table of contents
    uint64_t end = 0;

    for (const pfp_rlbwt::run_data_entry& e : file.sections()) {
        EXPECT_EQ(e.offset, end) << e.id;
        end = e.offset + e.size;
    }

    EXPECT_EQ(end + file.sections().size() * sizeof(pfp_rlbwt::run_data_entry) + sizeof(uint64_t),
        std::filesystem::file_size(runs));
    const std::vector<uint8_t> bytes = read_bytes(runs);

    // the extra sections
    for (const pfp_rlbwt::extra_section& extra : params.extra_sections) {
        std::ostringstream extra_bytes;
        extra.write(extra_bytes);
        const pfp_rlbwt::run_data_entry* e = file.find(extra.id);
        ASSERT_NE(e, nullptr) << extra.id;
        EXPECT_EQ(std::string(bytes.begin() + e->offset, bytes.begin() + e->offset + e->size), extra_bytes.str())
            << extra.id;
    }

    // the fields
    const pfp_rlbwt::run_data_fields& f = file.fields();
    EXPECT_EQ(f.n, n);
    EXPECT_EQ(f.r, r);
    EXPECT_EQ(f.max_symbol, *std::max_element(expected.sym.begin(), expected.sym.end()));
    EXPECT_EQ(f.max_run_length, *std::max_element(expected.len.begin(), expected.len.end()));
    EXPECT_EQ(f.sa_samples, sa);
    EXPECT_EQ(f.lcp_samples, lcp);
    EXPECT_EQ(f.position_bytes, position_bytes);
    EXPECT_EQ(f.symbol_bytes, sizeof(sym_t));

    // the plain format: the records with byte widths, read directly from the file
    if (plain) {
        const pfp_rlbwt::record_layout layout(f);
        const uint64_t records = file.section(pfp_rlbwt::section_records).offset;
        EXPECT_EQ(layout.bytes, position_bytes * (1 + 2 * sa + lcp) + sizeof(sym_t));

        for (uint64_t i = 0; i < r; i++) {
            const uint8_t* rec = bytes.data() + records + i * layout.bytes;
            ASSERT_EQ(pfp_rlbwt::read_bytes(rec, position_bytes), expected.start[i]) << "run " << i;
            ASSERT_EQ(pfp_rlbwt::read_bytes(rec + position_bytes, sizeof(sym_t)), expected.sym[i]) << "run " << i;
            const uint8_t* values = rec + position_bytes + sizeof(sym_t);

            if (sa) {
                ASSERT_EQ(pfp_rlbwt::read_bytes(values, position_bytes), expected.ssa[i]) << "run " << i;
                ASSERT_EQ(pfp_rlbwt::read_bytes(values + position_bytes, position_bytes), expected.esa[i]) << "run " << i;
            }

            if (lcp) ASSERT_EQ(pfp_rlbwt::read_bytes(values + 2 * sa * position_bytes, position_bytes), expected.lcp[i]);
        }
    }

    // the streams: all runs, and a random range
    for (bool whole : { true, false }) {
        const uint64_t first = whole ? 0 : random_in(0, r - 1);
        const uint64_t last = whole ? r : random_in(first, r);
        pfp_rlbwt::run_stream stream(file, first, last);
        pfp_rlbwt::run_record rec;

        for (uint64_t i = first; i < last; i++) {
            ASSERT_FALSE(stream.done());
            stream.next(rec);
            ASSERT_NO_FATAL_FAILURE(check_run(rec, expected, i, sa, lcp));
        }

        ASSERT_TRUE(stream.done());
    }

    // the readers on the disk and in RAM
    if (plain) {
        pfp_rlbwt::plain_disk_reader reader(file);
        ASSERT_NO_FATAL_FAILURE(check_random_access(reader, expected, 1, sa, lcp));
    } else {
        pfp_rlbwt::compressed_disk_reader reader(file);
        ASSERT_NO_FATAL_FAILURE(check_random_access(reader, expected, 1 + 2 * sa + lcp, sa, lcp));
        pfp_rlbwt::compressed_run_data in_ram(file, random_in(1, max_threads));
        ASSERT_NO_FATAL_FAILURE(check_random_access(in_ram, expected, 0, sa, lcp));
    }

    pfp_rlbwt::packed_run_data packed(file, random_in(1, max_threads));
    ASSERT_NO_FATAL_FAILURE(check_random_access(packed, expected, 0, sa, lcp));

    // the full BWT
    if (params.write_bwt) {
        const std::vector<uint8_t> full = read_bytes(prefix + pfp_rlbwt::ext_bwt);
        ASSERT_EQ(full.size(), n * sizeof(sym_t));

        for (uint64_t i = 0; i < n; i++) {
            sym_t c;
            std::memcpy(&c, full.data() + i * sizeof(sym_t), sizeof(sym_t));
            ASSERT_EQ(uint64_t(c), expected.bwt[i]) << "BWT[" << i << "]";
        }
    }

    // checks the full SA or LCP array in the file with extension ext (position_bytes bytes per value)
    auto check_full = [&](const char* ext, const std::vector<uint64_t>& expected_values) {
        const std::vector<uint8_t> full = read_bytes(prefix + ext);
        ASSERT_EQ(full.size(), n * position_bytes) << ext;

        for (uint64_t i = 0; i < n; i++) {
            ASSERT_EQ(pfp_rlbwt::read_bytes(full.data() + i * position_bytes, position_bytes), expected_values[i])
                << ext << " " << i;
        }
    };

    if (params.write_sa) ASSERT_NO_FATAL_FAILURE(check_full(pfp_rlbwt::ext_sa, expected.sa));
    if (params.write_lcp) ASSERT_NO_FATAL_FAILURE(check_full(pfp_rlbwt::ext_lcp, expected.lcp_full));
}

// draws the parameters of a construction
pfp_rlbwt::params random_params()
{
    pfp_rlbwt::params params;
    params.window = random_in(1, max_window);
    params.modulus = random_in(1, max_modulus);
    params.num_threads = random_in(1, max_threads);
    // the values at the run boundaries and the full arrays are written with probability 1/2 each
    params.write_sa_samples = random_in(0, 1) != 0;
    params.write_lcp_samples = random_in(0, 1) != 0;
    params.write_bwt = random_in(0, 1) != 0;
    params.write_sa = random_in(0, 1) != 0;
    params.write_lcp = random_in(0, 1) != 0;
    params.format = random_in(0, 1) == 0 ? pfp_rlbwt::file_format::plain : pfp_rlbwt::file_format::compressed;
    // the block sizes: multiples of 8 runs up to max_runs_per_block runs, or the defaults
    if (random_in(0, 1) == 0) {
        params.runs_per_block = pfp_rlbwt::byte_bits * random_in(1, max_runs_per_block / pfp_rlbwt::byte_bits);
        params.runs_per_symbol_sample = random_in(1, params.runs_per_block);
    }
    // 0 (the fewest bytes) or any width in [min_value_bytes..max_value_bytes]
    params.position_bytes = random_in(0, 1) == 0 ? 0 : random_in(pfp_rlbwt::min_value_bytes, pfp_rlbwt::max_value_bytes);

    // up to max_extra_sections extra sections with distinct identifiers and random bytes
    for (uint64_t i = random_in(0, max_extra_sections); i > 0; i--) {
        const uint64_t id = pfp_rlbwt::first_extra_section + params.extra_sections.size() + random_in(0, 1);
        if (!params.extra_sections.empty() && params.extra_sections.back().id >= id) continue;
        std::string bytes(random_in(0, max_extra_section_bytes), 0);
        for (char& c : bytes) c = char(random_in(0, UINT8_MAX));
        params.extra_sections.push_back({ id, [bytes](std::ostream& out) { out << bytes; } });
    }

    return params;
}

/* draws a text with symbols of type sym_t: random or repetitive over sigma symbols, which are min_symbol, ... or, for
   wider symbols, random values */
template <typename sym_t>
std::vector<sym_t> random_text(bool repetitive)
{
    const uint64_t sigma = random_in(1, max_sigma);
    std::vector<sym_t> alphabet;

    for (uint64_t c = 0; c < sigma; c++) {
        alphabet.push_back(sizeof(sym_t) == 1 || random_in(0, 1) == 0 ? sym_t(pfp_rlbwt::min_symbol + c)
            : sym_t(random_in(pfp_rlbwt::min_symbol, max_value_of<sym_t>())));
    }

    std::vector<sym_t> text;

    if (repetitive) {
        const std::string s = generate_repetitive_input<std::string>(
            random_repetitive_params(random_in(1, text_length), gen), char(0), char(sigma - 1));
        for (char c : s) text.push_back(alphabet[uint8_t(c)]);
    } else {
        text.resize(random_in(0, text_length));
        for (sym_t& c : text) c = alphabet[random_in(0, sigma - 1)];
    }

    return text;
}

// returns a functionality of run_fuzz that checks the construction on random texts with symbols of type sym_t
template <typename sym_t>
fuzz_functionality fuzz(const std::string& name, bool repetitive)
{
    return { name, [repetitive](uint64_t) { verify_build(random_text<sym_t>(repetitive), random_params()); }, false };
}

// returns a transformation of FASTA files with random options
pfp_rlbwt::fasta_transform random_fasta_transform()
{
    pfp_rlbwt::fasta_transform transform(random_in(0, 1) == 0 ? "ACGT" : "ACGTN",
        uint8_t(fasta_separators[random_in(0, fasta_separators.size() - 1)]), random_in(0, 1) == 0);
    if (random_in(0, 1) == 0) transform.in_alphabet.fill(1);
    transform.separator_before_first_record = random_in(0, 1) == 0;
    transform.lines_before_first_header = random_in(0, 1) == 0;
    transform.carriage_return = pfp_rlbwt::fasta_transform::carriage_return_t(
        random_in(pfp_rlbwt::fasta_transform::_dropped, pfp_rlbwt::fasta_transform::_kept));
    return transform;
}

// a line event of scan_fasta_lines: '>' header_start, 'b' header bytes, 's' sequence symbols, 'e' line_end
struct fasta_event {
    char kind;
    std::string bytes; // the header bytes ('b')
    uint64_t value = 0; // the number of symbols ('s'), 2 * header + after_cr ('e')
    bool operator==(const fasta_event&) const = default;
};

// appends an event; consecutive header bytes and consecutive sequence symbols are merged
void add_fasta_event(std::vector<fasta_event>& events, fasta_event event)
{
    if ((event.kind == 'b' && event.bytes.empty()) || (event.kind == 's' && event.value == 0)) return;

    if (!events.empty() && events.back().kind == event.kind && (event.kind == 'b' || event.kind == 's')) {
        events.back().bytes += event.bytes;
        events.back().value += event.value;
    } else {
        events.push_back(event);
    }
}

// a chunk of scan_fasta_lines that records its line events
struct fasta_event_chunk : pfp_rlbwt::fasta_chunk {
    std::vector<fasta_event> events;

    void header_start() { add_fasta_event(events, { '>' }); }
    void header_bytes(const uint8_t* b, uint64_t len) { add_fasta_event(events, { 'b', { (const char*) b, len } }); }
    void sequence_symbols(uint64_t count) { add_fasta_event(events, { 's', {}, count }); }
    void line_end(bool header, bool after_cr) { add_fasta_event(events, { 'e', {}, uint64_t(2 * header + after_cr) }); }
};

/* checks the parallel scan of a random FASTA file (with random chunk sizes and sample periods), its line events and
   the text that fasta_source reads from it (from the file and from memory, reversed or not, mapped or not, with the
   samples of the scan or with samples it computes itself) against the transformation applied to the bytes one after
   another */
void verify_fasta()
{
    const pfp_rlbwt::fasta_transform transform = random_fasta_transform();
    std::string file(random_in(0, text_length), 0);
    for (char& c : file) c = fasta_bytes[random_in(0, fasta_bytes.size() - 1)];
    if (random_in(0, 1) == 0 && !file.empty()) file[0] = char(pfp_rlbwt::fasta_header_marker);

    std::string text;
    std::array<uint64_t, 256> freq {};
    std::vector<fasta_event> events;
    pfp_rlbwt::fasta_transform::state_t state;

    for (uint64_t i = 0; i < file.size(); i++) {
        uint8_t c = 0;
        const uint8_t byte = uint8_t(file[i]);
        const bool line_feed_follows = i + 1 < file.size() && uint8_t(file[i + 1]) == pfp_rlbwt::fasta_line_feed;
        const bool dropped_cr = byte == pfp_rlbwt::fasta_carriage_return
            && (transform.carriage_return == pfp_rlbwt::fasta_transform::_dropped
                || (transform.carriage_return == pfp_rlbwt::fasta_transform::_dropped_before_line_feed
                    && line_feed_follows));
        const pfp_rlbwt::fasta_transform::state_t before = state;
        const pfp_rlbwt::fasta_transform::step_t step = transform.step(byte, state, c, [&] { return line_feed_follows; });

        if (byte == pfp_rlbwt::fasta_line_feed) {
            add_fasta_event(events, { 'e', {}, uint64_t(2 * before.in_header
                + (i > 0 && uint8_t(file[i - 1]) == pfp_rlbwt::fasta_carriage_return)) });
        } else if (state.in_header && !before.in_header) {
            add_fasta_event(events, { '>' });
        } else if (before.in_header && !dropped_cr) {
            add_fasta_event(events, { 'b', std::string(1, char(byte)) });
        }

        if (step == pfp_rlbwt::fasta_transform::_symbol) add_fasta_event(events, { 's', {}, 1 });
        if (step == pfp_rlbwt::fasta_transform::_separator) c = transform.separator;
        if (step == pfp_rlbwt::fasta_transform::_nothing) continue;
        text.push_back(char(c));
        freq[c]++;
    }

    const pfp_rlbwt::memory_file bytes(reinterpret_cast<const uint8_t*>(file.data()), file.size());
    const uint64_t num_threads = random_in(1, max_threads);

    std::vector<fasta_event> scanned_events;
    const pfp_rlbwt::fasta_scan scan = pfp_rlbwt::scan_fasta_lines<fasta_event_chunk>([&] { return bytes; }, transform,
        { .num_threads = num_threads, .sample_period = random_in(1, max_sample_period),
          .min_chunk_bytes = random_in(1, file.size() + 1) },
        [&](fasta_event_chunk& chunk, uint64_t) {
            for (const fasta_event& event : chunk.events) add_fasta_event(scanned_events, event);
        });
    ASSERT_EQ(scan.text_size, text.size());
    ASSERT_EQ(scan.freq, freq);
    ASSERT_TRUE(scanned_events == events);

    const std::string path = temporary_file_path("pfp-rlbwt-test-fasta-");
    const file_remover remover { { path } };
    std::ofstream(path, std::ios::binary).write(file.data(), file.size());
    // reverses the byte values
    std::array<uint8_t, 256> map;
    for (uint64_t c = 0; c < map.size(); c++) map[c] = uint8_t(map.size() - 1 - c);

    for (bool in_memory : { false, true }) {
        for (bool reverse : { false, true }) {
            const bool mapped = random_in(0, 1) == 0;
            const pfp_rlbwt::fasta_sample* samples = random_in(0, 1) == 0 ? scan.samples.data() : nullptr;
            std::string expected = text;
            if (mapped) for (char& c : expected) c = char(map[uint8_t(c)]);
            if (reverse) std::reverse(expected.begin(), expected.end());

            const pfp_rlbwt::fasta_source<uint8_t> src = in_memory
                ? pfp_rlbwt::fasta_source<uint8_t>(reinterpret_cast<const uint8_t*>(file.data()), file.size(),
                    text.size(), transform, samples, scan.samples.size(), reverse, mapped ? &map : nullptr,
                    num_threads)
                : pfp_rlbwt::fasta_source<uint8_t>(path, text.size(), transform, samples, scan.samples.size(),
                    reverse, mapped ? &map : nullptr, num_threads);
            if (text.empty()) continue;
            std::unique_ptr<pfp_rlbwt::text_source<uint8_t>::reader> reader = src.new_reader();
            std::string read(text.size(), 0);

            // consecutive ranges, then random ones
            for (uint64_t pos = 0; pos < text.size();) {
                const uint64_t len = random_in(1, text.size() - pos);
                reader->read(pos, len, reinterpret_cast<uint8_t*>(read.data() + pos));
                pos += len;
            }

            ASSERT_EQ(read, expected) << "in memory: " << in_memory << ", reverse: " << reverse;

            for (uint64_t k = 0; k < fasta_random_reads; k++) {
                const uint64_t pos = random_in(0, text.size() - 1);
                const uint64_t len = random_in(1, text.size() - pos);
                reader->read(pos, len, reinterpret_cast<uint8_t*>(read.data()));
                ASSERT_EQ(read.substr(0, len), expected.substr(pos, len)) << pos << " " << len;
            }
        }
    }

    // 16-bit symbols: the byte c is read as c + 256
    std::array<uint16_t, 256> wide_map;
    for (uint64_t c = 0; c < wide_map.size(); c++) wide_map[c] = uint16_t(c + wide_map.size());
    const pfp_rlbwt::fasta_source<uint16_t> wide(reinterpret_cast<const uint8_t*>(file.data()), file.size(),
        text.size(), transform, scan.samples.data(), scan.samples.size(), false, &wide_map);
    if (text.empty()) return;
    std::vector<uint16_t> wide_read(text.size());
    wide.new_reader()->read(0, text.size(), wide_read.data());
    for (uint64_t i = 0; i < text.size(); i++) ASSERT_EQ(wide_read[i], wide_map[uint8_t(text[i])]) << i;
}

} // namespace

TEST(test_pfp_rlbwt, fuzz)
{
    run_fuzz("pfp-rlbwt", {
        fuzz<uint8_t>("random 8 bit", false),
        fuzz<uint8_t>("repetitive 8 bit", true),
        fuzz<uint16_t>("repetitive 16 bit", true),
        fuzz<uint24_t>("repetitive 24 bit", true),
        fuzz<uint32_t>("repetitive 32 bit", true),
        fuzz<uint40_t>("random 40 bit", false),
        fuzz<uint48_t>("random 48 bit", false),
        fuzz<uint48_t>("repetitive 48 bit", true),
        fuzz<uint56_t>("repetitive 56 bit", true),
        fuzz<uint64_t>("random 64 bit", false),
        fuzz<uint64_t>("repetitive 64 bit", true),
    }, fuzz_iterations(), &text_length);
}

TEST(test_pfp_rlbwt, fasta_source)
{
    run_fuzz("pfp-rlbwt-fasta", { { "random FASTA files", [](uint64_t) { verify_fasta(); }, false } },
        fuzz_iterations(), &text_length);
}

TEST(test_pfp_rlbwt, file_source)
{
    const std::string text = generate_repetitive_input<std::string>(
        random_repetitive_params(max_text_length, gen), 'a', 'd');
    const std::string file = temporary_file_path("pfp-rlbwt-test-file-");
    const file_remover remover { { file } };
    std::ofstream(file, std::ios::binary).write(text.data(), text.size());

    for (bool reverse : { false, true }) {
        // swaps 'a' and 'd'
        std::array<uint8_t, 256> map;
        for (uint64_t c = 0; c < map.size(); c++) map[c] = uint8_t(c == 'a' ? 'd' : c == 'd' ? 'a' : c);

        for (bool mapped : { false, true }) {
            std::string expected = text;
            if (mapped) for (char& c : expected) c = char(map[uint8_t(c)]);
            if (reverse) std::reverse(expected.begin(), expected.end());
            pfp_rlbwt::file_source<uint8_t> src(file, reverse, mapped ? &map : nullptr);
            std::unique_ptr<pfp_rlbwt::text_source<uint8_t>::reader> reader = src.new_reader();
            std::string read(text.size(), 0);
            reader->read(0, text.size(), (uint8_t*) read.data());
            EXPECT_EQ(read, expected);
        }
    }

    // a file of 16-bit symbols, read reversed
    std::vector<uint16_t> wide(text.begin(), text.end());
    std::ofstream(file, std::ios::binary).write((const char*) wide.data(), wide.size() * sizeof(uint16_t));
    pfp_rlbwt::file_source<uint16_t> src(file, true);
    std::vector<uint16_t> read(wide.size());
    src.new_reader()->read(0, wide.size(), read.data());
    std::reverse(wide.begin(), wide.end());
    EXPECT_EQ(read, wide);

    // a file of 24-bit symbols, read reversed
    std::vector<uint24_t> packed;
    for (uint64_t i = 0; i < text.size(); i++) packed.emplace_back(uint64_t(uint8_t(text[i])) << (8 * (i % 3)));
    std::ofstream(file, std::ios::binary).write((const char*) packed.data(), packed.size() * sizeof(uint24_t));
    pfp_rlbwt::file_source<uint24_t> packed_src(file, true);
    std::vector<uint24_t> packed_read(packed.size());
    packed_src.new_reader()->read(0, packed.size(), packed_read.data());
    std::reverse(packed.begin(), packed.end());
    for (uint64_t i = 0; i < packed.size(); i++) ASSERT_EQ(uint64_t(packed_read[i]), uint64_t(packed[i])) << i;
}

// checks that a window that overflows the text length, invalid block sizes and invalid extra sections throw, and that
// the default block sizes work with all samples in both formats
TEST(test_pfp_rlbwt, parameter_limits)
{
    const std::vector<uint8_t> text = random_text<uint8_t>(true);
    pfp_rlbwt::memory_source<uint8_t> source(text.data(), text.size());
    const std::string prefix = temporary_file_path("pfp-rlbwt-test-");
    const file_remover remover { output_files(prefix) };

    // n + window does not fit into 64 bits
    pfp_rlbwt::params huge_window;
    huge_window.window = std::numeric_limits<uint64_t>::max();
    EXPECT_THROW(pfp_rlbwt::build(source, prefix, huge_window), std::invalid_argument);

    // blocks that are empty or no multiple of 8 runs, symbol samples beyond a block
    for (uint64_t runs_per_block : { uint64_t { 0 }, uint64_t { pfp_rlbwt::byte_bits + 1 } }) {
        pfp_rlbwt::params invalid;
        invalid.runs_per_block = runs_per_block;
        EXPECT_THROW(pfp_rlbwt::build(source, prefix, invalid), std::invalid_argument);
    }

    for (uint64_t period : { uint64_t { 0 }, pfp_rlbwt::params().runs_per_block + 1 }) {
        pfp_rlbwt::params invalid;
        invalid.runs_per_symbol_sample = period;
        EXPECT_THROW(pfp_rlbwt::build(source, prefix, invalid), std::invalid_argument);
    }

    // an extra section with an identifier of a section of build, and two extra sections with the same identifier
    const auto write_nothing = [](std::ostream&) { };
    pfp_rlbwt::params reserved_id;
    reserved_id.extra_sections.push_back({ pfp_rlbwt::section_lcp_blocks, write_nothing });
    EXPECT_THROW(pfp_rlbwt::build(source, prefix, reserved_id), std::invalid_argument);
    pfp_rlbwt::params equal_ids;
    equal_ids.extra_sections.push_back({ pfp_rlbwt::first_extra_section, write_nothing });
    equal_ids.extra_sections.push_back({ pfp_rlbwt::first_extra_section, write_nothing });
    EXPECT_THROW(pfp_rlbwt::build(source, prefix, equal_ids), std::invalid_argument);

    for (pfp_rlbwt::file_format format : { pfp_rlbwt::file_format::plain, pfp_rlbwt::file_format::compressed }) {
        pfp_rlbwt::params all;
        all.format = format;
        all.write_sa_samples = all.write_lcp_samples = all.write_sa = all.write_lcp = true;
        verify_build(text, all);
    }
}
