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

#ifndef PFP_RLBWT_PFP_RLBWT_TYPES_HPP
#define PFP_RLBWT_PFP_RLBWT_TYPES_HPP

#include <algorithm>
#include <bit>
#include <cstdint>
#include <cstring>
#include <functional>
#include <limits>
#include <memory>
#include <ostream>
#include <vector>

namespace pfp_rlbwt {

// the smallest symbol a text may contain (0 and 1 are used internally)
inline constexpr uint64_t min_symbol = 2;

// the symbol of the end of the text T$ in the BWT
inline constexpr uint64_t terminator = 0;

/* largest width (in bits) of a bit-packed value of the output: a value is read with one unaligned 64-bit load, which
   also holds the bits of the byte before it (see read_bits) */
inline constexpr uint8_t max_packed_bits =
    std::numeric_limits<uint64_t>::digits - (std::numeric_limits<uint8_t>::digits - 1);

/**
 * @brief an unsigned integer stored in its `bytes` least significant bytes (little endian); the type of the
 *        symbols of a text with 24-, 40- or 48-bit symbols
 */
template <uint8_t bytes>
struct packed_symbol {
    static_assert(std::endian::native == std::endian::little,
        "a packed_symbol stores the least significant bytes of a value, which assumes a little-endian layout");

    uint8_t data[bytes];

    packed_symbol() = default;
    packed_symbol(uint64_t value) { std::memcpy(data, &value, bytes); }

    operator uint64_t() const
    {
        uint64_t value = 0;
        std::memcpy(&value, data, bytes);
        return value;
    }
};

using sym24_t = packed_symbol<3>;
using sym40_t = packed_symbol<5>;
using sym48_t = packed_symbol<6>;

/**
 * @brief a text whose symbols are read in blocks by several threads at once (one reader per thread)
 * @tparam sym_t type of a symbol (uint8_t, uint16_t, sym24_t, uint32_t, sym40_t or sym48_t)
 */
template <typename sym_t = uint8_t>
class text_source {
public:
    // reads symbols of the text; positions are usually read in ascending order
    class reader {
    public:
        virtual ~reader() = default;

        /**
         * @brief copies the symbols T[pos..pos+len-1] to dst
         * @param pos [0..n-1] first position
         * @param len [1..n-pos] number of symbols
         * @param dst destination of len symbols
         */
        virtual void read(uint64_t pos, uint64_t len, sym_t* dst) = 0;
    };

    virtual ~text_source() = default;

    // returns the length n of the text (without the terminator)
    virtual uint64_t size() const = 0;

    // returns a new reader of the text
    virtual std::unique_ptr<reader> new_reader() const = 0;
};

// returns the fewest whole bytes (at least one) that hold the value x
inline constexpr uint8_t byte_width(uint64_t x) { return uint8_t(std::max<uint64_t>(1, (std::bit_width(x) + 7) / 8)); }

// the formats of the runs file (see ext_runs)
enum class file_format : uint8_t {
    plain, // one record per run, each field with a fixed number of bytes (see section_records)
    compressed // the runs, the SA samples and the LCP samples compressed, each readable with O(1) disk reads
};

/* the sections of the runs file (see ext_runs) that build writes; a caller adds further sections with identifiers
   from first_extra_section on (see params::extra_sections). Which sections a file stores tells its format. Values are
   stored little endian; bit-packed data stores the bits of a byte least significant bit first. */
enum run_data_section : uint64_t {
    section_fields, // the numbers that describe the file (see run_data_fields)
    /* plain: [0..r-1] the record of run i at offset i * record_bytes: the BWT position of its start and its symbol and,
       if sa_samples is set, the SA values at its start and at its end and, if lcp_samples is set, the LCP value at its
       start (see record_layout) */
    section_records,
    /* plain: [0..ceil(r / runs_per_block)-1] the start of run b * runs_per_block (position_bytes bytes each), so that a
       reader finds the run containing a BWT position with one disk read */
    section_record_samples,
    // compressed: the runs in blocks of runs_per_block runs (see run_block)
    section_run_blocks,
    /* compressed: [0..ceil(r / runs_per_block)-1] the offset in bytes of block b in section_run_blocks
       (run_block_offset_bits bits) and the start of its first run (position_bits bits) */
    section_run_block_index,
    section_symbol_code, // compressed: the canonical Huffman code of the run symbols (see symbol_code)
    section_sa_run_starts, // compressed: the SA value at the start of each run (position_bits bits each)
    section_sa_run_ends, // compressed: the SA value at the end of each run (position_bits bits each)
    section_lcp_blocks, // compressed: the LCP values at the run starts in blocks of runs_per_block values (see lcp_block)
    // compressed: [0..ceil(r / runs_per_block)-1] the offset in bytes of block b in section_lcp_blocks
    // (lcp_block_offset_bits bits each)
    section_lcp_block_index,
    first_extra_section // the smallest identifier of a section added by a caller
};

// an entry of the table of contents of the runs file (see ext_runs)
struct run_data_entry {
    uint64_t id; // identifier of the section (see run_data_section)
    uint64_t offset; // offset (in bytes) of the section in the file
    uint64_t size; // size (in bytes) of the section
};

static_assert(sizeof(run_data_entry) == 3 * sizeof(uint64_t));

// the numbers that describe a runs file (the section section_fields, one 64-bit integer each)
struct run_data_fields {
    uint64_t n = 0; // length of T$
    uint64_t r = 0; // number of runs
    uint64_t max_symbol = 0; // the largest run symbol
    uint64_t max_run_length = 0; // the largest run length
    uint64_t sa_samples = 0; // 1 <=> the SA values at the starts and at the ends of the runs are stored
    uint64_t lcp_samples = 0; // 1 <=> the LCP values at the starts of the runs are stored
    // plain: number of runs per record sample; compressed: number of runs per block and LCP values per LCP block
    uint64_t runs_per_block = 0;
    uint64_t runs_per_symbol_sample = 0; // compressed: number of runs per sample of the symbol codes (see run_block)
    uint64_t position_bytes = 0; // plain: bytes per run start, SA value and LCP value: byte_width(n)
    uint64_t symbol_bytes = 0; // plain: bytes per run symbol: the bytes of a symbol of the text
    uint64_t position_bits = 0; // compressed: bits per SA value, block start and LCP block base: bit_width(n)
    uint64_t run_block_offset_bits = 0; // compressed: bits per offset in section_run_block_index
    uint64_t lcp_block_offset_bits = 0; // compressed: bits per offset in section_lcp_block_index
};

static_assert(sizeof(run_data_fields) == 13 * sizeof(uint64_t));

// a section that a caller adds to the runs file (see params::extra_sections)
struct extra_section {
    uint64_t id; // its identifier (at least first_extra_section)
    std::function<void(std::ostream&)> write; // writes the bytes of the section
};

// parameters of the construction
struct params {
    uint64_t window = 10; // length w of a trigger string
    uint64_t modulus = 100; // a window is a trigger string with probability 1/modulus
    uint16_t num_threads = 1; // number of threads
    uint8_t position_bytes = 0; // bytes per text position in the arrays of the construction, in [4..6] (0: the
                                // fewest that hold |T| + window)
    file_format format = file_format::plain; // the format of the runs file
    // number of runs per block (compressed) or per record sample (plain); a multiple of 8
    uint64_t runs_per_block = 1024;
    uint64_t runs_per_symbol_sample = 64; // compressed: runs per sample of the symbol codes, in [1..runs_per_block]
    bool write_sa_samples = false; // write the SA values at the starts and at the ends of the BWT runs
    bool write_lcp_samples = false; // write the LCP values at the starts of the BWT runs
    bool write_bwt = false; // write the BWT with one symbol per position
    bool write_sa = false; // write the suffix array
    bool write_lcp = false; // write the LCP array
    // further sections of the runs file, written after the sections of build (distinct identifiers)
    std::vector<extra_section> extra_sections;
    bool log = false; // print the phases and their running times
};

// sizes of the parsing and of the output
struct result {
    uint64_t n_with_terminator = 0; // length of the text T$ (including the terminator)
    uint64_t r = 0; // number of runs in the BWT of T$
    uint64_t num_phrases = 0; // number of phrases of the prefix-free parsing
    uint64_t num_phrase_runs = 0; // number of maximal runs of equal consecutive phrases
    uint64_t num_words = 0; // number of distinct phrases (dictionary words)
    uint64_t dictionary_size = 0; // total length of the dictionary words (including one separator per word)
};

// file extensions of the output (appended to the prefix passed to build)
/* the runs: the sections (see run_data_section and params::extra_sections), then the table of contents: one
   run_data_entry per section and their number (64 bits) */
inline constexpr const char* ext_runs = ".bwt-runs";
inline constexpr const char* ext_bwt = ".bwt"; // BWT[0..n-1] (sizeof(sym_t) bytes per symbol)
inline constexpr const char* ext_sa = ".sa"; // SA[0..n-1] (byte_width(n) bytes per value)
inline constexpr const char* ext_lcp = ".lcp"; // LCP[0..n-1] (byte_width(n) bytes per value)
// temporary: the suffix array of the dictionary (see dictionary_sa)
inline constexpr const char* ext_dictionary_sa = ".dsa";
// temporary: the ends of the phrase runs and the symbols before the phrases after them (see phrase_arrays)
inline constexpr const char* ext_phrase_run_ends = ".run-ends";
inline constexpr const char* ext_phrase_run_symbols = ".run-symbols";
// temporary: the inverse suffix array of the parse of the phrase runs (see build_parse_sa)
inline constexpr const char* ext_parse_isa = ".isa";
// temporary: the records of the plain format, from which the compressed format is computed
inline constexpr const char* ext_temporary_records = ".records";

/**
 * @brief the layout of a record of the plain format (see section_records): the offsets of its fields in bytes
 */
struct record_layout {
    uint8_t position_bytes = 0; // bytes per run start, SA value and LCP value
    uint8_t symbol_bytes = 0; // bytes per run symbol
    uint8_t symbol = 0; // offset of the symbol (the start is at offset 0)
    uint8_t sa_start = 0, sa_end = 0; // offsets of the SA values at the start and at the end of the run
    uint8_t lcp = 0; // offset of the LCP value at the start of the run
    uint8_t bytes = 0; // bytes per record

    record_layout() = default;

    record_layout(uint8_t position_bytes, uint8_t symbol_bytes, bool sa_samples, bool lcp_samples)
        : position_bytes(position_bytes), symbol_bytes(symbol_bytes), symbol(position_bytes),
          sa_start(position_bytes + symbol_bytes), sa_end(sa_start + sa_samples * position_bytes),
          lcp(sa_end + sa_samples * position_bytes), bytes(lcp + lcp_samples * position_bytes)
    { }

    explicit record_layout(const run_data_fields& f)
        : record_layout(f.position_bytes, f.symbol_bytes, f.sa_samples != 0, f.lcp_samples != 0)
    { }
};

} // namespace pfp_rlbwt

#endif // PFP_RLBWT_PFP_RLBWT_TYPES_HPP
