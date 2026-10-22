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
 * @file readers.hpp
 * @brief readers of a runs file (see ext_runs): streamed readers (sequentially from any run), random-access readers on
 *        the disk and two representations in RAM (bit-packed, and the compressed format)
 */

#ifndef PFP_RLBWT_PFP_RLBWT_OUTPUT_READERS_HPP
#define PFP_RLBWT_PFP_RLBWT_OUTPUT_READERS_HPP

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <data_structures/interleaved_bit_aligned_vectors.hpp>
#include <misc/files.hpp>
#include <misc/memory.hpp>
#include <pfp_rlbwt/output/blocks.hpp>
#include <pfp_rlbwt/output/run_data_file.hpp>

namespace pfp_rlbwt {

// a run of the BWT of T$ with the values stored for it (0 for values the file does not store)
struct run_record {
    uint64_t start = 0; // the BWT position of its first symbol
    uint64_t length = 0; // its length
    uint64_t symbol = 0; // its symbol
    uint64_t sa_start = 0; // the SA value at its start
    uint64_t sa_end = 0; // the SA value at its end
    uint64_t lcp = 0; // the LCP value at its start (0 for run 0)
};

// the start and the symbol of a run
struct run_head {
    uint64_t start = 0;
    uint64_t symbol = 0;
};

/**
 * @brief reads ranges of a file with O_DIRECT (if the file system does not support it, through the page cache, whose
 *        pages are dropped): read() reads exactly the sectors of a range with one pread; read_ahead() returns the range
 *        from the buffer if it holds it, and otherwise reads at least stream_bytes bytes from the range on with one
 *        pread. A returned range is followed by at least read_padding readable bytes.
 */
class direct_block_reader {
public:
    static constexpr uint64_t sector = direct_io_alignment;
    static constexpr uint64_t stream_bytes = direct_reader::block; // bytes read_ahead reads at once (at least)

private:
    std::string path; // path of the file
    int fd = -1; // file descriptor
    bool direct = false; // true <=> the file is opened with O_DIRECT
    uint8_t* buffer = nullptr; // capacity bytes, aligned to sector
    uint64_t capacity = 0;
    uint64_t buffer_beg = 0, buffer_end = 0; // the bytes [buffer_beg, buffer_end) of the file are in the buffer
    uint64_t reads = 0; // number of preads

    // reads the sectors containing the bytes [beg, end) into the buffer (at least the bytes [beg, min_end))
    void fill(uint64_t beg, uint64_t end, uint64_t min_end);

public:
    direct_block_reader() = default;
    explicit direct_block_reader(const std::string& file_path) { open(file_path); }
    ~direct_block_reader() { close(); }

    direct_block_reader(const direct_block_reader&) = delete;
    direct_block_reader& operator=(const direct_block_reader&) = delete;

    direct_block_reader(direct_block_reader&& other) noexcept { *this = std::move(other); }

    direct_block_reader& operator=(direct_block_reader&& other) noexcept
    {
        std::swap(path, other.path);
        std::swap(fd, other.fd);
        std::swap(direct, other.direct);
        std::swap(buffer, other.buffer);
        std::swap(capacity, other.capacity);
        std::swap(buffer_beg, other.buffer_beg);
        std::swap(buffer_end, other.buffer_end);
        std::swap(reads, other.reads);
        return *this;
    }

    // opens a file; throws if it cannot be opened
    void open(const std::string& file_path);

    // closes the file and frees the buffer
    void close();

    // reads the bytes [beg, end) of the file (end > beg) with one pread and returns a pointer to them
    const uint8_t* read(uint64_t beg, uint64_t end)
    {
        fill(beg, end, end);
        return buffer + (beg - buffer_beg);
    }

    // returns a pointer to the bytes [beg, end) of the file (end > beg), read with at most one pread
    const uint8_t* read_ahead(uint64_t beg, uint64_t end)
    {
        if (beg < buffer_beg || end > buffer_end) fill(beg, std::max(end, beg + stream_bytes), end);
        return buffer + (beg - buffer_beg);
    }

    // returns the number of preads so far
    uint64_t num_reads() const { return reads; }
};

/**
 * @brief reads the runs of a file in the plain format in order, from any run on
 */
class plain_run_stream {
private:
    const run_data_file* file = nullptr;
    record_layout layout;
    direct_block_reader in;
    uint64_t offset = 0; // offset of section_records in the file
    uint64_t next_run = 0, end_run = 0; // the runs [next_run, end_run) are read

public:
    /**
     * @brief a stream of the runs [first, end)
     * @param file the runs file (it must stay valid)
     * @param first the first run
     * @param end the run after the last one (at most r)
     */
    plain_run_stream(const run_data_file& file, uint64_t first, uint64_t end)
        : file(&file), layout(file.fields()), in(file.file()), offset(file.section(section_records).offset),
          next_run(first), end_run(end)
    { }

    // returns the index of the next run
    uint64_t position() const { return next_run; }

    // returns whether all runs have been read
    bool done() const { return next_run >= end_run; }

    // reads the next run
    void next(run_record& rec)
    {
        const uint64_t i = next_run++;
        const bool last = i + 1 == file->r();
        const uint64_t beg = offset + i * layout.bytes;
        const uint8_t* p = in.read_ahead(beg, beg + layout.bytes + (last ? 0 : layout.position_bytes));
        rec.start = read_bytes(p, layout.position_bytes);
        rec.symbol = read_bytes(p + layout.symbol, layout.symbol_bytes);
        rec.length = (last ? file->n() : read_bytes(p + layout.bytes, layout.position_bytes)) - rec.start;

        if (layout.sa_end != layout.sa_start) {
            rec.sa_start = read_bytes(p + layout.sa_start, layout.position_bytes);
            rec.sa_end = read_bytes(p + layout.sa_end, layout.position_bytes);
        }

        if (layout.bytes != layout.lcp) rec.lcp = read_bytes(p + layout.lcp, layout.position_bytes);
    }
};

/**
 * @brief reads the runs of a file in the compressed format in order, from any run on; each block of runs and of LCP
 *        values is decoded completely
 */
class compressed_run_stream {
private:
    const run_data_file* file = nullptr;
    direct_block_reader runs_in, sa_starts_in, sa_ends_in, lcp_in; // a reader of each section
    bool with_sa = false, with_lcp = false; // true <=> the SA (LCP) values are read
    uint64_t runs_offset = 0, sa_starts_offset = 0, sa_ends_offset = 0, lcp_offset = 0; // offsets of the sections
    // the current block: [0..m-1] the start and the symbol of each run, starts[m] = the start of the next block
    no_init_vector<uint64_t> starts, symbols;
    no_init_vector<uint64_t> lcps; // [0..m-1] the LCP values of the current block
    uint64_t block = 0; // the current block
    uint64_t next_run = 0, end_run = 0; // the runs [next_run, end_run) are read

    // decodes block b
    void load_block(uint64_t b);

    // returns the SA value of run i in the section of in at offset
    uint64_t sa_value(direct_block_reader& in, uint64_t offset, uint64_t i)
    {
        const uint8_t bits = uint8_t(file->fields().position_bits);
        const uint64_t bit = i * bits;
        const uint8_t* p = in.read_ahead(offset + bit / byte_bits, offset + div_ceil<uint64_t>(bit + bits, byte_bits));
        return read_bits(p, bit % byte_bits, bits);
    }

public:
    /**
     * @brief a stream of the runs [first, end)
     * @param file the runs file (it must stay valid)
     * @param first the first run
     * @param end the run after the last one (at most r)
     * @param read_sa whether the SA values are read (if the file stores them)
     * @param read_lcp whether the LCP values are read (if the file stores them)
     */
    compressed_run_stream(const run_data_file& file, uint64_t first, uint64_t end, bool read_sa = true,
        bool read_lcp = true);

    // returns the index of the next run
    uint64_t position() const { return next_run; }

    // returns whether all runs have been read
    bool done() const { return next_run >= end_run; }

    // reads the next run
    void next(run_record& rec)
    {
        const uint64_t i = next_run++;
        const uint64_t b = i / file->fields().runs_per_block;
        if (b != block) load_block(b);
        const uint64_t k = i - b * file->fields().runs_per_block;
        rec.start = starts[k];
        rec.length = starts[k + 1] - starts[k];
        rec.symbol = symbols[k];

        if (with_sa) {
            rec.sa_start = sa_value(sa_starts_in, sa_starts_offset, i);
            rec.sa_end = sa_value(sa_ends_in, sa_ends_offset, i);
        }

        if (with_lcp) rec.lcp = lcps[k];
    }
};

/**
 * @brief reads the runs of a file in either format in order, from any run on
 */
class run_stream {
private:
    std::optional<plain_run_stream> plain;
    std::optional<compressed_run_stream> compressed;

public:
    /**
     * @brief a stream of the runs [first, end)
     * @param file the runs file (it must stay valid)
     * @param first the first run
     * @param end the run after the last one (default: r)
     * @param read_sa whether the SA values are read (compressed format; the plain format reads all values)
     * @param read_lcp whether the LCP values are read (compressed format)
     */
    explicit run_stream(const run_data_file& file, uint64_t first = 0,
        uint64_t end = std::numeric_limits<uint64_t>::max(), bool read_sa = true, bool read_lcp = true)
    {
        end = std::min(end, file.r());
        if (file.format() == file_format::plain) plain.emplace(file, first, end);
        else compressed.emplace(file, first, end, read_sa, read_lcp);
    }

    // returns the index of the next run
    uint64_t position() const { return plain ? plain->position() : compressed->position(); }

    // returns whether all runs have been read
    bool done() const { return plain ? plain->done() : compressed->done(); }

    // reads the next run
    void next(run_record& rec)
    {
        if (plain) plain->next(rec);
        else compressed->next(rec);
    }
};

/**
 * @brief random access to the runs of a file in the plain format on the disk; each access reads the file once (with
 *        O_DIRECT): run, length, record, run_of and each sample read one range of the file. One reader per thread.
 */
class plain_disk_reader {
private:
    const run_data_file* file = nullptr;
    record_layout layout;
    direct_block_reader in;
    uint64_t offset = 0; // offset of section_records in the file

    // returns the offset in the file of field of record i
    uint64_t at(uint64_t i, uint8_t field) const { return offset + i * layout.bytes + field; }

public:
    // a reader of a file (it must stay valid)
    explicit plain_disk_reader(const run_data_file& file)
        : file(&file), layout(file.fields()), in(file.file()), offset(file.section(section_records).offset)
    { }

    uint64_t num_runs() const { return file->r(); }
    uint64_t text_length() const { return file->n(); }

    // returns the number of reads of the file so far
    uint64_t num_reads() const { return in.num_reads(); }

    // returns run i with all values the file stores
    run_record record(uint64_t i);

    // returns the start and the symbol of run i
    run_head run(uint64_t i)
    {
        const uint8_t* p = in.read(at(i, 0), at(i, layout.sa_start));
        return { read_bytes(p, layout.position_bytes), read_bytes(p + layout.symbol, layout.symbol_bytes) };
    }

    // returns the length of run i
    uint64_t length(uint64_t i)
    {
        if (i + 1 == file->r()) return file->n() - read_bytes(in.read(at(i, 0), at(i, layout.symbol)), layout.position_bytes);
        const uint8_t* p = in.read(at(i, 0), at(i + 1, layout.symbol));
        return read_bytes(p + layout.bytes, layout.position_bytes) - read_bytes(p, layout.position_bytes);
    }

    // returns the run containing the BWT position j (j < n)
    uint64_t run_of(uint64_t j);

    // returns the SA value at the start of run i
    uint64_t sa_start(uint64_t i) { return read_bytes(in.read(at(i, layout.sa_start), at(i, layout.sa_end)), layout.position_bytes); }

    // returns the SA value at the end of run i
    uint64_t sa_end(uint64_t i) { return read_bytes(in.read(at(i, layout.sa_end), at(i, layout.lcp)), layout.position_bytes); }

    // returns the LCP value at the start of run i
    uint64_t lcp(uint64_t i) { return read_bytes(in.read(at(i, layout.lcp), at(i, layout.bytes)), layout.position_bytes); }
};

/**
 * @brief random access to the runs of a file in the compressed format; the bytes of a section are fetched from a
 *        source (the disk or RAM). run, length and run_of fetch one block of runs, each SA sample fetches its bytes
 *        and lcp fetches one block of LCP values: one read of the file each with a source on the disk.
 * @tparam source_t provides file() (the run_data_file) and get(section, beg, end), which returns a pointer to the bytes
 *         [beg, end) of a section, followed by read_padding readable bytes
 */
// the sections of the compressed format whose bytes a source of compressed_access fetches
enum fetched_section : uint8_t { fetched_run_blocks, fetched_sa_starts, fetched_sa_ends, fetched_lcp_blocks,
    num_fetched_sections };

template <typename source_t>
class compressed_access {
protected:
    source_t source;

    const run_data_file& file() const { return source.file(); }

    // returns the block of runs b
    run_block run_block_at(uint64_t b)
    {
        const run_data_file& f = file();
        const uint8_t* p = source.get(fetched_run_blocks, f.run_block_offset(b), f.run_block_offset(b + 1));
        return run_block(p, f.runs_in_block(b), f.block_start(b + 1) - f.block_start(b),
            f.fields().runs_per_symbol_sample, f.symbols());
    }

    // returns the SA value of run i in the section s
    uint64_t sa_value(fetched_section s, uint64_t i)
    {
        const uint8_t bits = uint8_t(file().fields().position_bits);
        const uint64_t bit = i * bits;
        return read_bits(source.get(s, bit / byte_bits, div_ceil<uint64_t>(bit + bits, byte_bits)), bit % byte_bits, bits);
    }

    // returns the start, the length and the symbol of run i (one fetch)
    run_record run_and_length(uint64_t i)
    {
        const run_data_file& f = file();
        const uint64_t R = f.fields().runs_per_block, b = i / R, k = i - b * R;
        const run_block blk = run_block_at(b);
        const uint64_t start = blk.start(k);
        const uint64_t next = k + 1 < f.runs_in_block(b) ? blk.start(k + 1) : f.block_start(b + 1) - f.block_start(b);
        return { .start = f.block_start(b) + start, .length = next - start, .symbol = blk.symbol(k) };
    }

public:
    template <typename... args_t>
    explicit compressed_access(args_t&&... args) : source(std::forward<args_t>(args)...) { }

    uint64_t num_runs() const { return file().r(); }
    uint64_t text_length() const { return file().n(); }

    // returns the start and the symbol of run i
    run_head run(uint64_t i)
    {
        const uint64_t R = file().fields().runs_per_block, b = i / R;
        const run_block blk = run_block_at(b);
        return { file().block_start(b) + blk.start(i - b * R), blk.symbol(i - b * R) };
    }

    // returns the length of run i
    uint64_t length(uint64_t i) { return run_and_length(i).length; }

    // returns the run containing the BWT position j (j < n)
    uint64_t run_of(uint64_t j)
    {
        const uint64_t b = file().block_of_position(j);
        return b * file().fields().runs_per_block + run_block_at(b).predecessor(j - file().block_start(b));
    }

    // returns the SA value at the start of run i
    uint64_t sa_start(uint64_t i) { return sa_value(fetched_sa_starts, i); }

    // returns the SA value at the end of run i
    uint64_t sa_end(uint64_t i) { return sa_value(fetched_sa_ends, i); }

    // returns the LCP value at the start of run i
    uint64_t lcp(uint64_t i)
    {
        const run_data_file& f = file();
        const uint64_t R = f.fields().runs_per_block, b = i / R;
        const uint8_t* p = source.get(fetched_lcp_blocks, f.lcp_block_offset(b), f.lcp_block_offset(b + 1));
        return lcp_block(p, f.runs_in_block(b), uint8_t(f.fields().position_bits)).value(i - b * R);
    }

    // returns run i with all values the file stores
    run_record record(uint64_t i)
    {
        run_record rec = run_and_length(i);

        if (file().has_sa_samples()) {
            rec.sa_start = sa_start(i);
            rec.sa_end = sa_end(i);
        }

        if (file().has_lcp_samples()) rec.lcp = lcp(i);
        return rec;
    }
};

// a source of compressed_access that reads the sections from the disk (one pread per get)
class compressed_disk_source {
private:
    const run_data_file* f = nullptr;
    direct_block_reader in;
    std::array<uint64_t, num_fetched_sections> offsets {}; // the offsets of the fetched sections in the file

public:
    explicit compressed_disk_source(const run_data_file& file);

    const run_data_file& file() const { return *f; }

    const uint8_t* get(uint8_t section, uint64_t beg, uint64_t end)
    {
        return in.read(offsets[section] + beg, offsets[section] + end);
    }

    uint64_t num_reads() const { return in.num_reads(); }
};

/**
 * @brief random access to the runs of a file in the compressed format on the disk (one reader per thread): run, length
 *        and run_of read one block of runs, sa_start and sa_end read the bytes of one value and lcp reads one block of
 *        LCP values, each with one pread (with O_DIRECT); record reads up to four times
 */
class compressed_disk_reader : public compressed_access<compressed_disk_source> {
public:
    // a reader of a file (it must stay valid)
    explicit compressed_disk_reader(const run_data_file& file) : compressed_access(file) { }

    // returns the number of reads of the file so far
    uint64_t num_reads() const { return source.num_reads(); }
};

// a source of compressed_access that holds the sections in RAM
class compressed_memory_source {
private:
    run_data_file f; // a copy of the runs file
    std::array<no_init_vector<uint8_t>, num_fetched_sections> bytes; // the fetched sections, each followed by read_padding bytes

public:
    compressed_memory_source() = default;

    /**
     * @brief reads the sections of a file into RAM
     * @param file the runs file (in the compressed format)
     * @param threads number of threads that read
     */
    compressed_memory_source(const run_data_file& file, uint16_t threads);

    const run_data_file& file() const { return f; }

    const uint8_t* get(uint8_t section, uint64_t beg, uint64_t) const { return bytes[section].data() + beg; }

    // returns the number of bytes in RAM
    uint64_t size_in_bytes() const;
};

/**
 * @brief the runs of a file in the compressed format, held in RAM in this format (RAM mode "compressed")
 */
class compressed_run_data : public compressed_access<compressed_memory_source> {
public:
    compressed_run_data() = default;

    /**
     * @brief reads a file in the compressed format
     * @param file the runs file
     * @param threads number of threads that read
     */
    compressed_run_data(const run_data_file& file, uint16_t threads = 1) : compressed_access(file, threads) { }

    // returns the number of bytes in RAM
    uint64_t size_in_bytes() const { return source.size_in_bytes(); }
};

/**
 * @brief the runs of a file in either format, held in RAM bit-packed (RAM mode "packed"): per run its start, its
 *        symbol and the values the file stores, interleaved, each with a fixed number of bits
 */
class packed_run_data {
public:
    // the vectors of runs
    enum field : uint8_t { start_field, symbol_field, sa_start_field, sa_end_field, lcp_field };

private:
    interleaved_bit_aligned_vectors runs;
    uint64_t n = 0; // length of T$

public:
    packed_run_data() = default;

    /**
     * @brief reads a file
     * @param file the runs file
     * @param threads number of threads that read
     */
    packed_run_data(const run_data_file& file, uint16_t threads = 1);

    uint64_t num_runs() const { return runs.size(); }
    uint64_t text_length() const { return n; }

    // returns the number of bytes in RAM
    uint64_t size_in_bytes() const { return runs.size_in_bytes(); }

    // returns the start and the symbol of run i
    run_head run(uint64_t i) const { return { runs.get<start_field>(i), runs.get<symbol_field>(i) }; }

    // returns the length of run i
    uint64_t length(uint64_t i) const
    {
        return (i + 1 == runs.size() ? n : runs.get<start_field>(i + 1)) - runs.get<start_field>(i);
    }

    // returns the run containing the BWT position j (j < n)
    uint64_t run_of(uint64_t j) const
    {
        uint64_t lo = 0, hi = runs.size(); // start(lo) <= j < start(hi)

        while (hi - lo > 1) {
            const uint64_t mid = lo + (hi - lo) / 2;
            if (runs.get<start_field>(mid) <= j) lo = mid;
            else hi = mid;
        }

        return lo;
    }

    uint64_t sa_start(uint64_t i) const { return runs.get<sa_start_field>(i); }
    uint64_t sa_end(uint64_t i) const { return runs.get<sa_end_field>(i); }
    uint64_t lcp(uint64_t i) const { return runs.get<lcp_field>(i); }

    // returns run i with all values the file stores (0 for the others)
    run_record record(uint64_t i) const
    {
        return { .start = runs.get<start_field>(i), .length = length(i), .symbol = runs.get<symbol_field>(i),
            .sa_start = runs.get<sa_start_field>(i), .sa_end = runs.get<sa_end_field>(i), .lcp = runs.get<lcp_field>(i) };
    }
};

} // namespace pfp_rlbwt

#include <pfp_rlbwt/output/readers.tpp>

#endif // PFP_RLBWT_PFP_RLBWT_OUTPUT_READERS_HPP
