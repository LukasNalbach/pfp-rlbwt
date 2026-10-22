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

#ifndef PFP_RLBWT_PFP_RLBWT_OUTPUT_READERS_TPP
#define PFP_RLBWT_PFP_RLBWT_OUTPUT_READERS_TPP

#include <pfp_rlbwt/output/readers.hpp>

#include <cstring>
#include <stdexcept>

#include <omp.h>

namespace pfp_rlbwt {

// ############################# direct_block_reader #############################

inline void direct_block_reader::open(const std::string& file_path)
{
    close();
    path = file_path;
    fd = open_direct_or_buffered(path, platform::file_access::read, direct);
    if (fd < 0) throw std::runtime_error("pfp_rlbwt: cannot read " + path);
    if (!direct) advise_uncached(fd);
}

inline void direct_block_reader::close()
{
    if (fd >= 0) platform::close(fd);
    if (buffer != nullptr) platform::free_aligned(buffer);
    fd = -1;
    buffer = nullptr;
    capacity = buffer_beg = buffer_end = 0;
}

inline void direct_block_reader::fill(uint64_t beg, uint64_t end, uint64_t min_end)
{
    const uint64_t read_beg = direct ? (beg / sector) * sector : beg;
    const uint64_t read_end = direct ? div_ceil(end, sector) * sector : end;
    const uint64_t len = read_end - read_beg;

    if (capacity < len + read_padding) {
        if (buffer != nullptr) platform::free_aligned(buffer);
        capacity = div_ceil(len + read_padding, sector) * sector;
        buffer = (uint8_t*) platform::allocate_aligned_counted(capacity, sector);
        if (buffer == nullptr) throw std::bad_alloc();
    }

    const uint64_t num_read = pread_fully(fd, buffer, len, read_beg);
    reads++;
    if (!direct) platform::drop_cached(fd, read_beg, num_read);
    buffer_beg = read_beg;
    buffer_end = read_beg + num_read;
    if (buffer_end < min_end) throw std::runtime_error("pfp_rlbwt: cannot read " + path + " at " + std::to_string(beg));
}

// ############################# plain_disk_reader #############################

inline run_record plain_disk_reader::record(uint64_t i)
{
    const bool last = i + 1 == file->r();
    const uint8_t* p = in.read(at(i, 0), last ? at(i, layout.bytes) : at(i + 1, layout.symbol));
    run_record rec;
    rec.start = read_bytes(p, layout.position_bytes);
    rec.symbol = read_bytes(p + layout.symbol, layout.symbol_bytes);
    rec.length = (last ? file->n() : read_bytes(p + layout.bytes, layout.position_bytes)) - rec.start;

    if (layout.sa_end != layout.sa_start) {
        rec.sa_start = read_bytes(p + layout.sa_start, layout.position_bytes);
        rec.sa_end = read_bytes(p + layout.sa_end, layout.position_bytes);
    }

    if (layout.bytes != layout.lcp) rec.lcp = read_bytes(p + layout.lcp, layout.position_bytes);
    return rec;
}

inline uint64_t plain_disk_reader::run_of(uint64_t j)
{
    // the records of the block of runs whose sampled start is the last one <= j
    const uint64_t b = file->block_of_position(j);
    const uint64_t first = b * file->fields().runs_per_block;
    const uint64_t m = file->runs_in_block(b);
    const uint8_t* p = in.read(at(first, 0), at(first + m - 1, layout.symbol));
    uint64_t lo = 0, hi = m; // start(lo) <= j < start(hi)

    while (hi - lo > 1) {
        const uint64_t mid = lo + (hi - lo) / 2;
        if (read_bytes(p + mid * layout.bytes, layout.position_bytes) <= j) lo = mid;
        else hi = mid;
    }

    return first + lo;
}

// ############################# compressed_run_stream #############################

inline compressed_run_stream::compressed_run_stream(const run_data_file& file, uint64_t first, uint64_t end,
    bool read_sa, bool read_lcp)
    : file(&file), runs_in(file.file()), with_sa(read_sa && file.has_sa_samples()),
      with_lcp(read_lcp && file.has_lcp_samples()), runs_offset(file.section(section_run_blocks).offset),
      block(file.num_blocks()), next_run(first), end_run(end)
{
    if (with_sa) {
        sa_starts_in.open(file.file());
        sa_ends_in.open(file.file());
        sa_starts_offset = file.section(section_sa_run_starts).offset;
        sa_ends_offset = file.section(section_sa_run_ends).offset;
    }

    if (with_lcp) {
        lcp_in.open(file.file());
        lcp_offset = file.section(section_lcp_blocks).offset;
    }
}

inline void compressed_run_stream::load_block(uint64_t b)
{
    const run_data_file& f = *file;
    const uint64_t m = f.runs_in_block(b);
    const uint64_t first = f.block_start(b), next = f.block_start(b + 1);
    const uint8_t* p = runs_in.read_ahead(runs_offset + f.run_block_offset(b), runs_offset + f.run_block_offset(b + 1));
    resize_no_init(starts, m + 1);
    resize_no_init(symbols, m);
    uint64_t k = 0;

    run_block(p, m, next - first, f.fields().runs_per_symbol_sample, f.symbols()).for_each([&](uint64_t x, uint64_t c) {
        starts[k] = first + x;
        symbols[k++] = c;
    });

    starts[m] = next;

    if (with_lcp) {
        const uint8_t* q = lcp_in.read_ahead(lcp_offset + f.lcp_block_offset(b), lcp_offset + f.lcp_block_offset(b + 1));
        resize_no_init(lcps, m);
        k = 0;
        lcp_block(q, m, uint8_t(f.fields().position_bits)).for_each([&](uint64_t v) { lcps[k++] = v; });
    }

    block = b;
}

// ############################# compressed sources #############################

inline compressed_disk_source::compressed_disk_source(const run_data_file& file) : f(&file), in(file.file())
{
    const uint64_t ids[num_fetched_sections] = { section_run_blocks, section_sa_run_starts, section_sa_run_ends,
        section_lcp_blocks };

    for (uint8_t s = 0; s < num_fetched_sections; s++) {
        if (file.contains(ids[s])) offsets[s] = file.section(ids[s]).offset;
    }
}

inline compressed_memory_source::compressed_memory_source(const run_data_file& file, uint16_t threads) : f(file)
{
    const uint64_t ids[num_fetched_sections] = { section_run_blocks, section_sa_run_starts, section_sa_run_ends,
        section_lcp_blocks };
    constexpr uint64_t chunk = direct_block_reader::stream_bytes;

    for (uint8_t s = 0; s < num_fetched_sections; s++) {
        if (!file.contains(ids[s])) continue;
        const run_data_entry& e = file.section(ids[s]);
        no_init_vector<uint8_t>& dst = bytes[s];
        resize_no_init(dst, e.size + read_padding);
        std::memset(dst.data() + e.size, 0, read_padding);
        const uint64_t chunks = div_ceil(e.size, chunk);

        #pragma omp parallel num_threads(std::max<uint64_t>(1, std::min<uint64_t>(threads, chunks)))
        {
            direct_block_reader in(file.file());

            #pragma omp for schedule(dynamic, 1)
            for (uint64_t c = 0; c < chunks; c++) {
                const uint64_t beg = c * chunk, len = std::min(chunk, e.size - beg);
                std::memcpy(dst.data() + beg, in.read(e.offset + beg, e.offset + beg + len), len);
            }
        }
    }
}

inline uint64_t compressed_memory_source::size_in_bytes() const
{
    uint64_t size = f.ram_bytes();
    for (const no_init_vector<uint8_t>& b : bytes) size += b.size();
    return size;
}

// ############################# packed_run_data #############################

inline packed_run_data::packed_run_data(const run_data_file& file, uint16_t threads) : n(file.n())
{
    const run_data_fields& f = file.fields();
    const uint8_t position_bits = std::bit_width(f.n);
    const uint8_t sa_bits = f.sa_samples != 0 ? position_bits : 0;
    runs = interleaved_bit_aligned_vectors({ position_bits, uint8_t(std::max<uint64_t>(1, std::bit_width(f.max_symbol))),
        sa_bits, sa_bits, uint8_t(f.lcp_samples != 0 ? position_bits : 0), 0, 0, 0 });
    const uint64_t r = f.r;
    runs.resize_no_init(r, threads);

    /* the threads read ranges of whole multiples of 64 runs, which start at whole 64-bit words; the runs that start
       less than 64 bits before the end of a range are set atomically */
    constexpr uint64_t word = std::numeric_limits<uint64_t>::digits;
    const uint64_t units = div_ceil(r, word);
    const uint64_t t = std::max<uint64_t>(1, std::min<uint64_t>(threads, units));
    const uint64_t guard = div_ceil(word, runs.bits_per_entry()) + 1;

    #pragma omp parallel for num_threads(t) schedule(static, 1)
    for (uint64_t i_t = 0; i_t < t; i_t++) {
        const uint64_t beg = std::min(r, ((i_t * units) / t) * word);
        const uint64_t end = std::min(r, (((i_t + 1) * units) / t) * word);
        run_stream in(file, beg, end);
        run_record rec;

        for (uint64_t i = beg; i < end; i++) {
            in.next(rec);

            if (i + guard < end) {
                runs.set<start_field>(i, rec.start);
                runs.set<symbol_field>(i, rec.symbol);
                if (sa_bits != 0) runs.set<sa_start_field>(i, rec.sa_start);
                if (sa_bits != 0) runs.set<sa_end_field>(i, rec.sa_end);
                if (f.lcp_samples != 0) runs.set<lcp_field>(i, rec.lcp);
            } else {
                runs.set_parallel<start_field>(i, rec.start);
                runs.set_parallel<symbol_field>(i, rec.symbol);
                if (sa_bits != 0) runs.set_parallel<sa_start_field>(i, rec.sa_start);
                if (sa_bits != 0) runs.set_parallel<sa_end_field>(i, rec.sa_end);
                if (f.lcp_samples != 0) runs.set_parallel<lcp_field>(i, rec.lcp);
            }
        }
    }
}

} // namespace pfp_rlbwt

#endif // PFP_RLBWT_PFP_RLBWT_OUTPUT_READERS_TPP
