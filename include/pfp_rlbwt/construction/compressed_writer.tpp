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

#ifndef PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_COMPRESSED_WRITER_TPP
#define PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_COMPRESSED_WRITER_TPP

#include <pfp_rlbwt/construction/compressed_writer.hpp>

#include <atomic>
#include <exception>
#include <filesystem>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>

#include <omp.h>

namespace pfp_rlbwt {

inline void write_compressed_runs(const std::string& records_name, const std::string& file_name, run_data_fields f,
    std::vector<std::pair<uint64_t, uint64_t>> frequencies, const params& prm)
{
    const symbol_code code(std::move(frequencies));
    const record_layout layout(f);
    const uint64_t r = f.r, R = f.runs_per_block, period = f.runs_per_symbol_sample;
    const uint8_t position_bits = uint8_t(f.position_bits);
    const bool with_sa = f.sa_samples != 0, with_lcp = f.lcp_samples != 0;
    const uint64_t num_blocks = div_ceil(r, R);
    const uint64_t sa_bytes = with_sa ? div_ceil<uint64_t>(r * position_bits, byte_bits) : 0;
    // the SA values are stored first, then the blocks of runs (the sizes of the SA values are known)
    const uint64_t runs_offset = 2 * sa_bytes;

    // the LCP blocks are written to a temporary file and appended after the blocks of runs
    const std::string lcp_name = file_name + "." + std::to_string(section_lcp_blocks);
    struct remover {
        const std::string& name;
        ~remover() { std::error_code ec; std::filesystem::remove(name, ec); }
    } lcp_remover { lcp_name };

    parallel_output_file out(file_name);
    std::unique_ptr<parallel_output_file> lcp_out;
    if (with_lcp) lcp_out = std::make_unique<parallel_output_file>(lcp_name);

    // [0..num_blocks-1] the start of the first run of each block, the offsets of its block of runs and of its LCP block
    no_init_vector<uint64_t> block_starts, run_offsets, lcp_offsets;
    resize_no_init(block_starts, num_blocks);
    resize_no_init(run_offsets, num_blocks);
    resize_no_init(lcp_offsets, with_lcp ? num_blocks : 0);

    const uint16_t p = prm.num_threads;
    // the number of blocks a thread encodes in a batch: their records have about stream_bytes bytes
    const uint64_t batch_blocks = std::max<uint64_t>(1, direct_block_reader::stream_bytes / (R * layout.bytes));
    std::vector<encoded_blocks> encoded(p);
    uint64_t runs_size = 0, lcps_size = 0; // the bytes of the blocks of runs and of the LCP blocks written so far
    std::exception_ptr error;
    std::atomic<bool> failed = false; // true <=> a thread has stored an exception in error
    std::mutex error_lock;

    auto store_error = [&] {
        std::lock_guard<std::mutex> guard(error_lock);
        if (!error) error = std::current_exception();
        failed = true;
    };

    #pragma omp parallel num_threads(p)
    {
        const uint64_t t = omp_get_thread_num();
        encoded_blocks& e = encoded[t];
        direct_block_reader in;
        // [0..R] the values of the runs of a block; starts[m] = the start of the next block
        no_init_vector<uint64_t> starts, symbols, sa_starts, sa_ends, lcps;
        resize_no_init(starts, R + 1);
        resize_no_init(symbols, R);
        resize_no_init(sa_starts, R);
        resize_no_init(sa_ends, R);
        resize_no_init(lcps, R);

        try {
            in.open(records_name);
        } catch (...) {
            store_error();
        }

        for (uint64_t batch = 0; batch < num_blocks; batch += p * batch_blocks) {
            const uint64_t first = std::min(num_blocks, batch + t * batch_blocks);
            const uint64_t last = std::min(num_blocks, first + batch_blocks);
            e.clear();

            try {
                for (uint64_t b = first; b < last && !failed; b++) {
                    const uint64_t i0 = b * R, m = std::min(R, r - i0);

                    for (uint64_t k = 0; k <= m && i0 + k < r; k++) {
                        const uint64_t beg = (i0 + k) * layout.bytes;
                        const uint8_t* rec = in.read_ahead(beg, beg + (k < m ? layout.bytes : layout.position_bytes));
                        starts[k] = read_bytes(rec, layout.position_bytes);
                        if (k == m) break;
                        symbols[k] = read_bytes(rec + layout.symbol, layout.symbol_bytes);

                        if (with_sa) {
                            sa_starts[k] = read_bytes(rec + layout.sa_start, layout.position_bytes);
                            sa_ends[k] = read_bytes(rec + layout.sa_end, layout.position_bytes);
                        }

                        if (with_lcp) lcps[k] = read_bytes(rec + layout.lcp, layout.position_bytes);
                    }

                    if (i0 + m == r) starts[m] = f.n;
                    block_starts[b] = starts[0];
                    for (uint64_t k = m + 1; k-- > 0;) starts[k] -= starts[0];
                    uint64_t before = e.runs.bytes();
                    run_block::encode(e.runs, starts.data(), symbols.data(), m, period, code);
                    e.run_sizes.push_back(e.runs.bytes() - before);

                    if (with_sa) {
                        for (uint64_t k = 0; k < m; k++) e.sa_starts.put(sa_starts[k], position_bits);
                        for (uint64_t k = 0; k < m; k++) e.sa_ends.put(sa_ends[k], position_bits);
                    }

                    if (with_lcp) {
                        before = e.lcps.bytes();
                        lcp_block::encode(e.lcps, lcps.data(), m, position_bits);
                        e.lcp_sizes.push_back(e.lcps.bytes() - before);
                    }
                }

                e.sa_starts.pad_to_byte();
                e.sa_ends.pad_to_byte();
            } catch (...) {
                store_error();
            }

            #pragma omp barrier
            #pragma omp single
            {
                for (encoded_blocks& x : encoded) {
                    x.runs_offset = runs_size;
                    runs_size += x.runs.bytes();
                    x.lcps_offset = lcps_size;
                    lcps_size += x.lcps.bytes();
                }
            }

            try {
                if (!failed && first < last) {
                    uint64_t run_offset = e.runs_offset, lcp_offset = e.lcps_offset;

                    for (uint64_t b = first; b < last; b++) {
                        run_offsets[b] = run_offset;
                        run_offset += e.run_sizes[b - first];

                        if (with_lcp) {
                            lcp_offsets[b] = lcp_offset;
                            lcp_offset += e.lcp_sizes[b - first];
                        }
                    }

                    out.write(e.runs.data(), e.runs.bytes(), runs_offset + e.runs_offset);
                    if (with_lcp) lcp_out->write(e.lcps.data(), e.lcps.bytes(), e.lcps_offset);

                    if (with_sa) {
                        // a block has a multiple of 8 runs, so the values of a range of blocks start at a byte
                        const uint64_t offset = (first * R * position_bits) / byte_bits;
                        out.write(e.sa_starts.data(), e.sa_starts.bytes(), offset);
                        out.write(e.sa_ends.data(), e.sa_ends.bytes(), sa_bytes + offset);
                    }
                }
            } catch (...) {
                store_error();
            }
        }
    }

    if (error) std::rethrow_exception(error);
    encoded.clear();

    // the table of contents and the offset after the last section
    std::vector<run_data_entry> contents;
    uint64_t offset = 0;

    // appends a section to the file
    auto append = [&](uint64_t id, const void* data, uint64_t size) {
        if (size > 0) out.write(data, size, offset);
        contents.push_back({ .id = id, .offset = offset, .size = size });
        offset += size;
    };

    if (with_sa) {
        contents.push_back({ .id = section_sa_run_starts, .offset = 0, .size = sa_bytes });
        contents.push_back({ .id = section_sa_run_ends, .offset = sa_bytes, .size = sa_bytes });
    }

    contents.push_back({ .id = section_run_blocks, .offset = runs_offset, .size = runs_size });
    offset = runs_offset + runs_size;

    if (with_lcp) {
        lcp_out->close();
        contents.push_back({ .id = section_lcp_blocks, .offset = offset, .size = lcps_size });
        direct_block_reader in(lcp_name);

        for (uint64_t pos = 0; pos < lcps_size; pos += direct_block_reader::stream_bytes) {
            const uint64_t len = std::min(direct_block_reader::stream_bytes, lcps_size - pos);
            out.write(in.read_ahead(pos, pos + len), len, offset + pos);
        }

        offset += lcps_size;
    }

    // the block indexes (see section_run_block_index and section_lcp_block_index)
    f.run_block_offset_bits = std::bit_width(runs_size);
    f.lcp_block_offset_bits = with_lcp ? std::bit_width(lcps_size) : 0;
    bit_buffer index;

    for (uint64_t b = 0; b < num_blocks; b++) {
        index.put(run_offsets[b], uint8_t(f.run_block_offset_bits));
        index.put(block_starts[b], position_bits);
    }

    append(section_run_block_index, index.data(), index.bytes());

    if (with_lcp) {
        index.clear();
        for (uint64_t b = 0; b < num_blocks; b++) index.put(lcp_offsets[b], uint8_t(f.lcp_block_offset_bits));
        append(section_lcp_block_index, index.data(), index.bytes());
    }

    const std::string code_bytes = code.serialize(uint8_t(f.symbol_bytes));
    append(section_symbol_code, code_bytes.data(), code_bytes.size());
    append(section_fields, &f, sizeof(f));

    for (const extra_section& extra : prm.extra_sections) {
        std::ostringstream bytes;
        extra.write(bytes);
        const std::string s = bytes.str();
        append(extra.id, s.data(), s.size());
    }

    const uint64_t num_sections = contents.size();
    out.write(contents.data(), num_sections * sizeof(run_data_entry), offset);
    out.write(&num_sections, sizeof(num_sections), offset + num_sections * sizeof(run_data_entry));
    out.close();
}

} // namespace pfp_rlbwt

#endif // PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_COMPRESSED_WRITER_TPP
