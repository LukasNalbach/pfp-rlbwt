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

#ifndef PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_COMPRESSED_WRITER_HPP
#define PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_COMPRESSED_WRITER_HPP

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include <pfp_rlbwt/construction/common.hpp>
#include <pfp_rlbwt/output/blocks.hpp>
#include <pfp_rlbwt/output/readers.hpp>
#include <pfp_rlbwt/types.hpp>

namespace pfp_rlbwt {

// the encoded blocks of a range of blocks of one thread (see write_compressed_runs); aligned to a cache line
struct alignas(cache_line_bytes) encoded_blocks {
    bit_buffer runs, lcps, sa_starts, sa_ends; // the blocks of runs, the LCP blocks and the SA values
    std::vector<uint64_t> run_sizes, lcp_sizes; // the size in bytes of each block of runs and of LCP values
    uint64_t runs_offset = 0, lcps_offset = 0; // the offsets of runs and lcps in their sections

    void clear()
    {
        runs.clear();
        lcps.clear();
        sa_starts.clear();
        sa_ends.clear();
        run_sizes.clear();
        lcp_sizes.clear();
    }
};

/**
 * @brief writes the runs file in the compressed format from the records of the plain format, with prm.num_threads
 *        threads. Its sections: section_sa_run_starts, section_sa_run_ends (if sa_samples), section_run_blocks,
 *        section_lcp_blocks (if lcp_samples), section_run_block_index, section_lcp_block_index (if lcp_samples),
 *        section_symbol_code, section_fields, the extra sections and the table of contents. Each thread encodes
 *        consecutive blocks, reading their records with O_DIRECT; the blocks are written with O_DIRECT at their
 *        offsets.
 * @param records_name the file of the records (see record_layout)
 * @param file_name the runs file
 * @param f the fields (the fields of the compressed format are set here)
 * @param frequencies the (symbol, number of runs) pairs of the run symbols, sorted by the symbols
 * @param prm parameters
 */
inline void write_compressed_runs(const std::string& records_name, const std::string& file_name, run_data_fields f,
    std::vector<std::pair<uint64_t, uint64_t>> frequencies, const params& prm);

} // namespace pfp_rlbwt

#include <pfp_rlbwt/construction/compressed_writer.tpp>

#endif // PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_COMPRESSED_WRITER_HPP
