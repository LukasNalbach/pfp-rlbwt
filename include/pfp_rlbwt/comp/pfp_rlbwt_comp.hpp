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
 * @file pfp_rlbwt_comp.hpp
 * @brief compression of a byte text into its run-length encoded BWT and decompression with the move data structure
 */

#ifndef PFP_RLBWT_PFP_RLBWT_COMP_PFP_RLBWT_COMP_HPP
#define PFP_RLBWT_PFP_RLBWT_COMP_PFP_RLBWT_COMP_HPP

#include <omp.h>

#include <cstdint>
#include <string>

#include <pfp_rlbwt/comp/byte_text.hpp>

namespace pfp_rlbwt::comp {

// default of compress_params::max_start_pairs
inline constexpr uint64_t default_max_start_pairs = 4096;

// default of compress_params::runs_per_block
inline constexpr uint64_t default_runs_per_block = uint64_t(1) << 20;

// the number of L1 data cache misses a core can have outstanding
inline constexpr uint64_t outstanding_cache_misses = 16;

// default of decompress_params::num_strands: one prefetched block of M_LF per outstanding cache miss
inline constexpr uint64_t default_num_strands = outstanding_cache_misses;

// parameters of the compression
struct compress_params {
    uint64_t num_threads = omp_get_max_threads(); // number of threads of the BWT construction
    // the file stores k = min(r, max_start_pairs) start pairs, at most one per text range of length n / k
    uint64_t max_start_pairs = default_max_start_pairs;
    // the runs are coded in blocks of at most runs_per_block runs, each with new models, by up to num_threads threads
    uint64_t runs_per_block = default_runs_per_block;
    // compress the text as a FASTA file: the sequences into the run-length encoded BWT, the headers and the line
    // layout separately
    bool fasta = false;
    bool log = false; // print the phases and their running times
};

// sizes of the text and of the parts of the compressed file
struct compress_result {
    uint64_t n = 0; // n, the length of the text
    // the length of the text T whose BWT is stored (the sequence text in the FASTA mode, else n)
    uint64_t sequence_length = 0;
    uint64_t sigma = 0; // the number of distinct bytes in T
    uint64_t r = 0; // the number of runs of the BWT of T$
    uint64_t num_start_pairs = 0; // the number of start pairs
    uint64_t num_run_blocks = 0; // the number of blocks the runs are coded in
    uint64_t num_records = 0; // the number of records of a FASTA file
    // size of the mode, the FASTA fields, the length of T, the alphabet, r, the run of $, the number of start pairs and
    // the size of the coded run symbols in bytes
    uint64_t header_bytes = 0;
    uint64_t start_pair_bytes = 0; // size of the start pairs in bytes
    uint64_t symbol_bytes = 0; // size of the coded run symbols in bytes
    uint64_t length_bytes = 0; // size of the coded run lengths in bytes
    uint64_t fasta_header_bytes = 0; // size of the coded header lines of a FASTA file in bytes
    uint64_t fasta_layout_bytes = 0; // size of the coded line layout of a FASTA file in bytes
    uint64_t file_bytes = 0; // size of the compressed file in bytes
};

// parameters of the decompression
struct decompress_params {
    uint64_t num_threads = omp_get_max_threads(); // number of threads
    uint64_t num_strands = default_num_strands; // number of pieces each thread reverts alternately (see revert)
    bool log = false; // print the phases and their running times
};

/**
 * @brief compresses a text: builds the BWT of T$ with the prefix-free parsing (pfp-rlbwt) and writes its runs and
 *        the start pairs to a file; the temporary files of pfp-rlbwt are stored next to the output file
 * @param text the text
 * @param output_file the compressed file
 * @param params parameters
 * @return the sizes of the text and of the compressed file
 */
inline compress_result compress(const byte_text& text, const std::string& output_file,
    const compress_params& params = {});

/**
 * @brief compresses a text stored in a file (see compress(const byte_text&, ...))
 * @param input_file the file storing the text
 * @param output_file the compressed file
 * @param params parameters
 * @return the sizes of the text and of the compressed file
 */
inline compress_result compress(const std::string& input_file, const std::string& output_file,
    const compress_params& params = {})
{
    return compress(byte_text(input_file), output_file, params);
}

/**
 * @brief decompresses a file written by compress into a file: builds M_LF from the runs and reverts the pieces
 *        of the text between the start pairs with several threads, which write them to the output file (in the
 *        FASTA mode, the pieces of the sequence text are passed in order to the restoration of the FASTA file,
 *        which writes them with the headers and the line ends)
 * @param input_file the compressed file
 * @param output_file the file the text is written to
 * @param params parameters
 * @return n, the length of the text
 */
inline uint64_t decompress(const std::string& input_file, const std::string& output_file,
    const decompress_params& params = {});

/**
 * @brief decompresses a file written by compress into memory (see decompress(..., output_file, ...))
 * @param input_file the compressed file
 * @param params parameters
 * @return the text
 */
inline std::string decompress(const std::string& input_file, const decompress_params& params = {});

} // namespace pfp_rlbwt::comp

#include <pfp_rlbwt/comp/compress.tpp>
#include <pfp_rlbwt/comp/decompress.tpp>

#endif // PFP_RLBWT_PFP_RLBWT_COMP_PFP_RLBWT_COMP_HPP
