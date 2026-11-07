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
 * @file types.hpp
 * @brief the parameters and the result of the construction
 */

#ifndef PFP_RLBWT_PFP_RLBWT_TYPES_HPP
#define PFP_RLBWT_PFP_RLBWT_TYPES_HPP

#include <cstdint>

#include <rlbwt/types.hpp>

namespace pfp_rlbwt {

// the text sources, the runs file and its readers of string-utils (rlbwt/)
using namespace ::rlbwt;

// parameters of the construction and of the output (see output_params)
struct params : output_params {
    uint64_t window = 10; // length w of a trigger string
    uint64_t modulus = 100; // a window is a trigger string with probability 1/modulus
    uint16_t num_threads = 1; // number of threads
    uint8_t position_bytes = 0; // bytes per text position in the arrays of the construction, in [4..6] (0: the
                                // fewest that hold |T| + window)
    /* recursive: the parse is suffix-sorted with the prefix-free parsing of the parse (whose window and modulus are
       recursion_window and recursion_modulus), if it has at least recursion_min_phrase_runs phrase runs, at most half
       as many as the text has symbols and its prefix-free parsing is at most half as large as it; and so on for the
       parse of the parse, so the number of levels depends on the text */
    bool recursive = true;
    uint64_t recursion_window = 6;
    uint64_t recursion_modulus = 20;
    uint64_t recursion_min_phrase_runs = uint64_t { 1 } << 20;
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

// temporary: the suffix array of the dictionary (see dictionary_sa)
inline constexpr const char* ext_dictionary_sa = ".dsa";
// temporary: the ends of the phrase runs (see phrase_arrays)
inline constexpr const char* ext_phrase_run_ends = ".run-ends";
// temporary: the inverse suffix array of the parse of the phrase runs (see build_parse_sa)
inline constexpr const char* ext_parse_isa = ".isa";

} // namespace pfp_rlbwt

#endif // PFP_RLBWT_PFP_RLBWT_TYPES_HPP
