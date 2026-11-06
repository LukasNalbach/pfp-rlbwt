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

// the construction of Psi as a move data structure (see psi_mds), which does not depend on the symbol type

#include <array>

#include <omp.h>

#include <pfp_rlbwt/construction/suffix_arrays.hpp>

namespace pfp_rlbwt {

psi_mds::psi_mds(mds_pairs&& pairs, uint64_t num_ranks, const interleaved_bit_aligned_vectors& words,
    const interleaved_bit_aligned_vectors& advances, uint16_t threads)
{
    // the balancing gives each piece of an interval the word and the advance of the interval
    const std::array<mds_pair_row, 2> pair_rows { mds_pair_row { .values = &words },
        mds_pair_row { .values = &advances } };
    const mds_params prm { .num_threads = threads, .pair_rows = pair_rows.data() };
    // the construction sets the number of threads of OpenMP (omp_set_num_threads), which is restored afterwards
    const int omp_threads = omp_get_max_threads();
    const move_data_structure<CACHE_LINES, UNI> M(std::move(pairs), num_ranks, prm,
        { uint8_t(words.width_of(0)), uint8_t(advances.width_of(0)) });
    omp_set_num_threads(omp_threads);
    k = M.num_intervals();
    rows.resize(k + 1);
    advise_huge_pages(rows.data(), rows.size() * sizeof(row));
    starts = run_start_bits(num_ranks, threads);

    // row x: p(x), q(x), idx(x) with the low bits of the word, and the advance with the high bits of the word
    #pragma omp parallel for num_threads(threads)
    for (uint64_t x = 0; x < k; x++) {
        const uint64_t p_x = M.p(x), word = M.row<0>(x);
        starts.mark(p_x);
        rows[x] = { p_x, M.q(x), M.idx(x) | (word << value_bits),
            M.row<1>(x) | ((word >> word_part_bits) << value_bits) };
    }

    // p(k) = J + 1
    rows[k] = { num_ranks, 0, 0, 0 };
    starts.finish();
}

} // namespace pfp_rlbwt
