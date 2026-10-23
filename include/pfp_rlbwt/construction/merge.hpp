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

#ifndef PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_MERGE_HPP
#define PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_MERGE_HPP

#include <algorithm>
#include <array>
#include <atomic>
#include <condition_variable>
#include <exception>
#include <filesystem>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <misc/sort.hpp>

#include <data_structures/elias_fano_bit_vector.hpp>
#include <misc/files.hpp>
#include <misc/memory.hpp>
#include <omp.h>

#include <pfp_rlbwt/construction/compressed_writer.hpp>
#include <pfp_rlbwt/construction/suffix_arrays.hpp>
#include <pfp_rlbwt/output/bits.hpp>
#include <pfp_rlbwt/types.hpp>

namespace pfp_rlbwt {

// number of chunks of the dictionary suffix array per thread
constexpr uint64_t chunks_per_thread = 64;

// number of chunks per thread that are computed but not yet written, at most
constexpr uint64_t pending_chunks_per_thread = 4;

/* a BWT position: the suffix of S starting there begins with the dictionary suffix
   dict.text[dict_pos..dict_pos+dict_suffix_len-1] of the phrase containing it (if back > 0, dict_pos lies in the word
   of that phrase), continues with the back symbols of the phrases after it in its phrase run (back = 0 for the last
   phrase of a run), and then with the suffix of the parse of the phrase runs whose rank is parse_rank (the one
   starting at the next run) */
struct bwt_pos {
    uint64_t dict_pos;
    uint64_t dict_suffix_len;
    uint64_t parse_rank;
    uint64_t back;
};

/* a parse suffix in bwt_merger::lcp_in_runs: it starts with rem symbols of copies of word after their first w symbols
   (rem = 0: none), and continues with the suffix of the parse of the phrase runs of rank next */
struct run_walk {
    uint64_t word;
    uint64_t rem;
    uint64_t next;
};

// the SA value and the LCP value at the start of a BWT run
using run_start_values = std::pair<uint64_t, uint64_t>;

/* the PLCP array of T$ sampled at the SA values of the starts of the BWT runs: PLCP[x] = values[j] - (x - positions[j])
   for the largest positions[j] <= x (from one text position to the next, PLCP decreases by one unless the BWT
   position of the next one starts a run) */
struct sampled_plcp {
    // marks the sampled positions in [0..n-1] (n = |T$|)
    elias_fano_bit_vector is_sampled;
    // the sampled positions in increasing order: [0..r-1] positions[j] = the j-th smallest SA value at a run start
    uint_array positions;
    // [0..r-1] values[j] = PLCP[positions[j]]
    uint_array values;

    // returns PLCP[x]
    inline uint64_t get(uint64_t x) const
    {
        const uint64_t j = is_sampled.rank_1(x + 1) - 1;
        return values[j] - (x - positions[j]);
    }
};

/**
 * @brief samples the PLCP array at the SA values of the starts of the runs and frees the values of the runs
 * @param runs the SA and LCP values at the starts of the runs
 * @param n length of T$
 * @param threads number of threads
 */
inline sampled_plcp build_sampled_plcp(std::vector<run_start_values>& runs, uint64_t n, uint16_t threads)
{
    const uint64_t r = runs.size();
    radix_sort(runs.begin(), runs.end(), key_of_first {}, threads);
    sampled_plcp plcp;
    plcp.positions = uint_array(r, value_bytes(n));
    plcp.values = uint_array(r, value_bytes(n));

    #pragma omp parallel for num_threads(threads)
    for (uint64_t j = 0; j < r; j++) {
        plcp.positions.set(j, runs[j].first);
        plcp.values.set(j, runs[j].second);
    }

    plcp.is_sampled = elias_fano_bit_vector([&](uint64_t j) { return runs[j].first; }, r, n, threads);
    runs.clear();
    runs.shrink_to_fit();
    return plcp;
}

/* the BWT positions of a chunk of the dictionary suffix array (see merge): their runs (maximal within the chunk), the
   records of the runs inside the chunk and the part of the full BWT, SA and LCP array that is not yet written; aligned
   to a cache line (the threads write different chunks) */
template <typename sym_t>
struct alignas(cache_line_bytes) chunk_runs {
    uint64_t length; // number of BWT positions
    uint64_t num_runs; // number of runs
    sym_t first_sym, last_sym; // symbols of the first and the last run
    uint64_t first_len, last_len; // lengths of the first and the last run
    bwt_pos first_pos, last_pos; // the first and the last BWT position
    uint64_t first_sa_end; // the SA value at the end of the first run (if num_runs > 1 and SA samples are written)
    uint64_t last_start; // the start of the last run relative to the first position of the chunk
    // the SA value and the LCP value at the start of the last run (if num_runs > 1)
    uint64_t last_sa_start, last_lcp;
    /* the records (see record_layout) of the runs 2..num_runs-1, with their starts relative to the first position of
       the chunk */
    no_init_vector<uint8_t> records;
    uint64_t max_inner_length; // the largest length of a run in records
    // the SA and LCP values at the starts of the runs 2..num_runs (only if collected, see collects_run_starts)
    std::vector<run_start_values> run_starts;
    // [0..full_bwt.size()-1] full_bwt[i] = BWT[full_pos+i], not yet written
    std::vector<sym_t> full_bwt;
    // the SA and LCP values of the same positions, position_bytes bytes per value (little endian), not yet written
    std::vector<uint8_t> full_sa, full_lcp;
    uint64_t full_pos; // BWT position of the first value in full_bwt, full_sa and full_lcp

    // empties the chunk; the first BWT position of the chunk is start
    void clear(uint64_t start)
    {
        length = num_runs = first_len = last_len = 0;
        first_sym = last_sym = 0;
        first_sa_end = last_start = last_sa_start = last_lcp = max_inner_length = 0;
        records.clear();
        run_starts.clear();
        full_bwt.clear();
        full_sa.clear();
        full_lcp.clear();
        full_pos = start;
    }
};

/* the ranks of the parse suffixes preceded by a word whose suffix belongs to a group of equal dictionary suffixes
   that are not processed yet: psi[cur..end-1]; member is the index of the dictionary suffix in the group */
template <typename idx_t>
struct rank_range {
    idx_t cur;
    idx_t end;
    uint32_t member;
};

/* the suffixes of P (the parse of the phrases) that start with a word whose suffix belongs to a group of equal
   dictionary suffixes with a repeated word: those with the offsets [cur..end-1] among the suffixes starting with
   word; next_rank is the rank among the suffixes of P of the suffix one phrase after the one with offset cur; rw is
   psa.find_repeated(word); member is the index of the dictionary suffix in the group */
struct repeated_range {
    uint64_t cur;
    uint64_t end;
    uint64_t next_rank;
    uint64_t word;
    const repeated_word_ranks* rw;
    uint64_t member;
};

// the heaps of a thread (see bwt_merger::compute_chunk)
template <typename idx_t>
struct merge_heaps {
    std::vector<rank_range<idx_t>> plain;
    std::vector<repeated_range> repeated;
};

// capacity of the heap of a thread above which the heap is freed after a group of equal dictionary suffixes
constexpr uint64_t kept_heap_entries = uint64_t { 1 } << 12;

/**
 * @brief computes the BWT positions of T$ in chunks of the dictionary suffix array. The positions whose suffixes start
 *        with a dictionary suffix of length > w are processed in the order of the dictionary suffix array:
 *        1. a single full word s: the BWT symbols are the symbols before the phrases equal to s, in the order of the
 *           ranks of the parse suffixes after them (psi[C[s+1]..C[s+2]-1]);
 *        2. a single proper suffix of a word: all its positions have the symbol before the suffix in the word;
 *        3. a group of equal proper suffixes of several words: the ranks of the parse suffixes after the words are
 *           merged with a heap (one run if all the words have the same symbol before the suffix and no full arrays
 *           are written).
 *        The phrases of a word with a phrase run of two or more phrases (a repeated word) are ordered by the
 *        suffixes of P (the parse of the phrases) starting at them (see repeated_word_ranks); for a group with such
 *        a word, the ranks of the suffixes of P after the phrases are merged.
 */
template <typename sym_t, typename idx_t>
class bwt_merger {
public:
    // ############################# TYPES AND CONSTANTS #############################

    using range_t = rank_range<idx_t>;

private:
    // ranks argument of append_positions for positions whose ranks are not consecutive in psi
    static constexpr uint64_t no_ranks = UINT64_MAX;

public:
    // ############################# DATA MEMBERS #############################

    const uint64_t n, w, J; // n = |T|, w and J (see parsing)
    const dictionary<sym_t>& dict;
    const dictionary_sa<idx_t>& dsa;
    const parse_sa<idx_t>& psa;
    const phrase_arrays<idx_t, sym_t>& arr;
    const params& prm;
    const uint8_t position_bytes; // bytes per value of the full SA and LCP array and per position of a record
    const record_layout layout; // the layout of the records of the runs (see section_records)
    const bool writes_run_data; // true <=> the runs file is written (see ext_runs)
    // true <=> the SA values at the run boundaries (the LCP values at the run starts) are written
    const bool writes_sa_samples, writes_lcp_samples;
    const bool collects_run_starts; // true <=> the values at the starts of the runs are collected
    const sampled_plcp* plcp; // the sampled PLCP array, if the full LCP array is written
    const bool has_repeated; // true <=> a word is repeated (see psa.repeated)

private:
    // the files of the full BWT, SA and LCP array (nullptr if not written)
    std::unique_ptr<parallel_output_file> bwt_file, sa_file, lcp_file;

    // ############################# INTERNAL METHODS #############################

    /* appends cnt BWT positions with symbol sym, of which first is the first and last the last, to out; sa_of(i)
       returns the SA value of position i (in [0..cnt-1]), which is only called with the full SA or LCP array */
    template <typename sa_of_t>
    inline void append_positions_with(chunk_runs<sym_t>& out, sym_t sym, uint64_t cnt, const bwt_pos& first,
        const bwt_pos& last, sa_of_t sa_of) const;

    /* appends cnt BWT positions with symbol sym, of which first is the first and last the last, to out; the ranks of
       their parse suffixes are psi[ranks..ranks+cnt-1] (no_ranks if they are not consecutive in psi, which is only
       allowed without the full SA and LCP array) */
    inline void append_positions(chunk_runs<sym_t>& out, sym_t sym, uint64_t cnt, const bwt_pos& first,
        const bwt_pos& last, uint64_t ranks) const
    {
        append_positions_with(out, sym, cnt, first, last,
            [&](uint64_t i) { return sa(first.dict_suffix_len, psa.psi[ranks + i]); });
    }

    /* returns the BWT position with the dictionary suffix dict.text[x..x+len-1] whose phrase starts the suffix of P
       with offset q among those starting with word s (rw = psa.find_repeated(s)) */
    inline bwt_pos position(uint64_t x, uint64_t len, uint64_t s, const repeated_word_ranks* rw, uint64_t q) const
    {
        if (rw == nullptr) return { x, len, psa.psi[psa.C[s + 1] + q], 0 };
        const auto [run, j] = rw->run_of(q);
        return { x, len, psa.psi[rw->first_run + run], (j - 1) * rw->advance };
    }

    // appends the BWT positions of the full word s at x of length len, a repeated word, to out
    [[gnu::noinline]] void append_repeated_word(chunk_runs<sym_t>& out, uint64_t x, uint64_t len, uint64_t s) const;

    // appends the BWT positions of the proper suffix dict.text[x..x+len-1] of the repeated word s to out
    [[gnu::noinline]] void append_repeated_suffix(chunk_runs<sym_t>& out, uint64_t x, uint64_t len, uint64_t s) const;

    /* appends the BWT positions of the group of equal dictionary suffixes group[0..num-1] (a range of the dictionary
       suffix array) of length len, of which one or more belong to a repeated word, to out */
    [[gnu::noinline]] void append_repeated_group(chunk_runs<sym_t>& out, const idx_t* group, uint64_t num,
        uint64_t len, std::vector<repeated_range>& heap) const;

    /* returns the LCP value of two BWT positions whose common prefix has length l up to the parse suffixes a and b
       (see run_walk) */
    [[gnu::noinline]] uint64_t lcp_in_runs(uint64_t l, run_walk a, run_walk b) const;

    // writes the buffered part of the full BWT, SA and LCP array of out
    void write_full(chunk_runs<sym_t>& out) const;

    // appends the record of the last run of out, whose SA value at its end is sa_end, to out.records
    void append_record(chunk_runs<sym_t>& out, uint64_t sa_end) const;

public:
    // ############################# CONSTRUCTORS #############################

    bwt_merger(const parsing<sym_t>& prs, const dictionary_sa<idx_t>& dsa, const parse_sa<idx_t>& psa,
        const phrase_arrays<idx_t, sym_t>& arr, const params& prm, const std::string& prefix,
        bool writes_run_data, bool collects_run_starts, const sampled_plcp* plcp)
        : n(prs.n), w(prs.w), J(prs.num_phrase_runs), dict(prs.dict), dsa(dsa), psa(psa), arr(arr), prm(prm),
          position_bytes(byte_width(prs.n + 1)),
          layout(position_bytes, sizeof(sym_t), prm.write_sa_samples, prm.write_lcp_samples),
          writes_run_data(writes_run_data),
          writes_sa_samples(writes_run_data && prm.write_sa_samples),
          writes_lcp_samples(writes_run_data && prm.write_lcp_samples), collects_run_starts(collects_run_starts),
          plcp(plcp), has_repeated(!psa.repeated.empty())
    {
        if (prm.write_bwt) bwt_file = std::make_unique<parallel_output_file>(prefix + ext_bwt);
        if (prm.write_sa) sa_file = std::make_unique<parallel_output_file>(prefix + ext_sa);
        if (plcp != nullptr) lcp_file = std::make_unique<parallel_output_file>(prefix + ext_lcp);
    }

    // writes the rest of the full BWT, SA and LCP array and closes their files (see parallel_output_file::close)
    void close_full_files()
    {
        for (std::unique_ptr<parallel_output_file>* f : { &bwt_file, &sa_file, &lcp_file }) {
            if (*f != nullptr) (*f)->close();
        }
    }

    // ############################# ACCESS METHODS #############################

    /* returns the SA value (in T$) of the position with dictionary suffix length len, parse suffix rank j and back
       (see bwt_pos) */
    inline uint64_t sa(uint64_t len, uint64_t j, uint64_t back = 0) const
    {
        const uint64_t x = arr.ends_by_rank[j] + 1 - len - back; // position in S
        return x == 0 ? n : x - 1;
    }

    inline uint64_t sa(const bwt_pos& b) const { return sa(b.dict_suffix_len, b.parse_rank, b.back); }

    // returns the BWT symbol of the full-word position whose parse suffix has rank j
    inline sym_t symbol_before_full_word(uint64_t j) const { return bwt_symbol<sym_t>(arr.symbols_by_rank[j]); }

    // returns the BWT symbol before the phrase run starting the parse suffix of rank x >= 1
    inline sym_t symbol_before_run_of_rank(uint64_t x) const { return symbol_before_full_word(psa.psi[x]); }

    // returns the length of the common prefix of the words s1 != s2 after their first w symbols
    inline uint64_t words_lcp(uint64_t s1, uint64_t s2) const
    {
        const sym_t* u = dict.text.get() + dict.word_start[s1] + w;
        const sym_t* v = dict.text.get() + dict.word_start[s2] + w;
        uint64_t q = 0;
        while (u[q] == v[q]) q++;
        return q;
    }

    // returns the LCP value of two BWT positions
    uint64_t lcp(const bwt_pos& a, const bwt_pos& b) const;

    // ############################# CONSTRUCTION METHODS #############################

    // returns the number of BWT positions of the dictionary suffixes chunk_sa[0..num-1] (a range of SA)
    uint64_t count_bwt_positions(const idx_t* chunk_sa, uint64_t num) const;

    /**
     * @brief computes the BWT positions of the dictionary suffixes SA[beg..end-1]
     * @param chunk_sa SA[beg..end-1]
     * @param beg first rank in the dictionary suffix array
     * @param end rank after the last one
     * @param start BWT position of the first of them
     * @param out the chunk (it is cleared first)
     * @param heaps memory for the heaps of the calling thread
     */
    void compute_chunk(const idx_t* chunk_sa, uint64_t beg, uint64_t end, uint64_t start, chunk_runs<sym_t>& out,
        merge_heaps<idx_t>& heaps) const;
};

/**
 * @brief writes the runs file (see ext_runs): appends the chunks in their order and joins a run that continues from one
 *        chunk into the next. The records are written to the runs file (plain format) or to a temporary file, from
 *        which finish computes the compressed format (see write_compressed_runs).
 */
template <typename sym_t, typename idx_t>
class output_writer {
private:
    // ############################# DATA MEMBERS #############################

    const bwt_merger<sym_t, idx_t>& m;
    const record_layout& layout; // the layout of a record
    const bool compressed; // true <=> the runs file is written in the compressed format
    // where the values at the starts of the runs are appended (nullptr: they are not collected)
    std::vector<run_start_values>* run_starts;
    direct_ofstream file; // the file of the records (not opened if the runs file is not written)
    std::string file_name; // name of the runs file
    std::string records_name; // name of the file of the records
    uint64_t pos = 0; // BWT position of the first position of the next chunk
    uint64_t num_runs = 0; // number of written runs
    uint64_t max_length = 0; // the largest run length
    uint64_t max_symbol = 0; // the largest run symbol
    std::vector<uint64_t> samples; // plain: the starts of the runs 0, runs_per_block, 2 * runs_per_block, ...
    // compressed: the number of runs of each symbol, for the symbols that fit into two bytes and for the others
    std::vector<uint64_t> dense_counts;
    std::unordered_map<uint64_t, uint64_t> sparse_counts;
    bool pending = false; // true <=> a run has been read but not written (it can continue in the next chunk)
    sym_t pending_sym = 0;
    uint64_t pending_len = 0;
    uint64_t pending_start = 0; // the start of the pending run
    uint64_t pending_sa_start = 0, pending_lcp = 0; // the SA value and the LCP value at its start
    bwt_pos last_read_pos {}; // the last BWT position read

    // ############################# INTERNAL METHODS #############################

    // counts the run with the record rec (samples, symbol counts, largest symbol)
    void count_run(const uint8_t* rec);

    // writes the pending run, whose SA value at its end is sa_end
    void write_pending(uint64_t sa_end);

    // returns the fields of the runs file
    run_data_fields fields() const;

public:
    // ############################# CONSTRUCTORS #############################

    output_writer(const bwt_merger<sym_t, idx_t>& m, const std::string& prefix,
        std::vector<run_start_values>* run_starts);

    ~output_writer();

    output_writer(const output_writer&) = delete;
    output_writer& operator=(const output_writer&) = delete;

    // ############################# METHODS #############################

    // appends a chunk
    void append_chunk(chunk_runs<sym_t>& chunk);

    // writes the last run, the other sections, the extra sections and the table of contents; returns the number of runs
    uint64_t finish();
};

/**
 * @brief computes the BWT of T$ in the order of the dictionary suffixes and writes it run-length encoded to the runs
 *        file, optionally together with the SA values at the starts and ends of the runs and the LCP values at the
 *        starts of the runs, and optionally with one symbol per position, the suffix array and the LCP array
 * @param prs the parsing (n, w, J and the dictionary)
 * @param dsa the generalized suffix array of the dictionary
 * @param psa the suffix array of the parsing
 * @param arr the phrase runs (ends_by_rank and symbols_by_rank)
 * @param prm parameters
 * @param prefix path prefix of the output files
 * @param writes_run_data true <=> the runs file is written (see ext_runs)
 * @param run_starts if not nullptr, the SA and LCP values at the starts of the runs are appended to it
 * @param plcp if not nullptr, the full LCP array is written, computed with this sampled PLCP array
 * @return the number r of runs
 */
template <typename sym_t, typename idx_t>
uint64_t merge(const parsing<sym_t>& prs, const dictionary_sa<idx_t>& dsa, const parse_sa<idx_t>& psa,
    const phrase_arrays<idx_t, sym_t>& arr, const params& prm, const std::string& prefix, bool writes_run_data,
    std::vector<run_start_values>* run_starts = nullptr, const sampled_plcp* plcp = nullptr);

} // namespace pfp_rlbwt

#include <pfp_rlbwt/construction/merge.tpp>

#endif // PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_MERGE_HPP
