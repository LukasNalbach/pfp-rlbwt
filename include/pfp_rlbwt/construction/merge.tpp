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

#ifndef PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_MERGE_TPP
#define PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_MERGE_TPP

#include <pfp_rlbwt/construction/merge.hpp>

namespace pfp_rlbwt {

// ############################# bwt_merger #############################

template <typename sym_t, typename idx_t>
uint64_t bwt_merger<sym_t, idx_t>::lcp(const bwt_pos& a, const bwt_pos& b) const
{
    const sym_t* t = dict.text.get();
    const uint64_t m = std::min(a.dict_suffix_len, b.dict_suffix_len);
    uint64_t l = 0;
    while (l < m && t[a.dict_pos + l] == t[b.dict_pos + l]) l++;
    if (l < m || a.dict_suffix_len != b.dict_suffix_len) return l;

    // the same dictionary suffix: the text continues with the phrases of the parse suffixes (overlapping it by w)
    uint64_t j1 = a.parse_rank, j2 = b.parse_rank;

    if ((a.back | b.back) != 0) [[unlikely]] {
        return lcp_in_runs(l, { a.back == 0 ? 0 : dsa.word_of(a.dict_pos), a.back, j1 },
            { b.back == 0 ? 0 : dsa.word_of(b.dict_pos), b.back, j2 });
    }

    while (j1 != j2 && j1 != 0 && j2 != 0) {
        const uint64_t s1 = psa.word_of_rank(j1);
        const uint64_t s2 = psa.word_of_rank(j2);

        if (s1 != s2) return l + words_lcp(s1, s2);

        if (has_repeated && psa.repeated_set.contains(s1)) [[unlikely]] {
            return lcp_in_runs(l, { 0, 0, j1 }, { 0, 0, j2 });
        }

        l += dict.length(s1) - w;
        j1 = psa.psi[j1];
        j2 = psa.psi[j2];
    }

    return l;
}

template <typename sym_t, typename idx_t>
uint64_t bwt_merger<sym_t, idx_t>::lcp_in_runs(uint64_t l, run_walk a, run_walk b) const
{
    // continues v with the phrase run of rank v.next, whose word is s
    auto load = [&](run_walk& v, uint64_t s) {
        v.word = s;

        if (psa.repeated_set.contains(s)) {
            const uint64_t i = psa.SA[v.next];
            v.rem = arr.ends[i] - (i == 0 ? w - 1 : arr.ends[i - 1]);
        } else {
            v.rem = dict.length(s) - w;
        }

        v.next = psa.psi[v.next];
    };

    while (true) {
        if (a.rem == 0 && b.rem == 0) {
            if (a.next == 0 || b.next == 0) return l;
            const uint64_t s = psa.word_of_rank(a.next);
            // the suffix of rank b.next starts with s iff b.next lies in the range of the suffixes starting with s
            if (b.next < psa.C[s + 1] || b.next >= psa.C[s + 2]) return l + words_lcp(s, psa.word_of_rank(b.next));
            load(a, s);
            load(b, s);
        } else {
            for (run_walk* v : { &a, &b }) {
                if (v->rem != 0) continue;
                if (v->next == 0) return l;
                load(*v, psa.word_of_rank(v->next));
            }

            if (a.word != b.word) return l + words_lcp(a.word, b.word);
        }

        const uint64_t c = std::min(a.rem, b.rem);
        l += c;
        a.rem -= c;
        b.rem -= c;
    }
}

template <typename sym_t, typename idx_t>
uint64_t bwt_merger<sym_t, idx_t>::count_bwt_positions(uint64_t beg, uint64_t end) const
{
    uint64_t cnt = 0;

    for (uint64_t i = beg; i < end; i++) {
        const uint64_t x = dsa.SA[i];
        const uint64_t s = dsa.word_of(x);
        if (dict.word_start[s + 1] - 1 - x > w) cnt += dict.num_occ[s];
    }

    return cnt;
}

template <typename sym_t, typename idx_t>
void bwt_merger<sym_t, idx_t>::write_full(chunk_runs<sym_t>& out) const
{
    if (prm.write_bwt && !out.full_bwt.empty()) {
        bwt_file->write(out.full_bwt.data(), out.full_bwt.size() * sizeof(sym_t), out.full_pos * sizeof(sym_t));
    }

    if (prm.write_sa && !out.full_sa.empty()) {
        sa_file->write(out.full_sa.data(), out.full_sa.size(), out.full_pos * position_bytes);
    }

    if (plcp != nullptr && !out.full_lcp.empty()) {
        lcp_file->write(out.full_lcp.data(), out.full_lcp.size(), out.full_pos * position_bytes);
    }

    out.full_pos += prm.write_bwt ? out.full_bwt.size() : std::max(out.full_sa.size(), out.full_lcp.size()) / position_bytes;
    out.full_bwt.clear();
    out.full_sa.clear();
    out.full_lcp.clear();
}

template <typename sym_t, typename idx_t>
void bwt_merger<sym_t, idx_t>::append_record(chunk_runs<sym_t>& out, uint64_t sa_end) const
{
    const uint64_t offset = out.records.size();
    resize_no_init(out.records, offset + layout.bytes);
    uint8_t* rec = out.records.data() + offset;
    write_bytes(rec, out.last_start, layout.position_bytes);
    write_bytes(rec + layout.symbol, uint64_t(out.last_sym), layout.symbol_bytes);

    if (writes_sa_samples) {
        write_bytes(rec + layout.sa_start, out.last_sa_start, layout.position_bytes);
        write_bytes(rec + layout.sa_end, sa_end, layout.position_bytes);
    }

    if (writes_lcp_samples) write_bytes(rec + layout.lcp, out.last_lcp, layout.position_bytes);
    out.max_inner_length = std::max(out.max_inner_length, out.last_len);
}

template <typename sym_t, typename idx_t>
template <typename sa_of_t>
inline void bwt_merger<sym_t, idx_t>::append_positions_with(chunk_runs<sym_t>& out, sym_t sym, uint64_t cnt,
    const bwt_pos& first, const bwt_pos& last, sa_of_t sa_of) const
{
    if (out.length == 0) {
        out.first_pos = first;
        out.num_runs = 1;
        out.last_sym = sym;
        out.last_len = cnt;
    } else if (sym == out.last_sym) {
        out.last_len += cnt;
    } else {
        // the run (last_sym, last_len) ends at last_pos
        const uint64_t sa_last = writes_sa_samples ? sa(out.last_pos) : 0;

        if (out.num_runs == 1) {
            out.first_sym = out.last_sym;
            out.first_len = out.last_len;
            out.first_sa_end = sa_last;
        } else if (writes_run_data) {
            append_record(out, sa_last);
        }

        const uint64_t sa_first = writes_sa_samples || collects_run_starts ? sa(first) : 0;
        uint64_t l = 0;

        if (writes_lcp_samples || collects_run_starts) {
            l = lcp(out.last_pos, first);
            if (collects_run_starts) out.run_starts.emplace_back(sa_first, l);
        }

        out.last_start = out.length;
        out.last_sa_start = sa_first;
        out.last_lcp = l;
        out.num_runs++;
        out.last_sym = sym;
        out.last_len = cnt;
    }

    out.last_pos = last;
    out.length += cnt;

    if (prm.write_bwt) out.full_bwt.insert(out.full_bwt.end(), cnt, sym);

    if (prm.write_sa || plcp != nullptr) {
        for (uint64_t i = 0; i < cnt; i++) {
            const uint64_t v = sa_of(i);

            if (prm.write_sa) {
                out.full_sa.insert(out.full_sa.end(), (const uint8_t*) &v, (const uint8_t*) &v + position_bytes);
            }

            if (plcp != nullptr) {
                const uint64_t l = plcp->get(v);
                out.full_lcp.insert(out.full_lcp.end(), (const uint8_t*) &l, (const uint8_t*) &l + position_bytes);
            }
        }
    }

    if (out.full_bwt.size() * sizeof(sym_t) >= write_buffer_bytes || out.full_sa.size() >= write_buffer_bytes
        || out.full_lcp.size() >= write_buffer_bytes) {
        write_full(out);
    }
}

template <typename sym_t, typename idx_t>
void bwt_merger<sym_t, idx_t>::append_repeated_word(chunk_runs<sym_t>& out, uint64_t x, uint64_t len, uint64_t s) const
{
    const repeated_word_ranks& rw = *psa.find_repeated(s);
    const uint64_t occ = rw.occ;
    const uint64_t first_run = psa.C[s + 1];
    auto pos_of = [&](uint64_t q) { return position(x, len, s, &rw, q); };
    // a phrase that is not the first one of its run is preceded by the symbol inner, the first one by its own symbol
    const sym_t inner = bwt_symbol<sym_t>(arr.symbol_before_next[psa.SA[first_run]]);
    uint64_t run_beg = 0; // the positions run_beg, ... have the symbol sym
    sym_t sym = inner;

    // appends the positions run_beg..q-1 if the positions q, ... have the symbol c != sym
    auto add = [&](uint64_t q, sym_t c) {
        if (q == 0) sym = c;
        if (c == sym) return;
        const uint64_t b = run_beg;
        append_positions_with(out, sym, q - b, pos_of(b), pos_of(q - 1), [&](uint64_t i) { return sa(pos_of(b + i)); });
        run_beg = q;
        sym = c;
    };

    uint64_t q = 0;

    for (uint64_t h = 0; h < rw.num_runs; h++) {
        const uint64_t offset = rw.head_offset[h];
        if (offset > q) add(q, inner);
        add(offset, symbol_before_run(psa.SA[first_run + h]));
        q = offset + 1;
    }

    if (occ > q) add(q, inner);
    const uint64_t b = run_beg;
    append_positions_with(out, sym, occ - b, pos_of(b), pos_of(occ - 1), [&](uint64_t i) { return sa(pos_of(b + i)); });
}

template <typename sym_t, typename idx_t>
void bwt_merger<sym_t, idx_t>::append_repeated_suffix(chunk_runs<sym_t>& out, uint64_t x, uint64_t len,
    uint64_t s) const
{
    const repeated_word_ranks* rw = psa.find_repeated(s);
    auto pos_of = [&](uint64_t q) { return position(x, len, s, rw, q); };
    append_positions_with(out, bwt_symbol<sym_t>(dict.text[x - 1]), rw->occ, pos_of(0), pos_of(rw->occ - 1),
        [&](uint64_t i) { return sa(pos_of(i)); });
}

template <typename sym_t, typename idx_t>
void bwt_merger<sym_t, idx_t>::append_repeated_group(chunk_runs<sym_t>& out, uint64_t i, uint64_t j, uint64_t len,
    std::vector<repeated_range>& heap) const
{
    const sym_t* t = dict.text.get();
    auto sym_of = [&](uint64_t member) { return bwt_symbol<sym_t>(t[dsa.SA[i + member] - 1]); };
    // orders the ranges in heap by the rank of the suffix one phrase after their next one, the smallest first
    auto later = [](const repeated_range& a, const repeated_range& b) { return a.next_rank > b.next_rank; };
    // the dictionary suffix of a position is the one of its own word (see lcp)
    auto pos_of = [&](const repeated_range& r, uint64_t q) {
        return position(dsa.SA[i + r.member], len, r.word, r.rw, q);
    };

    // appends the positions of r with the offsets [r.cur..e-1]
    auto append = [&](const repeated_range& r, uint64_t e) {
        append_positions_with(out, sym_of(r.member), e - r.cur, pos_of(r, r.cur), pos_of(r, e - 1),
            [&](uint64_t k) { return sa(pos_of(r, r.cur + k)); });
    };

    heap.clear();

    for (uint64_t q = i; q < j; q++) {
        const uint64_t sq = dsa.word_of(dsa.SA[q]);
        const repeated_word_ranks* rwq = psa.find_repeated(sq);
        heap.push_back({ 0, dict.num_occ[sq], psa.next_rank_in_P(sq, rwq, 0), sq, rwq, q - i });
    }

    const sym_t sym = sym_of(0);
    bool same_sym = true;
    for (uint64_t q = i + 1; q < j && same_sym; q++) same_sym = sym_of(q - i) == sym;

    if (same_sym && !prm.write_sa && plcp == nullptr) {
        // a single run, from the smallest first to the largest last rank
        const repeated_range* first = &heap.front();
        const repeated_range* last = &heap.front();
        uint64_t last_rank = 0, cnt = 0;

        for (const repeated_range& r : heap) {
            const uint64_t r_last = psa.next_rank_in_P(r.word, r.rw, r.end - 1);
            if (r.next_rank < first->next_rank) first = &r;

            if (r_last >= last_rank) {
                last_rank = r_last;
                last = &r;
            }

            cnt += r.end;
        }

        append_positions(out, sym, cnt, pos_of(*first, 0), pos_of(*last, last->end - 1), no_ranks);
        return;
    }

    std::make_heap(heap.begin(), heap.end(), later);

    while (heap.size() > 1) {
        // the offsets of the first range with next ranks smaller than the next rank of every other range
        repeated_range& top = heap.front();
        uint64_t bound = heap[1].next_rank;
        if (heap.size() > 2) bound = std::min(bound, heap[2].next_rank);
        uint64_t lo = top.cur + 1, hi = top.end;

        while (lo < hi) {
            const uint64_t mid = (lo + hi) / 2;
            if (psa.next_rank_in_P(top.word, top.rw, mid) < bound) lo = mid + 1;
            else hi = mid;
        }

        append(top, lo);
        top.cur = lo;

        if (top.cur == top.end) {
            std::pop_heap(heap.begin(), heap.end(), later);
            heap.pop_back();
        } else {
            top.next_rank = psa.next_rank_in_P(top.word, top.rw, top.cur);
            // move the first range down to its place
            uint64_t q = 0;

            while (true) {
                uint64_t c = 2 * q + 1;
                if (c >= heap.size()) break;
                if (c + 1 < heap.size() && later(heap[c], heap[c + 1])) c++;
                if (!later(heap[q], heap[c])) break;
                std::swap(heap[q], heap[c]);
                q = c;
            }
        }
    }

    append(heap.front(), heap.front().end);

    if (heap.capacity() > kept_heap_entries) {
        heap.clear();
        heap.shrink_to_fit();
    }
}

template <typename sym_t, typename idx_t>
void bwt_merger<sym_t, idx_t>::compute_chunk(uint64_t beg, uint64_t end, uint64_t start, chunk_runs<sym_t>& out,
    merge_heaps<idx_t>& heaps) const
{
    std::vector<range_t>& heap = heaps.plain;
    out.clear(start);
    const sym_t* t = dict.text.get();
    const index_array<idx_t>& psi = psa.psi;
    // orders the ranges in heap by their next rank, the smallest first
    auto later = [&](const range_t& a, const range_t& b) { return psi[a.cur] > psi[b.cur]; };
    uint64_t i = beg;

    while (i < end) {
        const uint64_t x = dsa.SA[i];
        const uint64_t s = dsa.word_of(x);
        const uint64_t len = dict.word_start[s + 1] - 1 - x;

        if (len <= w) {
            i++;
            continue;
        }

        // the dictionary suffixes SA[i..j-1] are equal
        uint64_t j = i + 1;
        while (j < end && dsa.equals_previous(j)) j++;
        const uint64_t list = psa.C[s + 1]; // the ranks of the parse suffixes preceded by word s are psi[list..]
        const uint64_t occ = psa.C[s + 2] - list;
        auto pos_of = [&](uint64_t r) { return bwt_pos { x, len, psi[r], 0 }; };

        if (has_repeated && psa.repeated_set.contains(s)) [[unlikely]] {
            if (j == i + 1 && x == dict.word_start[s]) append_repeated_word(out, x, len, s);
            else if (j == i + 1) append_repeated_suffix(out, x, len, s);
            else append_repeated_group(out, i, j, len, heaps.repeated);
            i = j;
            continue;
        }

        if (j == i + 1 && x == dict.word_start[s]) {
            // the full word s: the symbol preceding each of its phrases
            uint64_t run_beg = 0;
            sym_t sym = symbol_before_full_word(psi[list]);

            for (uint64_t q = 1; q < occ; q++) {
                const sym_t c = symbol_before_full_word(psi[list + q]);
                if (c == sym) continue;
                append_positions(out, sym, q - run_beg, pos_of(list + run_beg), pos_of(list + q - 1), list + run_beg);
                run_beg = q;
                sym = c;
            }

            append_positions(out, sym, occ - run_beg, pos_of(list + run_beg), pos_of(list + occ - 1), list + run_beg);
        } else if (j == i + 1) {
            append_positions(out, bwt_symbol<sym_t>(t[x - 1]), occ, pos_of(list), pos_of(list + occ - 1), list);
        } else {
            // the words containing the suffix as a proper suffix; their ranges of ranks are merged
            auto sym_of = [&](uint64_t member) { return bwt_symbol<sym_t>(t[dsa.SA[i + member] - 1]); };
            bool repeated = false; // true <=> a word of the group is repeated
            // the ranks of the parse suffixes preceded by the word of SA[q] are psi[first..last-1]
            auto range_of = [&](uint64_t q) {
                const uint64_t sq = dsa.word_of(dsa.SA[q]);
                if (has_repeated) repeated |= psa.repeated_set.contains(sq);
                return std::make_pair(psa.C[sq + 1], psa.C[sq + 2]);
            };

            const sym_t sym = sym_of(0);
            bool same_sym = true;
            for (uint64_t q = i + 1; q < j && same_sym; q++) same_sym = sym_of(q - i) == sym;

            if (same_sym && !prm.write_sa && plcp == nullptr) {
                // a single run, from the smallest first to the largest last rank
                uint64_t first = UINT64_MAX, last = 0, cnt = 0;

                for (uint64_t q = i; q < j; q++) {
                    const auto [b, e] = range_of(q);
                    first = std::min<uint64_t>(first, psi[b]);
                    last = std::max<uint64_t>(last, psi[e - 1]);
                    cnt += e - b;
                }

                if (repeated) [[unlikely]] append_repeated_group(out, i, j, len, heaps.repeated);
                else append_positions(out, sym, cnt, { x, len, first, 0 }, { x, len, last, 0 }, no_ranks);
            } else {
                heap.clear();

                for (uint64_t q = i; q < j; q++) {
                    const auto [b, e] = range_of(q);
                    heap.push_back({ idx_t(b), idx_t(e), uint32_t(q - i) });
                }

                if (repeated) [[unlikely]] {
                    append_repeated_group(out, i, j, len, heaps.repeated);
                    i = j;
                    continue;
                }

                std::make_heap(heap.begin(), heap.end(), later);

                while (heap.size() > 1) {
                    // the ranks of the first range that are smaller than the next rank of every other range
                    range_t& top = heap.front();
                    uint64_t bound = psi[heap[1].cur];
                    if (heap.size() > 2) bound = std::min(bound, psi[heap[2].cur]);
                    uint64_t lo = top.cur + 1, hi = top.end;

                    while (lo < hi) {
                        const uint64_t mid = (lo + hi) / 2;
                        if (psi[mid] < bound) lo = mid + 1;
                        else hi = mid;
                    }

                    append_positions(out, sym_of(top.member), lo - top.cur, pos_of(top.cur), pos_of(lo - 1), top.cur);
                    top.cur = idx_t(lo);

                    if (top.cur == top.end) {
                        std::pop_heap(heap.begin(), heap.end(), later);
                        heap.pop_back();
                    } else {
                        // move the first range down to its place
                        uint64_t q = 0;

                        while (true) {
                            uint64_t c = 2 * q + 1;
                            if (c >= heap.size()) break;
                            if (c + 1 < heap.size() && later(heap[c], heap[c + 1])) c++;
                            if (!later(heap[q], heap[c])) break;
                            std::swap(heap[q], heap[c]);
                            q = c;
                        }
                    }
                }

                const range_t& l = heap.front();
                append_positions(out, sym_of(l.member), l.end - l.cur, pos_of(l.cur), pos_of(l.end - 1), l.cur);

                if (heap.capacity() > kept_heap_entries) {
                    heap.clear();
                    heap.shrink_to_fit();
                }
            }
        }

        i = j;
    }

    if (out.num_runs == 1) {
        out.first_sym = out.last_sym;
        out.first_len = out.last_len;
    }

    write_full(out);
}

// ############################# output_writer #############################

template <typename sym_t, typename idx_t>
output_writer<sym_t, idx_t>::output_writer(const bwt_merger<sym_t, idx_t>& m, const std::string& prefix,
    std::vector<run_start_values>* run_starts)
    : m(m), layout(m.layout), compressed(m.prm.format == file_format::compressed), run_starts(run_starts)
{
    if (!m.writes_run_data) return;
    file_name = prefix + ext_runs;
    records_name = compressed ? file_name + ext_temporary_records : file_name;
    file.open(records_name);
    if (!file.good()) throw std::runtime_error("pfp_rlbwt: cannot write " + records_name);

    if (compressed && sizeof(sym_t) <= sizeof(uint16_t)) {
        dense_counts.assign(uint64_t { 1 } << (byte_bits * sizeof(sym_t)), 0);
    }
}

template <typename sym_t, typename idx_t>
output_writer<sym_t, idx_t>::~output_writer()
{
    if (!compressed || records_name.empty()) return;
    std::error_code ec;
    if (file.is_open()) file.close();
    std::filesystem::remove(records_name, ec);
}

template <typename sym_t, typename idx_t>
inline void output_writer<sym_t, idx_t>::count_run(const uint8_t* rec)
{
    const uint64_t sym = read_bytes(rec + layout.symbol, layout.symbol_bytes);
    max_symbol = std::max(max_symbol, sym);

    if (compressed) {
        if (sym < dense_counts.size()) dense_counts[sym]++;
        else sparse_counts[sym]++;
    } else if (num_runs % m.prm.runs_per_block == 0) {
        samples.push_back(read_bytes(rec, layout.position_bytes));
    }

    num_runs++;
}

template <typename sym_t, typename idx_t>
void output_writer<sym_t, idx_t>::write_pending(uint64_t sa_end)
{
    max_length = std::max(max_length, pending_len);

    if (!m.writes_run_data) {
        num_runs++;
        return;
    }

    uint8_t rec[sizeof(uint64_t) * 5];
    write_bytes(rec, pending_start, layout.position_bytes);
    write_bytes(rec + layout.symbol, uint64_t(pending_sym), layout.symbol_bytes);

    if (m.writes_sa_samples) {
        write_bytes(rec + layout.sa_start, pending_sa_start, layout.position_bytes);
        write_bytes(rec + layout.sa_end, sa_end, layout.position_bytes);
    }

    if (m.writes_lcp_samples) write_bytes(rec + layout.lcp, pending_lcp, layout.position_bytes);
    count_run(rec);
    file.write((const char*) rec, layout.bytes);
}

template <typename sym_t, typename idx_t>
void output_writer<sym_t, idx_t>::append_chunk(chunk_runs<sym_t>& chunk)
{
    if (chunk.length == 0) return;

    if (pending && pending_sym == chunk.first_sym) {
        pending_len += chunk.first_len;
    } else {
        // a run starts at the first position of the chunk
        if (pending) write_pending(m.writes_sa_samples ? m.sa(last_read_pos) : 0);
        const uint64_t sa_first = m.writes_sa_samples || run_starts != nullptr ? m.sa(chunk.first_pos) : 0;
        uint64_t l = 0;

        if (m.writes_lcp_samples || run_starts != nullptr) {
            l = pending ? m.lcp(last_read_pos, chunk.first_pos) : 0;
            if (run_starts != nullptr) run_starts->emplace_back(sa_first, l);
        }

        pending = true;
        pending_sym = chunk.first_sym;
        pending_len = chunk.first_len;
        pending_start = pos;
        pending_sa_start = sa_first;
        pending_lcp = l;
    }

    // the first run of the chunk ends inside the chunk: write it and the inner runs; the last run stays pending
    if (chunk.num_runs > 1) {
        write_pending(chunk.first_sa_end);
        max_length = std::max(max_length, chunk.max_inner_length);

        if (m.writes_run_data) {
            // the starts of the records become BWT positions
            for (uint64_t offset = 0; offset < chunk.records.size(); offset += layout.bytes) {
                uint8_t* rec = chunk.records.data() + offset;
                write_bytes(rec, pos + read_bytes(rec, layout.position_bytes), layout.position_bytes);
                count_run(rec);
            }

            file.write((const char*) chunk.records.data(), chunk.records.size());
        } else {
            num_runs += chunk.num_runs - 2;
        }

        if (run_starts != nullptr) {
            run_starts->insert(run_starts->end(), chunk.run_starts.begin(), chunk.run_starts.end());
        }

        pending_sym = chunk.last_sym;
        pending_len = chunk.last_len;
        pending_start = pos + chunk.last_start;
        pending_sa_start = chunk.last_sa_start;
        pending_lcp = chunk.last_lcp;
    }

    pos += chunk.length;
    last_read_pos = chunk.last_pos;
}

template <typename sym_t, typename idx_t>
run_data_fields output_writer<sym_t, idx_t>::fields() const
{
    run_data_fields f;
    f.n = m.n + 1;
    f.r = num_runs;
    f.max_symbol = max_symbol;
    f.max_run_length = max_length;
    f.sa_samples = m.writes_sa_samples;
    f.lcp_samples = m.writes_lcp_samples;
    f.runs_per_block = m.prm.runs_per_block;
    f.runs_per_symbol_sample = compressed ? m.prm.runs_per_symbol_sample : 0;
    f.position_bytes = layout.position_bytes;
    f.symbol_bytes = layout.symbol_bytes;
    f.position_bits = std::bit_width(f.n);
    return f;
}

template <typename sym_t, typename idx_t>
uint64_t output_writer<sym_t, idx_t>::finish()
{
    write_pending(m.writes_sa_samples ? m.sa(last_read_pos) : 0);
    if (!m.writes_run_data) return num_runs;
    run_data_fields f = fields();

    if (compressed) {
        file.close();
        if (file.fail()) throw std::runtime_error("pfp_rlbwt: cannot write " + records_name);
        std::vector<std::pair<uint64_t, uint64_t>> frequencies;
        for (uint64_t c = 0; c < dense_counts.size(); c++) if (dense_counts[c] != 0) frequencies.emplace_back(c, dense_counts[c]);
        frequencies.insert(frequencies.end(), sparse_counts.begin(), sparse_counts.end());
        std::sort(frequencies.begin(), frequencies.end());
        write_compressed_runs(records_name, file_name, f, std::move(frequencies), m.prm);
        return num_runs;
    }

    // the records, the record samples, the fields, the extra sections and the table of contents (see ext_runs)
    std::vector<run_data_entry> contents;
    contents.push_back({ .id = section_records, .offset = 0, .size = file.num_bytes_written() });
    uint64_t offset = file.num_bytes_written();
    for (uint64_t s : samples) file.write((const char*) &s, layout.position_bytes);
    contents.push_back({ .id = section_record_samples, .offset = offset, .size = file.num_bytes_written() - offset });
    offset = file.num_bytes_written();
    file.write((const char*) &f, sizeof(f));
    contents.push_back({ .id = section_fields, .offset = offset, .size = sizeof(f) });

    for (const extra_section& extra : m.prm.extra_sections) {
        offset = file.num_bytes_written();
        extra.write(file);
        contents.push_back({ .id = extra.id, .offset = offset, .size = file.num_bytes_written() - offset });
    }

    const uint64_t num_sections = contents.size();
    file.write((const char*) contents.data(), contents.size() * sizeof(run_data_entry));
    file.write((const char*) &num_sections, sizeof(num_sections));
    file.close();
    if (file.fail()) throw std::runtime_error("pfp_rlbwt: cannot write " + file_name);
    return num_runs;
}

// ############################# merge #############################

template <typename sym_t, typename idx_t>
uint64_t merge(const parsing<sym_t>& prs, const dictionary_sa<idx_t>& dsa, const parse_sa<idx_t>& psa,
    const phrase_arrays<idx_t, sym_t>& arr, const params& prm, const std::string& prefix, bool writes_run_data,
    std::vector<run_start_values>* run_starts, const sampled_plcp* plcp)
{
    const uint16_t p = prm.num_threads;
    const uint64_t size = prs.dict.size;
    bwt_merger<sym_t, idx_t> m(prs, dsa, psa, arr, prm, prefix, writes_run_data, run_starts != nullptr, plcp);
    output_writer<sym_t, idx_t> writer(m, prefix, run_starts);

    /* the chunks of the dictionary suffix array: [0..num_chunks] chunk c covers the ranks [bounds[c], bounds[c+1]);
       a group of equal suffixes is not split */
    const uint64_t num_chunks = std::max<uint64_t>(1, std::min<uint64_t>(size, p * chunks_per_thread));
    std::vector<uint64_t> bounds(num_chunks + 1);

    for (uint64_t c = 0; c <= num_chunks; c++) {
        uint64_t b = (c * size) / num_chunks;
        while (b < size && dsa.equals_previous(b)) b++;
        bounds[c] = b;
    }

    // [0..num_chunks] chunk_starts[c] = BWT position of the first position of chunk c (only for the full arrays)
    std::vector<uint64_t> chunk_starts(num_chunks + 1, 0);

    if (prm.write_bwt || prm.write_sa || plcp != nullptr) {
        #pragma omp parallel for num_threads(p) schedule(dynamic, 1)
        for (uint64_t c = 0; c < num_chunks; c++) chunk_starts[c + 1] = m.count_bwt_positions(bounds[c], bounds[c + 1]);

        for (uint64_t c = 0; c < num_chunks; c++) chunk_starts[c + 1] += chunk_starts[c];
    }

    /* a thread computes the next chunk; the thread completing the next chunk to write writes all completed ones.
       slots[c % num_slots] stores chunk c until it is written; done[c % num_slots] = 1 <=> chunk c is computed and not
       yet written */
    const uint64_t num_slots = p * pending_chunks_per_thread;
    std::vector<chunk_runs<sym_t>> slots(num_slots);
    std::vector<char> done(num_slots, 0);
    std::atomic<uint64_t> next_chunk { 0 };
    uint64_t written = 0; // number of written chunks
    bool writing = false; // true <=> a thread is writing chunks
    std::mutex lock;
    std::condition_variable slot_free;
    std::exception_ptr error;

    #pragma omp parallel num_threads(p)
    {
        merge_heaps<idx_t> heaps;

        while (true) {
            const uint64_t c = next_chunk.fetch_add(1);
            if (c >= num_chunks) break;

            // wait until the slot of chunk c is free
            {
                std::unique_lock<std::mutex> guard(lock);
                slot_free.wait(guard, [&] { return c < written + num_slots || error; });
                if (error) break;
            }

            chunk_runs<sym_t>& out = slots[c % num_slots];

            try {
                m.compute_chunk(bounds[c], bounds[c + 1], chunk_starts[c], out, heaps);
            } catch (...) {
                std::lock_guard<std::mutex> guard(lock);
                error = std::current_exception();
                slot_free.notify_all();
                break;
            }

            std::unique_lock<std::mutex> guard(lock);
            done[c % num_slots] = 1;
            if (writing) continue;
            writing = true;

            // write the completed chunks in their order
            while (written < num_chunks && done[written % num_slots] && !error) {
                chunk_runs<sym_t>& o = slots[written % num_slots];
                guard.unlock();

                try {
                    writer.append_chunk(o);
                } catch (...) {
                    guard.lock();
                    error = std::current_exception();
                    break;
                }

                guard.lock();
                done[written % num_slots] = 0;
                written++;
                slot_free.notify_all();
            }

            writing = false;
            slot_free.notify_all();
        }
    }

    if (error) std::rethrow_exception(error);
    const uint64_t r = writer.finish();
    m.close_full_files();
    return r;
}

} // namespace pfp_rlbwt

#endif // PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_MERGE_TPP
