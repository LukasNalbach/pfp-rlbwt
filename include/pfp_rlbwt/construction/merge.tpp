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
inline bool bwt_merger<sym_t, idx_t>::start_lcp(const bwt_pos& a, const bwt_pos& b, lcp_walk& wk) const
{
    const sym_t* t = dict.text.get();
    const uint64_t m = std::min(a.dict_suffix_len, b.dict_suffix_len);
    uint64_t l = 0;
    while (l < m && t[a.dict_pos + l] == t[b.dict_pos + l]) l++;
    wk.l = l;
    if (l < m || a.dict_suffix_len != b.dict_suffix_len) return false;

    /* the same dictionary suffix: the suffixes continue with the back symbols of the phrase run of the dictionary
       suffix and then with the phrase runs of the parse suffixes (overlapping them by w), whose words and advances the
       rows of Psi store (see psi_mds) */
    wk.u = { a.back == 0 ? 0 : dsa.word_of(a.dict_pos), a.back, a.parse_rank, 0 };
    wk.v = { b.back == 0 ? 0 : dsa.word_of(b.dict_pos), b.back, b.parse_rank, 0 };
    return true;
}

template <typename sym_t, typename idx_t>
inline bool bwt_merger<sym_t, idx_t>::step_lcp(lcp_walk& wk) const
{
    const psi_mds& psi = psa.psi;
    run_walk& u = wk.u;
    run_walk& v = wk.v;

    // continues r with the phrase run of rank r.next >= 1 and fetches the row that its next move reads
    auto load = [&](run_walk& r) {
        r.word = psi.word(r.x);
        r.rem = psi.advance(r.x);
        psi.move(r.next, r.x);
        psi.prefetch_move(r.x);
    };

    if (u.rem == 0 && v.rem == 0) {
        if (u.next == v.next || u.next == 0 || v.next == 0) return true;
        load(u);
        load(v);
    } else {
        for (run_walk* r : { &u, &v }) {
            if (r->rem != 0) continue;
            if (r->next == 0) return true;
            load(*r);
        }
    }

    if (u.word != v.word) {
        wk.l += words_lcp(u.word, v.word);
        return true;
    }

    const uint64_t c = std::min(u.rem, v.rem);
    wk.l += c;
    u.rem -= c;
    v.rem -= c;
    return false;
}

template <typename sym_t, typename idx_t>
uint64_t bwt_merger<sym_t, idx_t>::lcp(const bwt_pos& a, const bwt_pos& b) const
{
    lcp_walk wk;

    if (start_lcp(a, b, wk)) {
        wk.u.x = psa.psi.interval_of(wk.u.next);
        wk.v.x = psa.psi.interval_of(wk.v.next);
        while (!step_lcp(wk)) { }
    }

    return wk.l;
}

template <typename sym_t, typename idx_t>
void bwt_merger<sym_t, idx_t>::compute_lcps(chunk_runs<sym_t>& out) const
{
    const psi_mds& psi = psa.psi;
    std::vector<uint64_t>& values = out.lcps->values;
    std::vector<lcp_walk>& walks = out.lcps->walks;
    const uint64_t num_walks = walks.size();

    // starts walk i and fetches the cache lines of the input intervals of the walk interleaved_lcp_walks later
    auto start = [&](uint64_t i) {
        if (i + interleaved_lcp_walks < num_walks) {
            psi.prefetch_interval_of(walks[i + interleaved_lcp_walks].u.next);
            psi.prefetch_interval_of(walks[i + interleaved_lcp_walks].v.next);
        }

        for (run_walk* r : { &walks[i].u, &walks[i].v }) {
            r->x = psi.interval_of(r->next);
            psi.prefetch_move(r->x);
        }
    };

    // [0..active-1] lane[t] = the walk that lane t advances
    std::array<uint64_t, interleaved_lcp_walks> lane;
    uint64_t active = 0, next = 0;

    for (; active < interleaved_lcp_walks && next < num_walks; next++) {
        start(next);
        lane[active++] = next;
    }

    while (active > 0) {
        for (uint64_t t = 0; t < active;) {
            lcp_walk& wk = walks[lane[t]];

            if (!step_lcp(wk)) {
                t++;
                continue;
            }

            values[wk.index - out.lcps->first] = wk.l;

            if (next < num_walks) {
                start(next);
                lane[t++] = next++;
            } else {
                lane[t] = lane[--active];
            }
        }
    }

    walks.clear();
    // value j of the chunk is the LCP value at the start of run j + 2, whose record is record j (if it is not the last)
    const uint64_t first = out.lcps->first, num = values.size();
    if (num == 0) return;
    const uint64_t records = writes_run_data ? out.records.size() / layout.bytes : 0;

    if (writes_lcp_samples) {
        for (uint64_t j = first; j < std::min(first + num, records); j++) {
            write_bytes(out.records.data() + j * layout.bytes + layout.lcp, values[j - first], layout.position_bytes);
        }
    }

    if (collects_run_starts) {
        for (uint64_t j = 0; j < num; j++) out.run_starts[first + j].second = values[j];
    }

    out.last_lcp = values[num - 1];
    out.lcps->first = first + num - 1;
    values.front() = values[num - 1];
    values.resize(1);
}

template <typename sym_t, typename idx_t>
uint64_t bwt_merger<sym_t, idx_t>::count_bwt_positions(const idx_t* chunk_sa, uint64_t num) const
{
    uint64_t cnt = 0;

    for (uint64_t i = 0; i < num; i++) {
        const uint64_t x = chunk_sa[i];
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

        // the LCP value at the start of the run is computed with those of the other runs (see compute_lcps)
        if (writes_lcp_samples || collects_run_starts) {
            lcp_walk wk;

            if (start_lcp(out.last_pos, first, wk)) {
                wk.index = out.lcps->first + out.lcps->values.size();
                out.lcps->walks.push_back(wk);
            }

            out.lcps->values.push_back(wk.l);
            if (collects_run_starts) out.run_starts.emplace_back(sa_first, 0);
            if (out.lcps->values.size() == pending_lcp_values) compute_lcps(out);
        }

        out.last_start = out.length;
        out.last_sa_start = sa_first;
        out.last_lcp = 0;
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
    /* a phrase that is not the first one of its run is preceded by the symbol inner, the first one by its own symbol;
       inner precedes each phrase after a phrase of s, such as the first phrase of the run after the run of rank
       first_run */
    const sym_t inner = symbol_at_rank(psa.psi[first_run]);
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
        add(offset, symbol_at_rank(first_run + h));
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
void bwt_merger<sym_t, idx_t>::append_repeated_group(chunk_runs<sym_t>& out, const idx_t* group, uint64_t num,
    uint64_t len, std::vector<repeated_range>& heap) const
{
    const sym_t* t = dict.text.get();
    auto sym_of = [&](uint64_t member) { return bwt_symbol<sym_t>(t[uint64_t(group[member]) - 1]); };
    // orders the ranges in heap by the rank of the suffix one phrase after their next one, the smallest first
    auto later = [](const repeated_range& a, const repeated_range& b) { return a.next_rank > b.next_rank; };
    // the dictionary suffix of a position is the one of its own word (see lcp)
    auto pos_of = [&](const repeated_range& r, uint64_t q) {
        return position(group[r.member], len, r.word, r.rw, q);
    };

    // appends the positions of r with the offsets [r.cur..e-1]
    auto append = [&](const repeated_range& r, uint64_t e) {
        append_positions_with(out, sym_of(r.member), e - r.cur, pos_of(r, r.cur), pos_of(r, e - 1),
            [&](uint64_t k) { return sa(pos_of(r, r.cur + k)); });
    };

    heap.clear();

    for (uint64_t q = 0; q < num; q++) {
        const uint64_t sq = dsa.word_of(group[q]);
        const repeated_word_ranks* rwq = psa.find_repeated(sq);
        heap.push_back({ 0, dict.num_occ[sq], psa.next_rank_among_phrases(sq, rwq, 0), sq, rwq, q });
    }

    const sym_t sym = sym_of(0);
    bool same_sym = true;
    for (uint64_t q = 1; q < num && same_sym; q++) same_sym = sym_of(q) == sym;

    if (same_sym && !prm.write_sa && plcp == nullptr) {
        // a single run, from the smallest first to the largest last rank
        const repeated_range* first = &heap.front();
        const repeated_range* last = &heap.front();
        uint64_t last_rank = 0, cnt = 0;

        for (const repeated_range& r : heap) {
            const uint64_t r_last = psa.next_rank_among_phrases(r.word, r.rw, r.end - 1);
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
            if (psa.next_rank_among_phrases(top.word, top.rw, mid) < bound) lo = mid + 1;
            else hi = mid;
        }

        append(top, lo);
        top.cur = lo;

        if (top.cur == top.end) {
            std::pop_heap(heap.begin(), heap.end(), later);
            heap.pop_back();
        } else {
            top.next_rank = psa.next_rank_among_phrases(top.word, top.rw, top.cur);
            sift_down_first(heap, later);
        }
    }

    append(heap.front(), heap.front().end);

    if (heap.capacity() > kept_heap_entries) {
        heap.clear();
        heap.shrink_to_fit();
    }
}

template <typename sym_t, typename idx_t>
void bwt_merger<sym_t, idx_t>::compute_chunk(const idx_t* chunk_sa, uint64_t beg, uint64_t end, uint64_t start,
    chunk_runs<sym_t>& out, merge_heaps<idx_t>& heaps) const
{
    std::vector<range_t>& heap = heaps.plain;
    out.clear(start);
    out.lcps = &heaps.lcps;
    out.lcps->first = 0;
    out.lcps->values.clear();
    const sym_t* t = dict.text.get();
    const psi_mds& psi = psa.psi;
    // orders the ranges in heap by their next rank, the smallest first
    auto later = [](const range_t& a, const range_t& b) { return a.next_rank > b.next_rank; };
    // returns SA[k] for k in [beg, end)
    auto SA = [&](uint64_t k) { return uint64_t(chunk_sa[k - beg]); };
    /* the dictionary suffixes before ahead_words have their words fetched, those before ahead_data the data of their
       words, and those before ahead_rows the rows of Psi of the first and the last rank of their words */
    uint64_t ahead_words = beg, ahead_data = beg, ahead_rows = beg;
    uint64_t i = beg;

    while (i < end) {
        for (; ahead_words < std::min(end, i + prefetched_suffixes); ahead_words++) {
            dsa.prefetch_word_of(SA(ahead_words));
        }

        for (; ahead_data < std::min(end, i + prefetched_suffixes / 2); ahead_data++) {
            const uint64_t y = SA(ahead_data), u = dsa.word_of(y);
            dict.word_start.prefetch(u + 1);
            psa.C.prefetch(u + 1);
            psa.C_intervals.prefetch(u + 1);
            __builtin_prefetch(t + y - 1);
        }

        for (; ahead_rows < std::min(end, i + prefetched_suffixes / 4); ahead_rows++) {
            const uint64_t u = dsa.word_of(SA(ahead_rows));
            psi.prefetch_row(psa.C_intervals[u + 1]);
            psi.prefetch_row(psa.C_intervals[u + 2] - 1);
        }

        const uint64_t x = SA(i);
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
        // returns the BWT position whose parse suffix follows the phrase run of rank r (in the input interval z of Psi)
        auto pos_in = [&](uint64_t r, uint64_t z) { return bwt_pos { x, len, psi.at(r, z), 0, r }; };

        if (has_repeated && psa.repeated_set.contains(s)) [[unlikely]] {
            if (j == i + 1 && x == dict.word_start[s]) append_repeated_word(out, x, len, s);
            else if (j == i + 1) append_repeated_suffix(out, x, len, s);
            else append_repeated_group(out, chunk_sa + (i - beg), j - i, len, heaps.repeated);
            i = j;
            continue;
        }

        if (j == i + 1 && x == dict.word_start[s]) {
            /* the full word s: the symbol preceding each of its phrases, which is the same on a run of the BWT of the
               parse of the phrase runs (see symbol_at_rank) */
            const uint64_t last = list + occ;
            uint64_t k = psa.bwt_runs.run_of(list); // the run of the BWT of the parse containing the current rank
            uint64_t run_beg = list;
            sym_t sym = run_symbol(k);
            uint64_t z = psa.C_intervals[s + 1]; // the input interval of Psi containing the rank of the next position

            // returns the BWT position of rank r (r does not decrease from call to call)
            auto next_pos = [&](uint64_t r) {
                while (psi.p(z + 1) <= r) z++;
                return pos_in(r, z);
            };

            // appends the positions of the ranks [run_beg, e) with the symbol sym
            auto append_run = [&](uint64_t e) {
                const bwt_pos first_pos = next_pos(run_beg);
                append_positions(out, sym, e - run_beg, first_pos, next_pos(e - 1), run_beg);
            };

            for (uint64_t y = psa.bwt_runs.next_start(list); y < last; y = psa.bwt_runs.next_start(y)) {
                const sym_t c = run_symbol(++k);
                if (c == sym) continue;
                append_run(y);
                run_beg = y;
                sym = c;
            }

            append_run(last);
        } else if (j == i + 1) {
            append_positions(out, bwt_symbol<sym_t>(t[x - 1]), occ, pos_in(list, psa.C_intervals[s + 1]),
                pos_in(list + occ - 1, psa.C_intervals[s + 2] - 1), list);
        } else {
            // the words containing the suffix as a proper suffix; their ranges of ranks are merged
            auto sym_of = [&](uint64_t member) { return bwt_symbol<sym_t>(t[SA(i + member) - 1]); };
            bool repeated = false; // true <=> a word of the group is repeated
            // returns the word of SA[q], the ranks of the parse suffixes preceded by it are psi[C[sq+1]..C[sq+2]-1]
            auto word_of = [&](uint64_t q) {
                const uint64_t sq = dsa.word_of(SA(q));
                if (has_repeated) repeated |= psa.repeated_set.contains(sq);
                return sq;
            };

            const sym_t sym = sym_of(0);
            bool same_sym = true;
            for (uint64_t q = i + 1; q < j && same_sym; q++) same_sym = sym_of(q - i) == sym;

            if (same_sym && !prm.write_sa && plcp == nullptr) {
                // a single run, from the smallest first to the largest last parse rank
                uint64_t first = UINT64_MAX, last = 0, cnt = 0;
                uint64_t first_rank = 0, last_rank = 0; // the ranks with psi[first_rank] = first, psi[last_rank] = last

                for (uint64_t q = i; q < j; q++) {
                    const uint64_t sq = word_of(q);
                    const uint64_t b = psa.C[sq + 1], e = psa.C[sq + 2];
                    const uint64_t psi_b = psi.q(psa.C_intervals[sq + 1]);
                    const uint64_t psi_e = psi.at(e - 1, psa.C_intervals[sq + 2] - 1);

                    if (psi_b < first) {
                        first = psi_b;
                        first_rank = b;
                    }

                    if (psi_e >= last) {
                        last = psi_e;
                        last_rank = e - 1;
                    }

                    cnt += e - b;
                }

                if (repeated) [[unlikely]] {
                    append_repeated_group(out, chunk_sa + (i - beg), j - i, len, heaps.repeated);
                } else {
                    append_positions(out, sym, cnt, { x, len, first, 0, first_rank }, { x, len, last, 0, last_rank },
                        no_ranks);
                }
            } else {
                heap.clear();

                for (uint64_t q = i; q < j; q++) {
                    const uint64_t sq = word_of(q);
                    const uint64_t x_b = psa.C_intervals[sq + 1];
                    heap.push_back({ idx_t(psi.q(x_b)), idx_t(psa.C[sq + 1]), idx_t(psa.C[sq + 2]), idx_t(x_b),
                        idx_t(psa.C_intervals[sq + 2]), uint32_t(q - i) });
                }

                if (repeated) [[unlikely]] {
                    append_repeated_group(out, chunk_sa + (i - beg), j - i, len, heaps.repeated);
                    i = j;
                    continue;
                }

                /* returns the first input interval z > r.x of Psi with z = r.end_x or q(z) >= bound: Psi maps the
                   ranks of the intervals r.x, ..., z-1 before bound (whole intervals: Psi maps an interval into a run
                   of the BWT of the parse, and bound is a position of this BWT with another value); galloping search */
                auto skip = [&](const range_t& r, uint64_t bound) {
                    auto before = [&](uint64_t z) { return z < r.end_x && psi.q(z) < bound; };
                    uint64_t z = r.x + 1, hi = z;

                    for (uint64_t step = 1; hi < r.end_x && before(hi); step *= 2) {
                        z = hi + 1;
                        hi = std::min<uint64_t>(hi + step, r.end_x);
                    }

                    while (z < hi) {
                        const uint64_t mid = (z + hi) / 2;
                        if (before(mid)) z = mid + 1;
                        else hi = mid;
                    }

                    return z;
                };

                if (!prm.write_sa && plcp == nullptr) {
                    /* a run per stretch of ranks with the same symbol before the suffix: the ranges with the same
                       symbol form a class with a heap, the classes form a heap by their next ranks, and the class with
                       the smallest next rank takes all its ranks before the smallest next rank of the other classes */
                    std::vector<symbol_class>& classes = heaps.classes;
                    classes.clear();

                    std::sort(heap.begin(), heap.end(),
                        [&](const range_t& a, const range_t& b) { return sym_of(a.member) < sym_of(b.member); });

                    for (uint64_t b = 0; b < heap.size();) {
                        const sym_t c = sym_of(heap[b].member);
                        uint64_t e = b + 1;
                        while (e < heap.size() && sym_of(heap[e].member) == c) e++;
                        std::make_heap(heap.begin() + b, heap.begin() + e, later);
                        classes.push_back({ b, e - b, uint64_t(c) });
                        b = e;
                    }

                    // orders the classes by their next ranks, the smallest first
                    auto later_class = [&](const symbol_class& a, const symbol_class& b) {
                        return heap[a.beg].next_rank > heap[b.beg].next_rank;
                    };

                    std::make_heap(classes.begin(), classes.end(), later_class);

                    while (!classes.empty()) {
                        symbol_class& cl = classes.front();
                        uint64_t bound = std::numeric_limits<uint64_t>::max();
                        if (classes.size() > 1) bound = heap[classes[1].beg].next_rank;
                        if (classes.size() > 2) bound = std::min<uint64_t>(bound, heap[classes[2].beg].next_rank);
                        const auto first = heap.begin() + cl.beg;
                        const bwt_pos first_pos { x, len, first->next_rank, 0, first->cur };
                        // the number of ranks taken and the largest next rank taken, last = psi[last_rank]
                        uint64_t cnt = 0, last = 0, last_rank = 0;

                        while (cl.size > 0 && first->next_rank < bound) {
                            std::pop_heap(first, first + cl.size, later);
                            range_t& r = first[cl.size - 1];
                            const uint64_t z = skip(r, bound), e = psi.p(z);
                            cnt += e - r.cur;
                            const uint64_t psi_last = psi.at(e - 1, z - 1);

                            if (psi_last >= last) {
                                last = psi_last;
                                last_rank = e - 1;
                            }

                            if (e == r.end) {
                                cl.size--;
                            } else {
                                r.cur = idx_t(e);
                                r.x = idx_t(z);
                                r.next_rank = idx_t(psi.q(z));
                                std::push_heap(first, first + cl.size, later);
                            }
                        }

                        append_positions(out, sym_t(cl.symbol), cnt, first_pos, bwt_pos { x, len, last, 0, last_rank },
                            no_ranks);

                        if (cl.size == 0) {
                            std::pop_heap(classes.begin(), classes.end(), later_class);
                            classes.pop_back();
                        } else {
                            sift_down_first(classes, later_class);
                        }
                    }
                } else {
                    std::make_heap(heap.begin(), heap.end(), later);

                    while (heap.size() > 1) {
                        // the ranks of the first range that are smaller than the next rank of every other range
                        range_t& top = heap.front();
                        uint64_t bound = heap[1].next_rank;
                        if (heap.size() > 2) bound = std::min<uint64_t>(bound, heap[2].next_rank);
                        const uint64_t z = skip(top, bound), e = psi.p(z);
                        append_positions(out, sym_of(top.member), e - top.cur,
                            bwt_pos { x, len, top.next_rank, 0, top.cur }, pos_in(e - 1, z - 1), top.cur);
                        top.cur = idx_t(e);
                        top.x = idx_t(z);

                        if (top.cur == top.end) {
                            std::pop_heap(heap.begin(), heap.end(), later);
                            heap.pop_back();
                        } else {
                            top.next_rank = idx_t(psi.q(z));
                            sift_down_first(heap, later);
                        }
                    }

                    const range_t& l = heap.front();
                    append_positions(out, sym_of(l.member), l.end - l.cur, bwt_pos { x, len, l.next_rank, 0, l.cur },
                        pos_in(l.end - 1, l.end_x - 1), l.cur);
                }

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

    compute_lcps(out);
    out.lcps = nullptr;
    write_full(out);
}

// ############################# output_writer #############################

template <typename sym_t, typename idx_t>
output_writer<sym_t, idx_t>::output_writer(const bwt_merger<sym_t, idx_t>& m, const std::string& prefix,
    std::vector<run_start_values>* run_starts, level_output* level_out)
    : m(m), layout(m.layout), compressed(m.prm.format == file_format::compressed), run_starts(run_starts)
{
    if (!m.writes_run_data) return;

    if (level_out != nullptr) {
        memory_records = &level_out->records;
        level_out->layout = layout;
        return;
    }

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
    write_records(rec, layout.bytes);
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

            write_records(chunk.records.data(), chunk.records.size());
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
    f.lcp_values_per_block = compressed ? m.prm.lcp_values_per_block : 0;
    f.position_bytes = layout.position_bytes;
    f.symbol_bytes = layout.symbol_bytes;
    f.position_bits = std::bit_width(f.n);
    return f;
}

template <typename sym_t, typename idx_t>
uint64_t output_writer<sym_t, idx_t>::finish()
{
    write_pending(m.writes_sa_samples ? m.sa(last_read_pos) : 0);
    if (!m.writes_run_data || memory_records != nullptr) return num_runs;
    run_data_fields f = fields();

    if (compressed) {
        file.close();
        if (file.fail()) throw std::runtime_error("pfp_rlbwt: cannot write " + records_name);
        std::vector<std::pair<uint64_t, uint64_t>> frequencies;
        for (uint64_t c = 0; c < dense_counts.size(); c++) if (dense_counts[c] != 0) frequencies.emplace_back(c, dense_counts[c]);
        frequencies.insert(frequencies.end(), sparse_counts.begin(), sparse_counts.end());
        std::sort(frequencies.begin(), frequencies.end());
        write_compressed_runs(records_name, file_name, f, std::move(frequencies), m.prm, m.prm.num_threads);
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
    std::vector<run_start_values>* run_starts, const sampled_plcp* plcp, level_output* level_out)
{
    const uint16_t p = prm.num_threads;
    const uint64_t size = prs.dict.size;
    bwt_merger<sym_t, idx_t> m(prs, dsa, psa, arr, prm, prefix, writes_run_data, run_starts != nullptr, plcp,
        level_out);
    output_writer<sym_t, idx_t> writer(m, prefix, run_starts, level_out);

    /* the chunks of the dictionary suffix array: [0..num_chunks] chunk c covers the ranks [bounds[c], bounds[c+1]);
       a group of equal suffixes is not split */
    const uint64_t num_chunks = std::max<uint64_t>(1, std::min<uint64_t>(size, p * chunks_per_thread));
    std::vector<uint64_t> bounds(num_chunks + 1);

    for (uint64_t c = 0; c <= num_chunks; c++) {
        uint64_t b = (c * size) / num_chunks;
        while (b < size && dsa.equals_previous(b)) b++;
        bounds[c] = b;
    }

    // a reader of the dictionary suffix array reads the largest chunk, but at most temporary_read_bytes bytes, at once
    uint64_t max_chunk = 1;
    for (uint64_t c = 0; c < num_chunks; c++) max_chunk = std::max(max_chunk, bounds[c + 1] - bounds[c]);
    const uint64_t read_bytes = std::min<uint64_t>(max_chunk * sizeof(idx_t), temporary_read_bytes);

    // [0..num_chunks] chunk_starts[c] = BWT position of the first position of chunk c (only for the full arrays)
    std::vector<uint64_t> chunk_starts(num_chunks + 1, 0);

    if (prm.write_bwt || prm.write_sa || plcp != nullptr) {
        std::vector<std::exception_ptr> failures(p);

        #pragma omp parallel num_threads(p)
        {
            try {
                array_reader sa_file(dsa.sa_file, sizeof(idx_t), read_bytes);
                no_init_vector<idx_t> chunk_sa;

                #pragma omp for schedule(dynamic, 1)
                for (uint64_t c = 0; c < num_chunks; c++) {
                    dsa.read_sa(sa_file, bounds[c], bounds[c + 1], chunk_sa);
                    chunk_starts[c + 1] = m.count_bwt_positions(chunk_sa.data(), chunk_sa.size());
                }
            } catch (...) {
                failures[omp_get_thread_num()] = std::current_exception();
            }
        }

        rethrow_first_exception(failures);
        for (uint64_t c = 0; c < num_chunks; c++) chunk_starts[c + 1] += chunk_starts[c];
    }

    /* a thread computes the next chunk; the thread completing the next chunk to write writes all completed ones.
       slots[c % num_slots] stores chunk c until it is written (then its memory is freed); done[c % num_slots] = 1 <=>
       chunk c is computed and not yet written */
    const uint64_t num_slots = p * pending_chunks_per_thread;
    std::vector<chunk_runs<sym_t>> slots(num_slots);
    std::vector<char> done(num_slots, 0);
    std::atomic<uint64_t> next_chunk { 0 };
    uint64_t written = 0; // number of written chunks
    bool writing = false; // true <=> a thread is writing chunks
    std::mutex lock;
    std::condition_variable slot_free;
    std::exception_ptr error;
    // counts the BWT positions of the written chunks, if chunk_starts is computed, and else their ranks
    const bool counts_positions = chunk_starts[num_chunks] != 0;
    phase_progress prog(prm.log, counts_positions ? chunk_starts[num_chunks] : size);

    #pragma omp parallel num_threads(p)
    {
        merge_heaps<idx_t> heaps;
        std::unique_ptr<array_reader> sa_file;
        no_init_vector<idx_t> chunk_sa; // the values of the dictionary suffix array of the chunk

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
                if (sa_file == nullptr) {
                    sa_file = std::make_unique<array_reader>(dsa.sa_file, sizeof(idx_t), read_bytes);
                }

                dsa.read_sa(*sa_file, bounds[c], bounds[c + 1], chunk_sa);
                m.compute_chunk(chunk_sa.data(), bounds[c], bounds[c + 1], chunk_starts[c], out, heaps);
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
                    o.release();
                    prog.reached(counts_positions ? chunk_starts[written + 1] : bounds[written + 1]);
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
    prog.finish();
    const uint64_t r = writer.finish();
    m.close_full_files();
    return r;
}

} // namespace pfp_rlbwt

#endif // PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_MERGE_TPP
