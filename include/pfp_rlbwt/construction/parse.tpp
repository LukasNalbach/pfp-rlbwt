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
 * @file parse.tpp
 * @brief the methods of the prefix-free parsing: parsing ranges, the trigger seed and counting and renaming symbols
 */

#ifndef PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_PARSE_TPP
#define PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_PARSE_TPP

#include <pfp_rlbwt/construction/parse.hpp>

namespace pfp_rlbwt {

template <typename reader_t>
inline uint64_t read_varint(reader_t& in)
{
    uint64_t x = 0;

    for (uint8_t shift = 0;; shift += 7) {
        const uint8_t c = in.next();
        x |= uint64_t(c & 127) << shift;
        if (c < 128) return x;
    }
}

inline uint32_t word_id_vector::reader::next()
{
    if (i == v.blocks[b].size) {
        b++;
        i = 0;
    }

    return get(v.blocks[b], i++);
}

inline void word_id_vector::add_block(uint32_t id)
{
    const uint8_t width = id > short_mask || (!blocks.empty() && blocks.back().width == long_bytes)
        ? long_bytes : short_bytes;
    capacity = blocks.empty() ? block_vector<uint32_t>::min_block_size
                              : std::min<uint64_t>(2 * capacity, max_block_bytes / width);
    blocks.push_back({ new_malloc_array<uint8_t>(capacity * width + long_bytes - width), 0, width });
}

inline void word_id_vector::push_back(uint32_t id)
{
    if (blocks.empty() || blocks.back().size == capacity || (id > short_mask && blocks.back().width == short_bytes))
        [[unlikely]] {
        add_block(id);
    }

    block& b = blocks.back();
    std::memcpy(b.bytes.get() + b.size * b.width, &id, b.width);
    b.size++;
}

template <typename fnc_t>
inline void word_id_vector::move_to(fnc_t fnc, uint16_t threads)
{
    // offsets[b] = number of ids in the blocks before block b
    std::vector<uint64_t> offsets(blocks.size() + 1, 0);
    for (uint64_t b = 0; b < blocks.size(); b++) offsets[b + 1] = offsets[b] + blocks[b].size;
    std::vector<std::exception_ptr> failures(blocks.size());

    #pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
    for (uint64_t b = 0; b < blocks.size(); b++) {
        try {
            for (uint64_t i = 0; i < blocks[b].size; i++) fnc(offsets[b] + i, get(blocks[b], i));
        } catch (...) {
            failures[b] = std::current_exception();
        }

        blocks[b].bytes.reset();
    }

    blocks.clear();
    capacity = 0;
    rethrow_first_exception(failures);
}

inline uint64_t string_hash(const uint8_t* s, uint64_t len)
{
    uint64_t h = len;
    uint64_t i = 0;

    for (; i + sizeof(uint64_t) <= len; i += sizeof(uint64_t)) {
        uint64_t v;
        std::memcpy(&v, s + i, sizeof(uint64_t));
        h = mix(h ^ v);
    }

    uint64_t v = 0;
    std::memcpy(&v, s + i, len - i);
    return mix(h ^ v);
}

inline const void* word_arena::store(const void* s, uint64_t len)
{
    if (len > num_free) {
        block_size = blocks.empty() ? min_arena_block_size : std::min(2 * block_size, max_arena_block_size);
        const uint64_t size = std::max(len, block_size);
        blocks.emplace_back(new_malloc_array<uint8_t>(size));
        free_begin = blocks.back().get();
        num_free = size;
    }

    std::memcpy(free_begin, s, len);
    const void* copy = free_begin;
    free_begin += len;
    num_free -= len;
    return copy;
}

template <typename sym_t>
inline auto word_table<sym_t>::block_of(uint64_t id) -> entry_t*
{
    std::atomic<entry_t*>& b = blocks[id >> log_entries_per_block];
    entry_t* e = b.load(std::memory_order_acquire);
    if (e != nullptr) return e;
    std::lock_guard<std::mutex> guard(blocks_lock);
    e = b.load(std::memory_order_relaxed);

    if (e == nullptr) {
        e = new entry_t[entry_mask + 1];
        b.store(e, std::memory_order_release);
    }

    return e;
}

template <typename sym_t>
inline void word_table<sym_t>::grow(shard& sh)
{
    std::vector<slot> old(2 * sh.slots.size(), slot { 0, 0 });
    old.swap(sh.slots);
    const uint64_t mask = sh.slots.size() - 1;

    for (const slot& s : old) {
        if (s.id_plus_1 == 0) continue;
        uint64_t i = s.hash & mask;
        while (sh.slots[i].id_plus_1 != 0) i = (i + 1) & mask;
        sh.slots[i] = s;
    }
}

template <typename sym_t>
inline uint32_t word_table<sym_t>::count_occurrences(const sym_t* s, uint64_t len, uint64_t hash, uint64_t count, word_arena& arena)
{
    shard& sh = shards[shard_shift == 64 ? 0 : hash >> shard_shift];
    std::lock_guard<std::mutex> guard(sh.lock);
    const uint64_t mask = sh.slots.size() - 1;

    for (uint64_t i = hash & mask;; i = (i + 1) & mask) {
        slot& sl = sh.slots[i];

        if (sl.id_plus_1 == 0) {
            const uint64_t id = next_id.fetch_add(1, std::memory_order_relaxed);
            if (id >= max_words) throw std::runtime_error("pfp_rlbwt: too many distinct phrases");
            block_of(id)[id & entry_mask] =
                entry_t { (const sym_t*) arena.store(s, len * sizeof(sym_t)), len, count };
            sl = slot { hash, uint32_t(id + 1) };
            if (2 * ++sh.used > sh.slots.size()) grow(sh);
            return id;
        }

        if (sl.hash == hash) {
            entry_t& e = entry(sl.id_plus_1 - 1);

            if (e.len == len && std::memcmp(e.str, s, len * sizeof(sym_t)) == 0) {
                e.num_occ += count;
                return sl.id_plus_1 - 1;
            }
        }
    }
}

template <typename sym_t>
inline uint64_t parse_phrases_in_range(typename text_ranges<sym_t>::reader& in, bool first, const params& prm,
    uint64_t seed, word_table<sym_t>& table, word_arena& arena, thread_phrases& out, uint64_t block_len)
{
    const uint64_t w = prm.window;
    rolling_window_hash<sym_t> window(w, prm.modulus, seed);
    std::vector<sym_t> phrase; // S[b_k..x], where x is the position read last
    bool in_phrase = first;
    if (first) phrase.push_back(dollar);
    uint64_t x = 1; // next position to read
    const uint64_t check_from = x + w - 1; // the first position at which the window has w symbols
    // the position after the range (unknown_length until the reader has passed the end of the range)
    uint64_t range_end = text_ranges<sym_t>::unknown_length;
    uint64_t max_sym = 0;
    /* the w symbols before the current block (0 before the range): w_before[i] is removed from the window when
       block[i] is appended to it, for i < w (for i >= w, block[i - w] is removed) */
    std::vector<sym_t> w_before(w, 0);

    std::vector<sym_t> run_phrase; // the phrase of the last phrase run, if it is not yet appended to out
    uint64_t m = 0; // m, the number of phrases of that run (0: there is none)

    // a word of a phrase run and its occurrences that are not yet added to the word table
    struct counted_word {
        uint64_t hash = 0;
        uint64_t id = word_table<sym_t>::max_words; // max_words: none
        uint64_t pending = 0;
    };

    /* the words of the last two phrase runs appended to out, the last one first (consecutive phrase runs have
       different words, so only the older one can be the word of the next run) */
    std::array<counted_word, 2> recent;

    // appends the last phrase run to out
    auto append_run = [&]() {
        if (m == 0) return;
        const uint64_t len = run_phrase.size();
        const uint64_t hash = string_hash((const uint8_t*) run_phrase.data(), len * sizeof(sym_t));
        counted_word& older = recent[1];

        if (older.hash == hash && older.id != word_table<sym_t>::max_words
                && table.is_word(older.id, run_phrase.data(), len)) {
            older.pending += m;
            std::swap(recent[0], recent[1]);
        } else {
            if (older.pending != 0) table.add_occurrences(older.id, older.hash, older.pending);
            recent[1] = recent[0];
            recent[0] = { hash, table.count_occurrences(run_phrase.data(), len, hash, m, arena), 0 };
        }

        out.word_ids.push_back(uint32_t(recent[0].id));
        constexpr uint64_t word_bits = std::numeric_limits<uint64_t>::digits;

        if (m >= 2) {
            out.multi_bits |= uint64_t { 1 } << (out.num_runs % word_bits);
            append_varint(out.copies, m);
        }

        if (++out.num_runs % word_bits == 0) {
            out.multi_runs.push_back(out.multi_bits);
            out.multi_bits = 0;
        }

        out.length += m * (len - w);
        out.last_copies = m;
        if (m >= 2) out.repeated_ids.push_back(out.word_ids.back());
        m = 0;
    };

    // appends the last phrase run and adds the pending occurrences to the word table
    auto finish = [&]() {
        append_run();

        for (const counted_word& c : recent) {
            if (c.pending != 0) table.add_occurrences(c.id, c.hash, c.pending);
        }

        return max_sym;
    };

    // adds the phrase to the last phrase run or starts a new one with it
    auto add_phrase = [&]() {
        if (m == 0 || phrase.size() != run_phrase.size()
                || std::memcmp(phrase.data(), run_phrase.data(), phrase.size() * sizeof(sym_t)) != 0) {
            append_run();
            run_phrase.swap(phrase);
            phrase.assign(run_phrase.end() - w, run_phrase.end());
        }

        m++;
    };

    while (true) {
        uint64_t len;
        const sym_t* block = in.next(block_len, len);
        if (len == 0) break;
        const uint64_t range_len = in.range_length();
        if (range_len != text_ranges<sym_t>::unknown_length) range_end = 1 + range_len;
        sym_t min_blk = block[0], max_blk = block[0];

        for (uint64_t i = 1; i < len; i++) {
            min_blk = std::min(min_blk, block[i]);
            max_blk = std::max(max_blk, block[i]);
        }

        if (min_blk < min_symbol) throw std::invalid_argument("pfp_rlbwt: the text contains a symbol < min_symbol");
        max_sym = std::max<uint64_t>(max_sym, max_blk);
        uint64_t copied = 0; // the symbols block[0..copied-1] have been appended to phrase or are not in a phrase

        /* ends the current phrase with the trigger string that ends at block[i] and starts the next one with it;
           returns whether the next phrase starts after the range */
        auto trigger = [&](uint64_t i) {
            const uint64_t b = x + i + 1 - w; // start of the trigger string

            if (in_phrase) {
                phrase.insert(phrase.end(), block + copied, block + i + 1);
                add_phrase();
                if (b >= range_end) return true;
                phrase.erase(phrase.begin(), phrase.end() - w);
            } else {
                if (b >= range_end) return true;
                // the trigger string: the last w - i - 1 symbols before the block (if i + 1 < w) and block[..i]
                const uint64_t in_block = std::min(w, i + 1);
                phrase.assign(w_before.begin() + in_block, w_before.end());
                phrase.insert(phrase.end(), block + i + 1 - in_block, block + i + 1);
                in_phrase = true;
            }

            copied = i + 1;
            return false;
        };

        // the positions before check_from complete the first window; from in_block on, block[i - w] is removed from the
        // window
        const uint64_t first_check = std::min(len, check_from - std::min(check_from, x));
        const uint64_t in_block = std::max(first_check, std::min(len, w));
        for (uint64_t i = 0; i < first_check; i++) window.add(block[i], w_before[i]);

        for (uint64_t i = first_check; i < in_block; i++) {
            window.add(block[i], w_before[i]);
            if (window.is_trigger() && trigger(i)) return finish();
        }

        for (uint64_t i = in_block; i < len; i++) {
            window.add(block[i], block[i - w]);
            if (!window.is_trigger()) [[likely]] continue;
            if (trigger(i)) return finish();
        }

        if (in_phrase) phrase.insert(phrase.end(), block + copied, block + len);

        // the w symbols before the next block
        if (len >= w) {
            std::copy(block + len - w, block + len, w_before.begin());
        } else {
            std::move(w_before.begin() + len, w_before.end(), w_before.begin());
            std::copy(block, block + len, w_before.end() - len);
        }

        x += len;
    }

    // the last phrase ends with $^w
    if (in_phrase) {
        phrase.insert(phrase.end(), w, sym_t(dollar));
        add_phrase();
    }

    return finish();
}

template <typename sym_t>
inline uint64_t choose_trigger_seed(text_ranges<sym_t>& text, const params& prm)
{
    const uint16_t p = prm.num_threads;
    const uint64_t w = prm.window;
    std::vector<uint64_t> windows(p, 0), triggers(p, 0), retry_triggers(p, 0);
    const uint64_t block_len = block_symbols<sym_t>(); // evaluated in this thread (see current_level)

    run_on_threads(p, [&](uint64_t t, uint64_t team_size) {
        for (uint64_t i = t; i < p; i += team_size) {
            std::unique_ptr<typename text_ranges<sym_t>::reader> reader = text.open(i, nullptr);
            uint64_t len;
            const sym_t* block = reader->next(block_len, len);
            // the first block lies after the range only if the range is empty
            if (len == 0 || reader->range_length() == 0) continue;
            rolling_window_hash<sym_t> window(w, prm.modulus, 0);
            rolling_window_hash<sym_t> retry_window(w, prm.modulus, retry_trigger_seed);
            uint64_t num_triggers = 0;
            uint64_t num_retry_triggers = 0;

            for (uint64_t x = 0; x < len; x++) {
                const sym_t out = x >= w ? block[x - w] : sym_t(0);
                window.add(block[x], out);
                retry_window.add(block[x], out);
                if (x + 1 < w) continue;
                num_triggers += window.is_trigger();
                num_retry_triggers += retry_window.is_trigger();
            }

            windows[i] = len < w ? 0 : len + 1 - w;
            triggers[i] = num_triggers;
            retry_triggers[i] = num_retry_triggers;
        }
    });

    const uint64_t num_windows = std::accumulate(windows.begin(), windows.end(), uint64_t { 0 });
    const uint64_t num_triggers = std::accumulate(triggers.begin(), triggers.end(), uint64_t { 0 });
    const uint64_t num_retry_triggers = std::accumulate(retry_triggers.begin(), retry_triggers.end(), uint64_t { 0 });
    const bool degenerate = num_triggers < num_windows / (degenerate_parse_tolerance * prm.modulus);
    return degenerate && num_retry_triggers > num_triggers ? retry_trigger_seed : 0;
}

template <typename sym_t>
inline std::array<uint64_t, 256> count_symbols(word_table<sym_t>& table, uint64_t w, uint16_t p)
{
    const uint64_t D = table.size();
    std::array<uint64_t, 256> freq {};

    #pragma omp parallel num_threads(p)
    {
        std::array<uint64_t, 256> counts {}; // the counts of the words of this thread

        #pragma omp for
        for (uint64_t s = 0; s < D; s++) {
            const word_entry<sym_t>& e = table.entry(s);
            for (uint64_t k = 0; k + w < e.len; k++) counts[e.str[k]] += e.num_occ;
        }

        #pragma omp critical
        for (uint64_t c = 0; c < counts.size(); c++) freq[c] += counts[c];
    }

    freq[dollar] = 0;
    return freq;
}

template <typename sym_t>
inline void rename_symbols(word_table<sym_t>& table, const std::array<uint8_t, 256>& map, uint64_t& max_symbol, uint16_t p)
{
    const uint64_t D = table.size();
    uint64_t min_renamed = std::numeric_limits<uint64_t>::max();
    uint64_t max_renamed = 0;

    #pragma omp parallel for num_threads(p) reduction(min : min_renamed) reduction(max : max_renamed)
    for (uint64_t s = 0; s < D; s++) {
        const word_entry<sym_t>& e = table.entry(s);
        sym_t* word = const_cast<sym_t*>(e.str); // the words lie in the arenas of the parsing

        for (uint64_t k = 0; k < e.len; k++) {
            if (word[k] == dollar) continue;
            word[k] = sym_t(map[word[k]]);
            min_renamed = std::min<uint64_t>(min_renamed, word[k]);
            max_renamed = std::max<uint64_t>(max_renamed, word[k]);
        }
    }

    if (min_renamed < min_symbol) throw std::invalid_argument("pfp_rlbwt: a symbol is renamed to one < min_symbol");
    max_symbol = max_renamed;
}

template <typename sym_t>
inline parsing<sym_t> parse(text_ranges<sym_t>& text, const params& prm, uint64_t seed, const std::function<bool(const result&)>* parsed)
{
    const uint16_t p = prm.num_threads;
    parsing<sym_t> prs;
    prs.w = prm.window;
    prs.position_bytes = prm.position_bytes;
    prs.phrases_of_thread.resize(p);
    word_table<sym_t> table(p);
    std::vector<word_arena> arenas(p);
    std::vector<uint64_t> max_sym(p, 0);
    std::vector<uint64_t> range_len(p, 0); // [0..p-1] the number of symbols of range i

    phase_progress prog(prm.log, text.max_size());
    const uint64_t block_len = block_symbols<sym_t>(); // evaluated in this thread (see current_level)

    // in part i, the phrases starting in range i are parsed, in part 0 also the one starting at 0
    run_on_threads(p, [&](uint64_t t, uint64_t team_size) {
        for (uint64_t i = t; i < p; i += team_size) {
            std::unique_ptr<typename text_ranges<sym_t>::reader> in = text.open(i, &prog);
            max_sym[i] = parse_phrases_in_range<sym_t>(*in, i == 0, prm, seed, table, arenas[i],
                prs.phrases_of_thread[i], block_len);
            range_len[i] = in->range_length();
        }
    });

    prog.finish();
    prs.n = std::accumulate(range_len.begin(), range_len.end(), uint64_t { 0 });
    const uint64_t n = prs.n;
    for (uint16_t i = 0; i < p; i++) prs.max_symbol = std::max(prs.max_symbol, max_sym[i]);

    std::array<uint64_t, 256> freq {}; // the frequencies of the symbols of T, if text_ranges needs them
    bool counted = false;

    if constexpr (sizeof(sym_t) == 1) {
        if (text.needs_frequencies()) {
            freq = count_symbols(table, prm.window, p);
            counted = true;
        }
    }

    if (const std::array<uint8_t, 256>* map = text.ranges_read(counted ? &freq : nullptr); map != nullptr) {
        if constexpr (sizeof(sym_t) == 1) rename_symbols(table, *map, prs.max_symbol, p);
        else throw std::invalid_argument("pfp_rlbwt: only symbols of one byte can be renamed");
    }

    prs.continues.assign(p, 0);
    uint64_t next_first = word_table<sym_t>::max_words; // the word id of the first phrase run of the next thread

    for (uint64_t i = p; i-- > 0;) {
        thread_phrases& runs = prs.phrases_of_thread[i];
        if (runs.num_runs == 0) continue;
        prs.continues[i] = runs.word_ids.back() == next_first;
        if (prs.continues[i]) runs.repeated_ids.push_back(runs.word_ids.back());
        next_first = runs.word_ids.front();
    }

    for (uint64_t i = 0; i < p; i++) {
        thread_phrases& runs = prs.phrases_of_thread[i];
        prs.num_phrase_runs += runs.num_runs - prs.continues[i];
        if (runs.num_runs % std::numeric_limits<uint64_t>::digits != 0) runs.multi_runs.push_back(runs.multi_bits);
        runs.word_ids.shrink_to_fit();
        runs.multi_runs.shrink_to_fit();
        runs.copies.shrink_to_fit();
    }

    dictionary<sym_t>& dict = prs.dict;
    const uint64_t D = table.size();
    dict.num_words = D;
    dict.size = 0;

    for (uint64_t s = 0; s < D; s++) {
        dict.size += table.entry(s).len + 1;
        prs.num_phrases += table.entry(s).num_occ;
    }

    if (parsed != nullptr) {
        result sizes;
        sizes.n_with_terminator = n + 1;
        sizes.num_phrases = prs.num_phrases;
        sizes.num_phrase_runs = prs.num_phrase_runs;
        sizes.num_words = D;
        sizes.dictionary_size = dict.size;
        if (!(*parsed)(sizes)) return prs;
    }

    // sort the words
    malloc_array_t<uint32_t> order = new_malloc_array<uint32_t>(D);

    #pragma omp parallel for num_threads(p)
    for (uint64_t i = 0; i < D; i++) order[i] = uint32_t(i);

    comparison_sort(order.get(), order.get() + D, [&](uint32_t a, uint32_t b) {
        const word_entry<sym_t>& x = table.entry(a);
        const word_entry<sym_t>& y = table.entry(b);
        return std::lexicographical_compare(x.str, x.str + x.len, y.str, y.str + y.len);
    }, p);

    // store the words in lexicographic order
    prs.word_rank = new_malloc_array<uint32_t>(D);
    dict.word_start = uint_array(D + 1, value_bytes(dict.size));
    dict.num_occ = uint_array(D, value_bytes(prs.num_phrases));
    uint64_t pos = 0;

    for (uint64_t s = 0; s < D; s++) {
        const word_entry<sym_t>& e = table.entry(order[s]);
        dict.word_start.set(s, pos);
        dict.num_occ.set(s, e.num_occ);
        pos += e.len + 1;
    }

    dict.word_start.set(D, pos);
    dict.text = new_malloc_array<sym_t>(dict.size);

    #pragma omp parallel for num_threads(p)
    for (uint64_t s = 0; s < D; s++) {
        const word_entry<sym_t>& e = table.entry(order[s]);
        std::memcpy(dict.text.get() + dict.word_start[s], e.str, e.len * sizeof(sym_t));
        dict.text[dict.word_start[s + 1] - 1] = separator;
        prs.word_rank[order[s]] = uint32_t(s);
    }

    for (thread_phrases& runs : prs.phrases_of_thread) {
        for (uint32_t id : runs.repeated_ids) prs.repeated_words.push_back(prs.word_rank[id]);
        runs.repeated_ids = std::vector<uint32_t>();
    }

    std::sort(prs.repeated_words.begin(), prs.repeated_words.end());
    prs.repeated_words.erase(std::unique(prs.repeated_words.begin(), prs.repeated_words.end()),
        prs.repeated_words.end());

    return prs;
}

template <typename idx_t, typename sym_t>
inline void store_parse(phrase_arrays<idx_t, sym_t>& arr, uint64_t J, uint16_t threads)
{
    // P[i] = 1 + the letter of phrase run i
    arr.P = new_malloc_array<idx_t>(J + 1);
    arr.P[J] = 0;

    for (uint64_t u = 0; u + 1 < arr.run_offsets.size(); u++) {
        const uint64_t first = arr.run_offsets[u], size = arr.run_offsets[u + 1] - first;

        arr.runs_of_part[u].word_ids.move_to([&](uint64_t j, uint32_t id) {
            if (j < size) arr.P[first + j] = idx_t(1 + arr.letter_of(first + j, id));
        }, threads);
    }

    arr.free_word_ids();
}

template <typename idx_t, typename sym_t>
inline phrase_arrays<idx_t, sym_t> build_phrase_arrays(parsing<sym_t>& prs, const std::string& prefix, uint16_t threads, bool recursive)
{
    const uint64_t J = prs.num_phrase_runs;
    const uint64_t w = prs.w;
    const dictionary<sym_t>& dict = prs.dict;
    const uint64_t D = dict.num_words;
    const uint8_t bytes = prs.position_bytes;
    const uint64_t parts = prs.phrases_of_thread.size();
    const std::vector<uint64_t>& repeated = prs.repeated_words;
    const uint64_t R = repeated.size();
    constexpr uint64_t word_bits = std::numeric_limits<uint64_t>::digits;
    phrase_arrays<idx_t, sym_t> arr;
    arr.position_bytes = bytes;
    if (R != 0) arr.repeated_set = ranked_bits([&](uint64_t y) { return repeated[y]; }, R, D);
    const ranked_bits& is_repeated = arr.repeated_set;
    // the number of phrase runs of part u
    auto runs_of = [&](uint64_t u) { return prs.phrases_of_thread[u].num_runs; };
    // the symbol before the last w symbols of the word with id `id`
    auto symbol = [&](uint64_t id) {
        const uint64_t s = prs.word_rank[id];
        return dict.text[dict.word_start[s] + dict.length(s) - w - 1];
    };

    /* per word id: [0..D-1] id_advance[id] = the advance of the word; bit id % 64 of id_repeated[id / 64] = 1 <=> the
       word is repeated */
    uint64_t max_advance = 0;
    for (uint64_t s = 0; s < D; s++) max_advance = std::max(max_advance, dict.length(s) - w);
    uint_array id_advance(D, value_bytes(max_advance));
    std::vector<uint64_t> id_repeated(D / word_bits + 1, 0);

    #pragma omp parallel for num_threads(threads)
    for (uint64_t b = 0; b < id_repeated.size(); b++) {
        for (uint64_t id = word_bits * b; id < std::min(D, word_bits * (b + 1)); id++) {
            const uint64_t s = prs.word_rank[id];
            id_advance.set(id, dict.length(s) - w);
            if (R != 0 && is_repeated.contains(s)) id_repeated[b] |= uint64_t { 1 } << (id % word_bits);
        }
    }

    auto advance = [&](uint64_t id) { return id_advance[id]; };

    /* the phrase runs of each part are stored after those of the part before: [0..parts] offsets[u] = number of runs
       the parts before part u store; the last run of a part whose last phrase run continues in the next part with
       phrase runs is not stored: that part stores the run as its first one */
    std::vector<uint64_t> offsets(parts + 1, 0);
    for (uint64_t u = 0; u < parts; u++) offsets[u + 1] = offsets[u] + runs_of(u) - prs.continues[u];
    /* [0..parts] length_before[u] = the total length of the phrase runs of the parts before part u; [0..parts-1]
       last_copies[u] = the number of phrases of the last phrase run of part u */
    std::vector<uint64_t> length_before(parts + 1, 0), last_copies(parts, 0);

    for (uint64_t u = 0; u < parts; u++) {
        length_before[u + 1] = prs.phrases_of_thread[u].length;
        last_copies[u] = prs.phrases_of_thread[u].last_copies;
    }

    for (uint64_t u = 0; u < parts; u++) length_before[u + 1] += length_before[u];

    /* [0..parts-1] carry[u] = the number of phrases of the parts before part u in the phrase run that part u stores as
       its first one; next_rank[u] = the rank of the word of the first phrase run after part u (no_rank: none) */
    constexpr uint64_t no_rank = std::numeric_limits<uint64_t>::max();
    std::vector<uint64_t> carry(parts, 0), next_rank(parts, no_rank);

    for (uint64_t u = 0, c = 0; u < parts; u++) {
        if (runs_of(u) == 0) continue;
        carry[u] = c;
        c = prs.continues[u] ? last_copies[u] + (runs_of(u) == 1 ? c : 0) : 0;
    }

    for (uint64_t u = parts, after = no_rank; u-- > 0;) {
        next_rank[u] = after;
        if (runs_of(u) != 0) after = prs.word_rank[prs.phrases_of_thread[u].word_ids.front()];
    }

    for (uint64_t u = 0; u < parts; u++) {
        if (runs_of(u) != 0) arr.last_symbol = symbol(prs.phrases_of_thread[u].word_ids.back());
    }

    /* the ends are written before P is allocated: the end of a run is the end of the run before it plus its number of
       phrases times the advance of its word; the phrase runs of the repeated words are collected per part as (i, index
       in repeated, key), and the extra phrases of the repeated words per part */
    std::vector<std::vector<std::tuple<uint64_t, uint64_t, uint64_t>>> repeated_runs(parts);
    std::vector<std::vector<uint64_t>> extra_of_part(parts, std::vector<uint64_t>(R, 0));

    {
        arr.ends_file = temporary_array(prefix + ext_phrase_run_ends, temporary_arrays_in_memory(J * bytes));
        const std::unique_ptr<positioned_output> ends_out = arr.ends_file.output(J * bytes);

        run_on_threads(threads, [&](uint64_t t, uint64_t team_size) {
            // the ends of a part are collected in buffers of about run_buffer_bytes bytes
            const uint64_t buffer_runs = std::max<uint64_t>(1, run_buffer_bytes / bytes);
            std::vector<uint8_t> ends(buffer_runs * bytes + sizeof(uint64_t));

            for (uint64_t u = t; u < parts; u += team_size) {
                thread_phrases& runs = prs.phrases_of_thread[u];
                word_id_vector::reader ids(runs.word_ids);
                block_vector<uint64_t>::reader multi_runs(runs.multi_runs);
                uint64_t multi_bits = 0; // the bits of multi_runs of the next run and the runs after it in its word
                block_vector<uint8_t>::reader copies(runs.copies);
                const uint64_t size = offsets[u + 1] - offsets[u];
                uint64_t end = w - 1 + length_before[u];
                uint64_t id = size == 0 ? 0 : ids.next(); // the id of the next run

                for (uint64_t first = 0; first < size; first += buffer_runs) {
                    const uint64_t num = std::min(buffer_runs, size - first);

                    for (uint64_t j = 0; j < num; j++) {
                        const uint64_t i = first + j; // the index of the run in the part
                        if (i % word_bits == 0) multi_bits = multi_runs.next();
                        const uint64_t m = (multi_bits >> (i % word_bits)) & 1 ? read_varint(copies) : 1;
                        end += m * advance(id);
                        std::memcpy(ends.data() + j * bytes, &end, bytes);
                        const uint64_t next_id = i + 1 < runs_of(u) ? ids.next() : 0;

                        if ((id_repeated[id / word_bits] >> (id % word_bits)) & 1) [[unlikely]] {
                            // the run with all its phrases (also those of the parts before), and the word after it
                            const uint64_t s = prs.word_rank[id];
                            const uint64_t all = m + (i == 0 ? carry[u] : 0);
                            const uint64_t after = i + 1 < runs_of(u) ? prs.word_rank[next_id] : next_rank[u];
                            const uint64_t y = is_repeated.rank(s);
                            const bool down = after == no_rank || after < s;
                            repeated_runs[u].emplace_back(offsets[u] + i, y, run_key(down, all));
                            extra_of_part[u][y] += all - 1;
                        }

                        id = next_id;
                    }

                    ends_out->write(ends.data(), num * bytes, (offsets[u] + first) * bytes);
                }

                runs.multi_runs = block_vector<uint64_t>();
                runs.copies = block_vector<uint8_t>();
            }
        });

        ends_out->close();
    }

    // the letters of the repeated words: one per distinct key of their phrase runs
    if (R != 0) {
        std::vector<std::vector<uint64_t>> keys(R);

        for (uint64_t u = 0; u < parts; u++) {
            for (const auto& [i, y, key] : repeated_runs[u]) keys[y].push_back(key);
        }

        arr.repeated_words.resize(R);
        arr.letters_before.assign(R + 1, 0);

        for (uint64_t y = 0; y < R; y++) {
            repeated_word& rw = arr.repeated_words[y];
            std::sort(keys[y].begin(), keys[y].end());
            keys[y].erase(std::unique(keys[y].begin(), keys[y].end()), keys[y].end());
            rw.word = repeated[y];
            rw.first_letter = rw.word + arr.letters_before[y];
            rw.keys = std::move(keys[y]);
            for (uint64_t u = 0; u < parts; u++) rw.extra_phrases += extra_of_part[u][y];
            arr.letters_before[y + 1] = arr.letters_before[y] + rw.keys.size() - 1;
        }

        for (uint64_t u = 0; u < parts; u++) {
            for (const auto& [i, y, key] : repeated_runs[u]) {
                const repeated_word& rw = arr.repeated_words[y];
                const uint64_t index = std::lower_bound(rw.keys.begin(), rw.keys.end(), key) - rw.keys.begin();
                arr.repeated_run_letters.emplace_back(i, rw.first_letter + index);
            }

            repeated_runs[u] = std::vector<std::tuple<uint64_t, uint64_t, uint64_t>>();
        }
    }

    arr.num_letters = D + (R == 0 ? 0 : arr.letters_before[R]);
    id_advance = uint_array();

    // the letters per word id (num_letters for the repeated words, whose phrase runs have letters of their own)
    arr.id_letter = uint_array(D, value_bytes(arr.num_letters));

    #pragma omp parallel for num_threads(threads)
    for (uint64_t id = 0; id < D; id++) {
        const uint64_t s = prs.word_rank[id];
        const bool rep = R != 0 && is_repeated.contains(s);
        arr.id_letter.set(id, rep ? arr.num_letters : R == 0 ? s : s + arr.letters_before[is_repeated.rank(s)]);
    }

    arr.runs_of_part = std::move(prs.phrases_of_thread);
    arr.word_rank = std::move(prs.word_rank);
    arr.run_offsets = std::move(offsets);
    if (!recursive) store_parse(arr, J, threads);
    return arr;
}

template <typename idx_t, typename sym_t, typename text_t>
inline void letter_source<idx_t, sym_t, text_t>::letter_reader::read(uint64_t pos, uint64_t len, text_t* dst)
{
    // the part storing phrase run pos
    uint64_t u = std::upper_bound(arr.run_offsets.begin(), arr.run_offsets.end(), pos)
        - arr.run_offsets.begin() - 1;
    auto ids = std::make_unique<word_id_vector::reader>(arr.runs_of_part[u].word_ids,
        pos - arr.run_offsets[u]);

    for (uint64_t k = 0; k < len; k++) {
        const uint64_t i = pos + k;

        while (i == arr.run_offsets[u + 1]) {
            u++;
            ids = std::make_unique<word_id_vector::reader>(arr.runs_of_part[u].word_ids, 0);
        }

        dst[k] = text_t(2 + arr.letter_of(i, ids->next()));
    }
}

} // namespace pfp_rlbwt

#endif // PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_PARSE_TPP
