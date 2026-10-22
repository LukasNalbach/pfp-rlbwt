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

#ifndef PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_PARSE_HPP
#define PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_PARSE_HPP

#include <algorithm>
#include <array>
#include <atomic>
#include <exception>
#include <limits>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <misc/sort.hpp>

#include <omp.h>

#include <misc/memory.hpp>
#include <pfp_rlbwt/construction/common.hpp>
#include <pfp_rlbwt/construction/phrase_runs.hpp>
#include <pfp_rlbwt/types.hpp>

namespace pfp_rlbwt {

/* The prefix-free parsing of S = $ T $^w: a window S[b..b+w-1] with 1 <= b and b+w-1 <= n is a trigger string if
   its hash is selected (see rolling_window_hash); with b_1 < ... < b_{K-1} the starts of the trigger strings,
   b_0 = 0 and b_K = n+1, phrase k is S[b_k..e_k] with e_k = b_{k+1}+w-1. Consecutive phrases overlap by w symbols.
   Each position x in [0..n] of S lies in exactly one range [b_k..e_k-w]; the suffix of S starting at x corresponds
   to the suffix of T$ starting at x-1 (x >= 1) and to the suffix $ (x = 0). A phrase run is a maximal run of equal
   consecutive phrases; the parsing stores the J phrase runs, each with the end of its last phrase. Phrase run i
   consists of m_i phrases equal to a word of length w + a, where a is the advance of the word, and ends at
   F_i = F_{i-1} + m_i * a (F_{-1} = w-1). */

// the dictionary: the distinct phrases (words) in lexicographic order
template <typename sym_t>
struct dictionary {
    uint64_t num_words = 0; // D, the number of words
    uint64_t size = 0; // total length of the words and their separators
    // the words: [0..size-1] the words in lexicographic order, each followed by separator
    malloc_array_t<sym_t> text;
    // the starting positions of the words: [0..D] word_start[s] = position of word s in text; word_start[D] = size
    uint_array word_start;
    // the numbers of occurrences of the words: [0..D-1] num_occ[s] = number of phrases equal to word s
    uint_array num_occ;

    // returns the length of word s
    inline uint64_t length(uint64_t s) const { return word_start[s + 1] - word_start[s] - 1; }
};

/* the phrase runs one thread has found, in text order; for j in [0..word_ids.size()-1], entry j belongs to the j-th
   of them, phrases k..k'; aligned to a cache line (the threads append to different ones) */
template <typename sym_t>
struct alignas(cache_line_bytes) thread_phrases {
    // word_ids[j] = id of the word equal to the phrases k..k' (see word_table)
    block_vector<uint32_t> word_ids;
    // symbol_before_next[j] = S[e_k-w], the symbol before phrase k+1 (the same for each phrase of the run)
    block_vector<sym_t> symbol_before_next;
    // the ends: bytes [j * position_bytes..(j + 1) * position_bytes - 1] store e_k' (little endian, see parsing)
    block_vector<uint8_t> ends;
    // the ids of the words of the phrase runs with at least two phrases (with repetitions)
    std::vector<uint32_t> repeated_ids;
};

// the parsing as the threads have found it
template <typename sym_t>
struct parsing {
    uint64_t n = 0; // n, the length of T
    uint64_t w = 0; // w, the length of a trigger string
    uint64_t num_phrases = 0; // K, the number of phrases
    uint64_t num_phrase_runs = 0; // J, the number of phrase runs
    uint64_t max_symbol = 0; // the largest symbol of T
    uint8_t position_bytes = 0; // bytes per position of S in the arrays (see uint_array)
    dictionary<sym_t> dict;
    // the phrases found by the threads: [0..p-1] phrases_of_thread[i] = the phrase runs found by thread i
    std::vector<thread_phrases<sym_t>> phrases_of_thread;
    /* [0..p-1] continues[i] = 1 <=> the last phrase run of thread i continues in the next thread with phrase runs
       (that thread stores the run) */
    std::vector<char> continues;
    // the ranks of the words with a phrase run of at least two phrases, in increasing order
    std::vector<uint64_t> repeated_words;
    // the ranks of the words: [0..D-1] word_rank[id] = lexicographic rank of the word with id `id` (see word_table)
    malloc_array_t<uint32_t> word_rank;
};

/* a word with a phrase run of at least two phrases; its phrase runs have the letters first_letter, ...,
   first_letter + keys.size() - 1 (see build_phrase_arrays) */
struct repeated_word {
    uint64_t word = 0; // the rank of the word
    uint64_t extra_phrases = 0; // the number of its phrases minus the number of its phrase runs
    uint64_t first_letter = 0;
    std::vector<uint64_t> keys; // the distinct keys of its phrase runs in increasing order (see run_key)
};

// the phrase runs as arrays
template <typename idx_t, typename sym_t>
struct phrase_arrays {
    /* the parse of the phrase runs: [0..J] P[i] = 1 + letter of phrase run i; P[J] = 0. The letters are ordered by
       word; the phrase runs of a word in repeated_words get one letter per key (see run_key), the others one
       letter per word */
    malloc_array_t<idx_t> P;
    uint64_t num_letters = 0; // the number of letters
    // [0..J-1] symbol_before_next[i] = the symbol before the phrase after a phrase of run i
    malloc_array_t<sym_t> symbol_before_next;
    // the ends of the phrase runs: [0..J-1] ends[i] = F_i, the end of the last phrase of run i
    uint_array ends;
    // the words with a phrase run of at least two phrases, in increasing order
    std::vector<repeated_word> repeated_words;
    ranked_bits repeated_set; // the set of these words
};

/* returns the key of a phrase run with m phrases that is followed by a smaller (down) or larger suffix: the keys of
   the down runs increase with m, those of the up runs decrease with m and are larger */
inline uint64_t run_key(bool down, uint64_t m) { return down ? m : std::numeric_limits<uint64_t>::max() - m; }

// floor(2^64 / golden ratio), the multiplier of Fibonacci hashing
constexpr uint64_t golden_ratio_64 = 0x9E3779B97F4A7C15;

/* a parse with fewer than n / (degenerate_parse_tolerance * modulus) phrases is degenerate (it is expected to have about
   n / modulus phrases); the seed of the trigger strings (see rolling_window_hash) is then switched between 0 and
   retry_trigger_seed: before the parse on a sample (see choose_trigger_seed), and after it, if the parse is degenerate */
constexpr uint64_t degenerate_parse_tolerance = 2;
constexpr uint64_t retry_trigger_seed = golden_ratio_64;

// number of symbols read from the text at once
constexpr uint64_t read_block_size = uint64_t { 1 } << 20;

// number of bytes of the first and of the largest blocks of a word_arena
constexpr uint64_t min_arena_block_size = uint64_t { 1 } << 16;
constexpr uint64_t max_arena_block_size = uint64_t { 1 } << 20;

// base-2 logarithm of the number of entries per block of the entries of the words (see word_table)
constexpr uint8_t log_entries_per_block = 16;

// number of shards of the word table per thread
constexpr uint64_t shards_per_thread = 64;

// initial number of slots of a shard
constexpr uint64_t shard_initial_slots = 64;

// spreads the bits of x over all bits of the result
inline uint64_t mix(uint64_t x)
{
    x ^= x >> 32;
    x *= golden_ratio_64;
    return x ^ (x >> 32);
}

// returns a hash of the bytes s[0..len-1]
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

/**
 * @brief the polynomial hash (mod 2^64, base golden_ratio_64) of the last w symbols of S, updated symbol by symbol;
 *        the window is a trigger string if the mixed hash, xor the seed, is smaller than UINT64_MAX / modulus
 */
template <typename sym_t>
class rolling_window_hash {
private:
    uint64_t w; // w, the window length
    // the last w symbols: ring[x & mask] = S[x] for the last w positions x added; bit_ceil(w) entries
    std::vector<sym_t> ring;
    uint64_t mask; // bit_ceil(w) - 1
    uint64_t threshold; // UINT64_MAX / modulus
    uint64_t pow_w = 1; // golden_ratio_64^w
    uint64_t hash = 0; // the hash of the window
    uint64_t seed; // selects the trigger strings (see is_trigger)

public:
    rolling_window_hash(uint64_t w, uint64_t modulus, uint64_t seed)
        : w(w), ring(std::bit_ceil(w), 0), mask(std::bit_ceil(w) - 1), threshold(UINT64_MAX / modulus), seed(seed)
    {
        for (uint64_t i = 0; i < w; i++) pow_w *= golden_ratio_64;
    }

    // appends S[pos] = c to the window (the windows must be appended at consecutive positions)
    inline void add(uint64_t pos, sym_t c)
    {
        const sym_t out = ring[(pos - w) & mask];
        ring[pos & mask] = c;
        hash = hash * golden_ratio_64 + uint64_t(c) - uint64_t(out) * pow_w;
    }

    // returns whether the window is a trigger string
    inline bool is_trigger() const { return mix(hash ^ seed) < threshold; }

    // replaces dst with the window S[pos-w+1..pos]
    void copy_to(uint64_t pos, std::vector<sym_t>& dst) const
    {
        dst.clear();
        for (uint64_t x = pos + 1 - w; x <= pos; x++) dst.push_back(ring[x & mask]);
    }
};

// the strings of the words a thread has added to the dictionary, in blocks that grow from min_arena_block_size to
// max_arena_block_size bytes
class word_arena {
private:
    /* the blocks: blocks[b] has min(min_arena_block_size * 2^b, max_arena_block_size) bytes, or the length of the
       word that did not fit into the block before if that is larger */
    std::vector<malloc_array_t<uint8_t>> blocks;
    uint8_t* free_begin = nullptr; // the first free byte of the last block
    uint64_t num_free = 0; // number of free bytes in the last block
    uint64_t block_size = 0; // size of the last block (unless a word did not fit into it)

public:
    // copies the bytes s[0..len-1] into the arena and returns the copy
    const void* store(const void* s, uint64_t len)
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
};

// a word and its number of occurrences
template <typename sym_t>
struct word_entry {
    const sym_t* str; // the symbols of the word (in a word_arena)
    uint64_t len; // its length
    uint64_t num_occ; // its number of occurrences
};

/**
 * @brief the words found by the parsing threads: each distinct word gets an id (0, 1, ... in the order in which the
 *        words are found) and its occurrences are counted. A word is looked up by its hash in one of the shards,
 *        each a hash table (open addressing, at most half full) with its own mutex.
 */
template <typename sym_t>
class word_table {
public:
    // ############################# TYPES AND CONSTANTS #############################

    using entry_t = word_entry<sym_t>;

    // the largest number of words (an id + 1 is stored in 32 bits)
    static constexpr uint64_t max_words = std::numeric_limits<uint32_t>::max();

private:
    static constexpr uint64_t entry_mask = (uint64_t { 1 } << log_entries_per_block) - 1;

    struct slot {
        uint64_t hash; // the hash of the word
        uint32_t id_plus_1; // 0 <=> the slot is empty
    };

    struct shard {
        std::mutex lock;
        // the hash table (linear probing, starting at slot h % slots.size() for hash h); slots.size() is a power of 2
        std::vector<slot> slots;
        uint64_t used = 0; // number of non-empty slots
    };

    // ############################# DATA MEMBERS #############################

    // the shards: shards[h >> shard_shift] stores the words with hash h; a power of two of them
    std::vector<shard> shards;
    uint8_t shard_shift; // 64 - log2(shards.size()); 64 if there is one shard
    /* the entries of the words in blocks of 2^log_entries_per_block entries: blocks[b][j] = entry of the word with
       id b * 2^log_entries_per_block + j; a block is allocated with its first entry (nullptr before) */
    std::vector<std::atomic<entry_t*>> blocks;
    std::mutex blocks_lock; // locked while a block is allocated
    std::atomic<uint64_t> next_id { 0 }; // the id of the next new word

    // ############################# INTERNAL METHODS #############################

    // returns the block of entries containing the entry of id, allocating it if necessary
    entry_t* block_of(uint64_t id)
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

    // doubles the number of slots of a shard
    static void grow(shard& sh)
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

public:
    // ############################# CONSTRUCTORS #############################

    explicit word_table(uint16_t threads)
        : shards(std::bit_ceil(uint64_t { threads } * shards_per_thread)),
          shard_shift(64 - std::bit_width(shards.size() - 1)),
          blocks((max_words >> log_entries_per_block) + 1)
    {
        for (shard& sh : shards) sh.slots.assign(shard_initial_slots, slot { 0, 0 });
        for (auto& b : blocks) b.store(nullptr, std::memory_order_relaxed);
    }

    ~word_table()
    {
        for (auto& b : blocks) delete[] b.load(std::memory_order_relaxed);
    }

    // ############################# METHODS #############################

    // returns the number of words
    uint64_t size() const { return next_id.load(); }

    // returns the entry of word id
    inline entry_t& entry(uint64_t id)
    {
        return blocks[id >> log_entries_per_block].load(std::memory_order_relaxed)[id & entry_mask];
    }

    // returns whether the word with id `id` is s[0..len-1] (id must have been returned to the calling thread)
    inline bool is_word(uint64_t id, const sym_t* s, uint64_t len)
    {
        const entry_t& e = entry(id);
        return e.len == len && std::memcmp(e.str, s, len * sizeof(sym_t)) == 0;
    }

    // adds count occurrences to the word with id `id` and hash `hash`
    void add_occurrences(uint64_t id, uint64_t hash, uint64_t count)
    {
        shard& sh = shards[shard_shift == 64 ? 0 : hash >> shard_shift];
        std::lock_guard<std::mutex> guard(sh.lock);
        entry(id).num_occ += count;
    }

    /**
     * @brief counts occurrences of the word s[0..len-1]; a new word is stored in arena
     * @param s the word
     * @param len its length
     * @param hash its hash
     * @param count the number of occurrences
     * @param arena the arena of the calling thread
     * @return the id of the word
     */
    uint32_t count_occurrences(const sym_t* s, uint64_t len, uint64_t hash, uint64_t count, word_arena& arena)
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
};

/**
 * @brief parses the phrases of S starting in [lo, hi) and appends their phrase runs to out
 * @param text the text T
 * @param lo first starting position of S
 * @param hi starting position of S after the last one
 * @param first true <=> the phrase starting at 0 is parsed too
 * @param prm parameters (window, modulus, position_bytes)
 * @param seed the seed of the trigger strings (see rolling_window_hash)
 * @param table the word table
 * @param arena the arena of the calling thread
 * @param out the phrases of the calling thread
 * @return the largest symbol read
 */
template <typename sym_t>
uint64_t parse_phrases_in_range(const text_source<sym_t>& text, uint64_t lo, uint64_t hi, bool first,
    const params& prm, uint64_t seed, word_table<sym_t>& table, word_arena& arena, thread_phrases<sym_t>& out)
{
    const uint64_t n = text.size();
    const uint64_t w = prm.window;
    std::unique_ptr<typename text_source<sym_t>::reader> reader = text.new_reader();
    malloc_array_t<sym_t> block = new_malloc_array<sym_t>(read_block_size);
    rolling_window_hash<sym_t> window(w, prm.modulus, seed);
    std::vector<sym_t> phrase; // S[b_k..x], where x is the position read last
    bool in_phrase = first;
    if (first) phrase.push_back(dollar);
    uint64_t x = first ? 1 : lo; // next position of S to read
    const uint64_t check_from = x + w - 1; // the first position at which the window has w symbols
    uint64_t max_sym = 0;

    std::vector<sym_t> run_phrase; // the phrase of the last phrase run, if it is not yet appended to out
    uint64_t run_copies = 0; // the number of phrases of that run (0: there is none)
    uint64_t run_end = 0; // the end of its last phrase

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
        if (run_copies == 0) return;
        const uint64_t len = run_phrase.size();
        const uint64_t hash = string_hash((const uint8_t*) run_phrase.data(), len * sizeof(sym_t));
        counted_word& older = recent[1];

        if (older.hash == hash && older.id != word_table<sym_t>::max_words
                && table.is_word(older.id, run_phrase.data(), len)) {
            older.pending += run_copies;
            std::swap(recent[0], recent[1]);
        } else {
            if (older.pending != 0) table.add_occurrences(older.id, older.hash, older.pending);
            recent[1] = recent[0];
            recent[0] = { hash, table.count_occurrences(run_phrase.data(), len, hash, run_copies, arena), 0 };
        }

        out.word_ids.push_back(uint32_t(recent[0].id));
        out.symbol_before_next.push_back(run_phrase[len - w - 1]);
        for (uint8_t b = 0; b < prm.position_bytes; b++) out.ends.push_back(uint8_t(run_end >> (8 * b)));
        if (run_copies >= 2) out.repeated_ids.push_back(out.word_ids.back());
        run_copies = 0;
    };

    // appends the last phrase run and adds the pending occurrences to the word table
    auto finish = [&]() {
        append_run();

        for (const counted_word& c : recent) {
            if (c.pending != 0) table.add_occurrences(c.id, c.hash, c.pending);
        }

        return max_sym;
    };

    // adds the phrase ending at position end of S to the last phrase run or starts a new one with it
    auto add_phrase = [&](uint64_t end) {
        if (run_copies == 0 || phrase.size() != run_phrase.size()
                || std::memcmp(phrase.data(), run_phrase.data(), phrase.size() * sizeof(sym_t)) != 0) {
            append_run();
            run_phrase.swap(phrase);
            phrase.assign(run_phrase.end() - w, run_phrase.end());
        }

        run_copies++;
        run_end = end;
    };

    while (x <= n) {
        const uint64_t len = std::min(read_block_size, n + 1 - x);
        reader->read(x - 1, len, block.get());
        const auto [min_blk, max_blk] = std::minmax_element(block.get(), block.get() + len);
        if (*min_blk < min_symbol) throw std::invalid_argument("pfp_rlbwt: the text contains a symbol < min_symbol");
        max_sym = std::max<uint64_t>(max_sym, *max_blk);
        uint64_t copied = 0; // the symbols block[0..copied-1] have been appended to phrase or are not in a phrase

        for (uint64_t i = 0; i < len; i++) {
            const uint64_t pos = x + i;
            window.add(pos, block[i]);
            if (pos < check_from || !window.is_trigger()) [[likely]] continue;
            const uint64_t b = pos + 1 - w; // start of the trigger string

            if (in_phrase) {
                phrase.insert(phrase.end(), block.get() + copied, block.get() + i + 1);
                add_phrase(pos);

                if (b >= hi) return finish();
                phrase.erase(phrase.begin(), phrase.end() - w);
            } else {
                if (b >= hi) return finish();
                window.copy_to(pos, phrase);
                in_phrase = true;
            }

            copied = i + 1;
        }

        if (in_phrase) phrase.insert(phrase.end(), block.get() + copied, block.get() + len);
        x += len;
    }

    // the last phrase ends with $^w
    if (in_phrase) {
        phrase.insert(phrase.end(), w, sym_t(dollar));
        add_phrase(n + w);
    }

    return finish();
}

/**
 * @brief chooses the seed of the trigger strings (see rolling_window_hash) from a sample of T, the first read_block_size
 *        symbols of each of the p ranges the text is parsed in: 0, unless the sample has fewer than
 *        1 / (degenerate_parse_tolerance * modulus) trigger strings per window with it and more with retry_trigger_seed
 * @param text the text T
 * @param prm parameters (window, modulus, num_threads)
 * @return the seed
 */
template <typename sym_t>
uint64_t choose_trigger_seed(const text_source<sym_t>& text, const params& prm)
{
    const uint16_t p = prm.num_threads;
    const uint64_t n = text.size();
    const uint64_t w = prm.window;
    std::vector<uint64_t> windows(p, 0), triggers(p, 0), retry_triggers(p, 0);

    run_on_threads(p, [&](uint64_t t, uint64_t team_size) {
        malloc_array_t<sym_t> block = new_malloc_array<sym_t>(read_block_size);
        std::unique_ptr<typename text_source<sym_t>::reader> reader = text.new_reader();

        for (uint64_t i = t; i < p; i += team_size) {
            const uint64_t lo = (i * n) / p;
            const uint64_t len = std::min(read_block_size, ((i + 1) * n) / p - lo);
            reader->read(lo, len, block.get());
            rolling_window_hash<sym_t> window(w, prm.modulus, 0);
            rolling_window_hash<sym_t> retry_window(w, prm.modulus, retry_trigger_seed);
            uint64_t num_triggers = 0;
            uint64_t num_retry_triggers = 0;

            for (uint64_t x = 0; x < len; x++) {
                window.add(x, block[x]);
                retry_window.add(x, block[x]);
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

/**
 * @brief computes the prefix-free parsing of S = $ T $^w and its dictionary
 * @param text the text T
 * @param prm parameters
 * @param seed the seed of the trigger strings (see rolling_window_hash)
 * @return the parsing
 */
template <typename sym_t>
parsing<sym_t> parse(const text_source<sym_t>& text, const params& prm, uint64_t seed)
{
    const uint16_t p = prm.num_threads;
    const uint64_t n = text.size();
    parsing<sym_t> prs;
    prs.n = n;
    prs.w = prm.window;
    prs.position_bytes = prm.position_bytes;
    prs.phrases_of_thread.resize(p);
    word_table<sym_t> table(p);
    std::vector<word_arena> arenas(p);
    std::vector<uint64_t> max_sym(p, 0);

    // part i parses the phrases starting in [1 + i*n/p, 1 + (i+1)*n/p) and part 0 also the one starting at 0
    run_on_threads(p, [&](uint64_t t, uint64_t team_size) {
        for (uint64_t i = t; i < p; i += team_size) {
            max_sym[i] = parse_phrases_in_range(text, 1 + (i * n) / p, 1 + ((i + 1) * n) / p, i == 0, prm, seed,
                table, arenas[i], prs.phrases_of_thread[i]);
        }
    });

    for (uint16_t i = 0; i < p; i++) prs.max_symbol = std::max(prs.max_symbol, max_sym[i]);
    prs.continues.assign(p, 0);
    uint64_t next_first = word_table<sym_t>::max_words; // the word id of the first phrase run of the next thread

    for (uint64_t i = p; i-- > 0;) {
        thread_phrases<sym_t>& b = prs.phrases_of_thread[i];
        if (b.word_ids.size() == 0) continue;
        prs.continues[i] = b.word_ids.back() == next_first;
        if (prs.continues[i]) b.repeated_ids.push_back(b.word_ids.back());
        next_first = b.word_ids.front();
    }

    for (uint64_t i = 0; i < p; i++) prs.num_phrase_runs += prs.phrases_of_thread[i].word_ids.size() - prs.continues[i];

    // sort the words
    dictionary<sym_t>& dict = prs.dict;
    const uint64_t D = table.size();
    dict.num_words = D;
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
    dict.size = 0;
    for (uint64_t s = 0; s < D; s++) {
        dict.size += table.entry(s).len + 1;
        prs.num_phrases += table.entry(s).num_occ;
    }

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

    for (thread_phrases<sym_t>& b : prs.phrases_of_thread) {
        for (uint32_t id : b.repeated_ids) prs.repeated_words.push_back(prs.word_rank[id]);
        b.repeated_ids = std::vector<uint32_t>();
    }

    std::sort(prs.repeated_words.begin(), prs.repeated_words.end());
    prs.repeated_words.erase(std::unique(prs.repeated_words.begin(), prs.repeated_words.end()),
        prs.repeated_words.end());

    return prs;
}

/**
 * @brief copies the phrase runs into arrays, assigns their letters and frees the phrases of prs
 * @param prs the parsing
 * @param threads number of threads
 * @return the arrays
 */
template <typename idx_t, typename sym_t>
phrase_arrays<idx_t, sym_t> build_phrase_arrays(parsing<sym_t>& prs, uint16_t threads)
{
    const uint64_t J = prs.num_phrase_runs;
    const uint64_t w = prs.w;
    const dictionary<sym_t>& dict = prs.dict;
    const uint64_t D = dict.num_words;
    phrase_arrays<idx_t, sym_t> arr;
    arr.P = new_malloc_array<idx_t>(J + 1);
    arr.P[J] = 0;
    // without repeated words, the letter of a run is the rank of its word
    const uint64_t letter_offset = prs.repeated_words.empty() ? 1 : 0;
    /* the values of each thread are stored after those of the thread before; the last value of a thread whose last
       phrase run continues in the next thread is overwritten by the first value of that thread */
    uint64_t offset = 0;

    for (uint64_t t = 0; t < prs.phrases_of_thread.size(); t++) {
        thread_phrases<sym_t>& b = prs.phrases_of_thread[t];
        const uint64_t size = b.word_ids.size();
        b.word_ids.move_to(arr.P.get() + offset, [&](uint32_t id) { return idx_t(letter_offset + prs.word_rank[id]); },
            threads);
        offset += size - prs.continues[t];
    }

    prs.word_rank.reset();
    arr.symbol_before_next = new_malloc_array<sym_t>(J);
    offset = 0;

    for (uint64_t t = 0; t < prs.phrases_of_thread.size(); t++) {
        thread_phrases<sym_t>& b = prs.phrases_of_thread[t];
        const uint64_t size = b.symbol_before_next.size();
        b.symbol_before_next.move_to(arr.symbol_before_next.get() + offset, [](sym_t c) { return c; }, threads);
        offset += size - prs.continues[t];
    }

    const uint8_t bytes = prs.position_bytes;
    arr.ends = uint_array(J, bytes);
    offset = 0;

    for (uint64_t t = 0; t < prs.phrases_of_thread.size(); t++) {
        thread_phrases<sym_t>& b = prs.phrases_of_thread[t];
        const uint64_t size = b.ends.size() / bytes;
        b.ends.move_to(arr.ends.raw() + offset * bytes, [](uint8_t c) { return c; }, threads);
        offset += size - prs.continues[t];
    }

    prs.phrases_of_thread.clear();
    const std::vector<uint64_t>& repeated = prs.repeated_words;
    const uint64_t R = repeated.size();

    if (R == 0) {
        arr.num_letters = D;
        return arr;
    }

    arr.repeated_set = ranked_bits([&](uint64_t y) { return repeated[y]; }, R, D);
    const ranked_bits& is_repeated = arr.repeated_set;

    // returns the number of phrases of run i, whose word is s
    auto copies = [&](uint64_t i, uint64_t s) {
        return (arr.ends[i] - (i == 0 ? w - 1 : arr.ends[i - 1])) / (dict.length(s) - w);
    };

    /* the distinct keys of the phrase runs of the repeated words, (index in repeated, key), and the extra phrases of
       the repeated words, per thread */
    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> keys_of_thread(threads);
    std::vector<std::vector<uint64_t>> extra_of_thread(threads, std::vector<uint64_t>(R, 0));

    #pragma omp parallel num_threads(threads)
    {
        const uint64_t t = omp_get_thread_num();
        std::vector<std::pair<uint64_t, uint64_t>>& keys = keys_of_thread[t];

        #pragma omp for
        for (uint64_t i = 0; i < J; i++) {
            const uint64_t s = uint64_t(arr.P[i]);
            if (!is_repeated.contains(s)) continue;
            const uint64_t y = is_repeated.rank(s);
            const uint64_t m = copies(i, s);
            const bool down = i + 1 == J || uint64_t(arr.P[i + 1]) < s;
            keys.emplace_back(y, run_key(down, m));
            extra_of_thread[t][y] += m - 1;
        }

        std::sort(keys.begin(), keys.end());
        keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
    }

    std::vector<std::pair<uint64_t, uint64_t>> keys;

    for (std::vector<std::pair<uint64_t, uint64_t>>& k : keys_of_thread) {
        keys.insert(keys.end(), k.begin(), k.end());
        k = std::vector<std::pair<uint64_t, uint64_t>>();
    }

    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
    arr.repeated_words.resize(R);
    // [0..R] letters_before[y] = the number of letters of the words in repeated[0..y-1] minus y
    std::vector<uint64_t> letters_before(R + 1, 0);

    for (uint64_t y = 0, z = 0; y < R; y++) {
        repeated_word& rw = arr.repeated_words[y];
        rw.word = repeated[y];
        rw.first_letter = rw.word + letters_before[y];
        for (; z < keys.size() && keys[z].first == y; z++) rw.keys.push_back(keys[z].second);
        for (uint64_t t = 0; t < threads; t++) rw.extra_phrases += extra_of_thread[t][y];
        letters_before[y + 1] = letters_before[y] + rw.keys.size() - 1;
    }

    arr.num_letters = D + letters_before[R];

    // part u replaces the words of the runs [u * J / threads, (u + 1) * J / threads) by their letters; the word after
    // its last run is read first
    std::vector<uint64_t> word_after(threads);
    for (uint64_t u = 0; u < threads; u++) word_after[u] = uint64_t(arr.P[((u + 1) * J) / threads]);

    run_on_threads(threads, [&](uint64_t t, uint64_t team_size) {
        for (uint64_t u = t; u < threads; u += team_size) {
            const uint64_t beg = (u * J) / threads, end = ((u + 1) * J) / threads;

            for (uint64_t i = beg; i < end; i++) {
                const uint64_t s = uint64_t(arr.P[i]);
                uint64_t letter = 0;

                if (!is_repeated.contains(s)) {
                    letter = s + letters_before[is_repeated.rank(s)];
                } else {
                    const repeated_word& rw = arr.repeated_words[is_repeated.rank(s)];
                    const uint64_t after = i + 1 == end ? word_after[u] : uint64_t(arr.P[i + 1]);
                    const bool down = i + 1 == J || after < s;
                    const uint64_t key = run_key(down, copies(i, s));
                    letter = rw.first_letter
                        + (std::lower_bound(rw.keys.begin(), rw.keys.end(), key) - rw.keys.begin());
                }

                arr.P[i] = idx_t(1 + letter);
            }
        }
    });

    return arr;
}

} // namespace pfp_rlbwt

#endif // PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_PARSE_HPP
