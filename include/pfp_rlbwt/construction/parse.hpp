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
 * @file parse.hpp
 * @brief the prefix-free parsing: the dictionary, the word ids and the phrase runs
 */

#ifndef PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_PARSE_HPP
#define PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_PARSE_HPP

#include <algorithm>
#include <array>
#include <atomic>
#include <cstring>
#include <exception>
#include <functional>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <omp.h>

#include <misc/hashing.hpp>
#include <misc/log.hpp>
#include <misc/memory.hpp>
#include <misc/sort.hpp>
#include <pfp_rlbwt/construction/common.hpp>
#include <pfp_rlbwt/construction/phrase_runs.hpp>
#include <pfp_rlbwt/text_ranges.hpp>
#include <pfp_rlbwt/types.hpp>

namespace pfp_rlbwt {

/* The prefix-free parsing of S = $ T $^w: a window S[b..b+w-1] with 1 <= b and b+w-1 <= n is a trigger string if
   its hash is below a threshold (see rolling_window_hash); with b_1 < ... < b_{K-1} the starts of the trigger strings,
   b_0 = 0 and b_K = n+1, phrase k is S[b_k..e_k] with e_k = b_{k+1}+w-1. Consecutive phrases overlap by w symbols.
   Each position x in [0..n] of S lies in exactly one range [b_k..e_k-w]; the suffix of S starting at x corresponds
   to the suffix of T$ starting at x-1 (x >= 1) and to the suffix $ (x = 0). A phrase run is a maximal run of equal
   consecutive phrases; the parsing stores the J phrase runs by their words and numbers of phrases. Phrase run i
   consists of m_i phrases equal to a word of length w + a, where a is the advance of the word, and ends at
   F_i = F_{i-1} + m_i * a (F_{-1} = w-1). P, the parse of the phrase runs, has one symbol per phrase run (see
   phrase_arrays::P). */

// the dictionary: the distinct phrases (words) in lexicographic order
template <typename sym_t>
struct dictionary {
    uint64_t num_words = 0; // D, the number of words
    uint64_t size = 0; // total length of the words and their separators
    // the words: [0..size-1] the words in lexicographic order, each followed by the separator
    malloc_array_t<sym_t> text;
    // the starting positions of the words: [0..D] word_start[s] = position of word s in text; word_start[D] = size
    uint_array word_start;
    // the numbers of occurrences of the words: [0..D-1] num_occ[s] = number of phrases equal to word s
    uint_array num_occ;

    // returns the length of word s
    inline uint64_t length(uint64_t s) const { return word_start[s + 1] - word_start[s] - 1; }
};

// appends x as a varint to a sequence of bytes: 7 bits per byte, the low bits first, the high bit set in each byte but
// the last
template <typename seq_t>
inline void append_varint(seq_t& out, uint64_t x)
{
    for (; x >= 128; x >>= 7) out.push_back(uint8_t((x & 127) | 128));
    out.push_back(uint8_t(x));
}

// returns the next varint of a reader of a sequence of bytes (see append_varint)
template <typename reader_t>
inline uint64_t read_varint(reader_t& in);

/**
 * @brief the word ids of the phrase runs a thread has found (see word_table), in blocks that grow as in block_vector; a
 *        block stores 3 bytes per id while the ids are below 2^24; with a larger id, a block with 4 bytes per id is
 *        started, and all blocks after it also have 4 bytes per id
 */
class word_id_vector {
private:
    // bytes per id below 2^(8 short_bytes), and per larger id
    static constexpr uint8_t short_bytes = 3, long_bytes = sizeof(uint32_t);
    static constexpr uint64_t short_mask = (uint64_t { 1 } << (8 * short_bytes)) - 1;

    struct block {
        // the ids, width bytes each, followed by long_bytes - width bytes (get reads long_bytes bytes per id)
        malloc_array_t<uint8_t> bytes;
        uint64_t size; // the number of ids
        uint8_t width; // bytes per id
    };

public:
    // reads the ids in order without freeing them
    class reader {
    private:
        const word_id_vector& v;
        uint64_t b = 0; // the block of the next id
        uint64_t i = 0; // the index of the next id in its block

    public:
        // reads from id first on
        explicit reader(const word_id_vector& v, uint64_t first = 0) : v(v)
        {
            while (b < v.blocks.size() && first >= v.blocks[b].size) first -= v.blocks[b++].size;
            if (b == v.blocks.size() && b != 0) i = v.blocks[--b].size; // first = size(): at the end of the last block
            else i = first;
        }

        // returns the next id (there has to be one)
        inline uint32_t next();
    };

private:
    /* the blocks: [0..blocks.size()-1] blocks[b] stores the ids after those of the blocks before it; block b has room
       for min(2^b block_vector<uint32_t>::min_block_size, max_block_bytes / width) ids */
    std::vector<block> blocks;
    uint64_t capacity = 0; // the number of ids the last block has room for

    // returns id i of block b
    static inline uint32_t get(const block& b, uint64_t i)
    {
        uint32_t id;
        std::memcpy(&id, b.bytes.get() + i * b.width, sizeof(uint32_t));
        return b.width == long_bytes ? id : uint32_t(id & short_mask);
    }

    // appends a block with room for more ids than the last one, with 4 bytes per id if id or the last block needs them
    void add_block(uint32_t id);

public:
    inline void push_back(uint32_t id);

    uint64_t size() const
    {
        uint64_t s = 0;
        for (const block& b : blocks) s += b.size;
        return s;
    }

    // frees the unused part of the last block
    void shrink_to_fit()
    {
        if (blocks.empty() || blocks.back().size == capacity) return;
        block& b = blocks.back();
        resize_malloc_array(b.bytes, b.size * b.width + long_bytes - b.width);
        capacity = b.size;
    }

    // returns the first id (the sequence must not be empty)
    uint32_t front() const { return get(blocks.front(), 0); }

    // returns the last id (the sequence must not be empty)
    uint32_t back() const { return get(blocks.back(), blocks.back().size - 1); }

    /**
     * @brief calls fnc(i, id) for the i-th id of the sequence, for each i, and frees the blocks; the blocks are
     *        processed in parallel and in any order
     * @param fnc function called for each id
     * @param threads number of threads
     */
    template <typename fnc_t>
    void move_to(fnc_t fnc, uint16_t threads);
};

/* the phrase runs one thread has found, in text order (the j-th of them consists of the phrases k..k'); aligned to a
   cache line (the threads append to different ones) */
struct alignas(cache_line_bytes) thread_phrases {
    // [0..num_runs-1] word_ids[j] = id of the word equal to the phrases k..k' of the j-th phrase run (see word_table)
    word_id_vector word_ids;
    /* bit j % 64 of multi_runs[j / 64] = 1 <=> the j-th phrase run has two phrases or more; the bits of the phrase runs
       after the last full word are in multi_bits */
    block_vector<uint64_t> multi_runs;
    uint64_t multi_bits = 0;
    // the numbers of phrases k' - k + 1 of the phrase runs with two phrases or more, one varint each (see read_varint)
    block_vector<uint8_t> copies;
    uint64_t num_runs = 0; // the number of phrase runs
    uint64_t length = 0; // the total length of the phrase runs (their numbers of phrases times their advances)
    uint64_t last_copies = 0; // the number of phrases of the last phrase run
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
    std::vector<thread_phrases> phrases_of_thread;
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
    /* the parse of the phrase runs (if it is not suffix-sorted recursively): [0..J] P[i] = 1 + letter of phrase run i;
       P[J] = 0. The letters are ordered by word; the phrase runs of a word in repeated_words are assigned one letter
       per key (see run_key), the others one letter per word */
    malloc_array_t<idx_t> P;
    uint64_t num_letters = 0; // the number of letters
    uint8_t position_bytes = 0; // bytes per end of a phrase run
    /* the ends of the phrase runs in text order (in a file, see temporary_array), until they are ordered by rank:
       value i in [0..J-1] of ends_file = F_i, the end of the last phrase of run i (position_bytes bytes, little
       endian) */
    temporary_array ends_file;
    sym_t last_symbol = 0; // the symbol before the last w symbols of the word of the last run
    /* the ends in the order of the suffix array SA of the parse of the phrase runs: [0..J] ends_by_rank[x] =
       F_{SA[x]-1}, the end of the run before run SA[x] (F_{-1} = w-1 for SA[x] = 0) */
    uint_array ends_by_rank;
    /* the BWT symbols of the values of P: [0..num_letters] value_symbols[v] = the symbol before the last w symbols of
       the word of value v (see symbols_of_values) */
    std::vector<sym_t> value_symbols;
    // the R words with a phrase run of at least two phrases, in increasing order
    std::vector<repeated_word> repeated_words;
    ranked_bits repeated_set; // the set of these words
    // [0..R] letters_before[y] = the number of letters of the words repeated_words[0..y-1] minus y
    std::vector<uint64_t> letters_before;
    // (i, letter of phrase run i) for the phrase runs i of the repeated words, in increasing order of i
    std::vector<std::pair<uint64_t, uint64_t>> repeated_run_letters;
    /* if the parse is suffix-sorted recursively, instead of P: the phrase runs as the parsing threads found them (only
       their word ids, see thread_phrases), the ranks of their words (see parsing::word_rank) and, for each part,
       [0..parts] run_offsets[u] = the number of phrase runs the parts before part u store (see build_phrase_arrays) */
    std::vector<thread_phrases> runs_of_part;
    malloc_array_t<uint32_t> word_rank;
    std::vector<uint64_t> run_offsets;

    /* the letters per word id: [0..D-1] id_letter[id] = the letter of the phrase runs of the word with id `id`, or
       num_letters if the word is repeated (see repeated_run_letters) */
    uint_array id_letter;

    // returns the letter of phrase run i, whose word has the id `id`
    inline uint64_t letter_of(uint64_t i, uint64_t id) const
    {
        const uint64_t letter = id_letter[id];
        if (letter != num_letters) [[likely]] return letter;
        return std::lower_bound(repeated_run_letters.begin(), repeated_run_letters.end(),
            std::make_pair(i, uint64_t { 0 }))->second;
    }

    // frees the word ids of the phrase runs, the ranks of the words and the letters per word id (see runs_of_part)
    void free_word_ids()
    {
        runs_of_part = std::vector<thread_phrases>();
        word_rank.reset();
        id_letter = uint_array();
    }
};

/* returns the key of a phrase run with m phrases that is followed by a smaller (down) or larger suffix: the keys of
   the down runs increase with m, those of the up runs decrease with m and are larger */
inline uint64_t run_key(bool down, uint64_t m) { return down ? m : std::numeric_limits<uint64_t>::max() - m; }

// returns the number of phrases m of a phrase run with the key run_key(down, m) (m < 2^63)
inline uint64_t run_copies(uint64_t key)
{
    return key < (uint64_t { 1 } << 63) ? key : std::numeric_limits<uint64_t>::max() - key;
}

/* a parse with fewer than n / (degenerate_parse_tolerance * modulus) phrases is degenerate (it is expected to have
   about n / modulus phrases); the seed of the trigger strings (see rolling_window_hash) is then switched between 0
   and retry_trigger_seed: before the parse on a sample (see choose_trigger_seed), and after it, if the parse is
   degenerate */
constexpr uint64_t degenerate_parse_tolerance = 2;
constexpr uint64_t retry_trigger_seed = golden_ratio_64;

/* number of bytes a parsing thread reads from the text at once; on the levels after the first of the recursive
   prefix-free parsing, whose texts are in RAM, level_block_bytes */
constexpr uint64_t parse_block_bytes = uint64_t { 1 } << 20;
constexpr uint64_t level_block_bytes = uint64_t { 1 } << 16;

// number of symbols of type sym_t a parsing thread reads from the text at once
template <typename sym_t>
inline uint64_t block_symbols()
{
    return (current_level == 1 ? parse_block_bytes : level_block_bytes) / sizeof(sym_t);
}

// number of bytes of the ends of phrase runs a thread collects before it writes them
constexpr uint64_t run_buffer_bytes = uint64_t { 1 } << 18;

// number of bytes of the first and of the largest blocks of a word_arena
constexpr uint64_t min_arena_block_size = uint64_t { 1 } << 16;
constexpr uint64_t max_arena_block_size = uint64_t { 1 } << 20;

// base-2 logarithm of the number of word entries per block of word_table::blocks
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
inline uint64_t string_hash(const uint8_t* s, uint64_t len);

/**
 * @brief the polynomial hash (mod 2^64, base golden_ratio_64) of the last w symbols of S, updated symbol by symbol;
 *        the window is a trigger string if mix(hash xor seed) is smaller than UINT64_MAX / modulus
 */
template <typename sym_t>
class rolling_window_hash {
private:
    uint64_t threshold; // UINT64_MAX / modulus
    uint64_t pow_w = 1; // golden_ratio_64^w
    uint64_t hash = 0; // the hash of the window
    uint64_t seed; // the seed on which the trigger strings depend (see is_trigger)

public:
    rolling_window_hash(uint64_t w, uint64_t modulus, uint64_t seed) : threshold(UINT64_MAX / modulus), seed(seed)
    {
        for (uint64_t i = 0; i < w; i++) pow_w *= golden_ratio_64;
    }

    // appends c to the window; out is the symbol w positions before c, which is removed from it (0 if there is none)
    inline void add(sym_t c, sym_t out) { hash = hash * golden_ratio_64 + uint64_t(c) - uint64_t(out) * pow_w; }

    // returns whether the window is a trigger string
    inline bool is_trigger() const { return mix(hash ^ seed) < threshold; }
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
    const void* store(const void* s, uint64_t len);
};

// a word and its number of occurrences
template <typename sym_t>
struct word_entry {
    const sym_t* str; // the symbols of the word (in a word_arena)
    uint64_t len; // its length
    uint64_t num_occ; // its number of occurrences
};

/**
 * @brief the words found by the parsing threads: each distinct word is assigned an id (0, 1, ... in the order in
 *        which the words are found) and its occurrences are counted. A word is looked up by its hash in one of the
 *        shards, each a hash table (open addressing, at most half full) with its own mutex.
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

    // a slot of a shard
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
    entry_t* block_of(uint64_t id);

    // doubles the number of slots of a shard
    static void grow(shard& sh);

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
    uint32_t count_occurrences(const sym_t* s, uint64_t len, uint64_t hash, uint64_t count, word_arena& arena);
};

/**
 * @brief parses the phrases of S starting in a range of T (see text_ranges) and appends their phrase runs to out;
 *        positions are counted from the range on: its first symbol is at position 1 (in the first range, as in S)
 * @param in reads T from the start of the range on
 * @param first true <=> the range is the first one; the phrase starting at 0 is parsed too
 * @param prm parameters (window and modulus)
 * @param seed the seed of the trigger strings (see rolling_window_hash)
 * @param table the word table
 * @param arena the arena of the calling thread
 * @param out the phrases of the calling thread
 * @param block_len number of symbols read at once (see block_symbols)
 * @return the largest symbol read
 */
template <typename sym_t>
uint64_t parse_phrases_in_range(typename text_ranges<sym_t>::reader& in, bool first, const params& prm,
    uint64_t seed, word_table<sym_t>& table, word_arena& arena, thread_phrases& out, uint64_t block_len);

/**
 * @brief chooses the seed of the trigger strings (see rolling_window_hash) from a sample of T, the first
 *        block_symbols<sym_t>() symbols of each of the p ranges the text is parsed in: 0, unless the sample has fewer
 *        than 1 / (degenerate_parse_tolerance * modulus) trigger strings per window with it and more with
 *        retry_trigger_seed
 * @param text the text T
 * @param prm parameters (window, modulus, num_threads)
 * @return the seed
 */
template <typename sym_t>
uint64_t choose_trigger_seed(text_ranges<sym_t>& text, const params& prm);

/**
 * @brief returns the frequencies of the symbols of T (of one byte) from the words of the parsing: each position of $ T
 *        lies in the first len - w symbols of exactly one phrase of length len (see parse; the dollar is not counted)
 * @param table the word table
 * @param w the length of a trigger string
 * @param p number of threads
 */
template <typename sym_t>
std::array<uint64_t, 256> count_symbols(word_table<sym_t>& table, uint64_t w, uint16_t p);

/**
 * @brief renames the symbols of T in the words of the parsing (the dollars are not renamed)
 * @param table the word table
 * @param map [0..255] map[c] = the new name of the symbol c (at least min_symbol)
 * @param max_symbol (output) the largest symbol of the renamed words
 * @param p number of threads
 */
template <typename sym_t>
void rename_symbols(word_table<sym_t>& table, const std::array<uint8_t, 256>& map, uint64_t& max_symbol, uint16_t p);

/**
 * @brief computes the prefix-free parsing of S = $ T $^w and its dictionary
 * @param text the text T in prm.num_threads ranges (see text_ranges::split); if text_ranges::ranges_read returns a
 *        map, the symbols of T are renamed by it in the words before they are sorted
 * @param prm parameters
 * @param seed the seed of the trigger strings (see rolling_window_hash)
 * @param parsed if not nullptr, called with the sizes of the parsing (see result) once the text is not read any more,
 *        before the dictionary is built; if it returns false, the parsing is returned without its dictionary
 * @return the parsing
 */
template <typename sym_t>
parsing<sym_t> parse(text_ranges<sym_t>& text, const params& prm, uint64_t seed,
    const std::function<bool(const result&)>* parsed = nullptr);

/**
 * @brief stores the letters of the phrase runs in P, the parse of the phrase runs, and frees the word ids of the phrase
 *        runs (see phrase_arrays::free_word_ids)
 * @param arr the phrase runs
 * @param J number of phrase runs
 * @param threads number of threads
 */
template <typename idx_t, typename sym_t>
void store_parse(phrase_arrays<idx_t, sym_t>& arr, uint64_t J, uint16_t threads);

/**
 * @brief writes the ends of the phrase runs to a file (see phrase_arrays::ends_file) and assigns the letters of the
 *        phrase runs; without recursion, stores the letters in P and frees the word ids of the phrase runs, else keeps
 *        the word ids, from which the letters are read on the next level (see letter_source)
 * @param prs the parsing
 * @param prefix path prefix of the file (see ext_phrase_run_ends)
 * @param threads number of threads
 * @param recursive true <=> P is suffix-sorted recursively (see sort_parse_recursively); then P is not built
 * @return the arrays
 */
template <typename idx_t, typename sym_t>
phrase_arrays<idx_t, sym_t> build_phrase_arrays(parsing<sym_t>& prs, const std::string& prefix, uint16_t threads,
    bool recursive);

/**
 * @brief the text T' = (P[0] + 1) (P[1] + 1) ... (P[J-1] + 1) of the next level of the recursive prefix-free parsing
 *        (see sort_parse_recursively), read from the word ids of the phrase runs (see phrase_arrays::runs_of_part)
 * @tparam text_t type of a symbol of T'
 */
template <typename idx_t, typename sym_t, typename text_t>
class letter_source : public text_source<text_t> {
private:
    // reads T' from the word ids of the phrase runs
    class letter_reader : public text_source<text_t>::reader {
    private:
        const phrase_arrays<idx_t, sym_t>& arr; // the phrase runs

    public:
        explicit letter_reader(const phrase_arrays<idx_t, sym_t>& arr) : arr(arr) { }

        void read(uint64_t pos, uint64_t len, text_t* dst) override;
    };

    const phrase_arrays<idx_t, sym_t>& arr; // the phrase runs
    uint64_t J; // the number of phrase runs

public:
    letter_source(const phrase_arrays<idx_t, sym_t>& arr, uint64_t J) : arr(arr), J(J) { }

    uint64_t size() const override { return J; }

    std::unique_ptr<typename text_source<text_t>::reader> new_reader() const override
    {
        return std::make_unique<letter_reader>(arr);
    }
};

} // namespace pfp_rlbwt

#include <pfp_rlbwt/construction/parse.tpp>

#endif // PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_PARSE_HPP
