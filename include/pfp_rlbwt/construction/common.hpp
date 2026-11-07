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
 * @file common.hpp
 * @brief common parts of the construction: arrays with a fixed byte width per value, block vectors and the output
 */

#ifndef PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_COMMON_HPP
#define PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_COMMON_HPP

#include <algorithm>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <filesystem>
#include <memory>
#include <mutex>
#include <new>
#include <ostream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <type_traits>
#include <utility>
#include <vector>

#include <omp.h>

#include <misc/files.hpp>
#include <misc/memory.hpp>
#include <misc/platform.hpp>
#include <misc/threads.hpp>

#include <pfp_rlbwt/types.hpp>

namespace pfp_rlbwt {

/* the text is parsed as S = $ T $^w, where $ = dollar is smaller than every symbol of T; in the dictionary, each
   word is followed by separator */
inline constexpr uint64_t separator = 0;
inline constexpr uint64_t dollar = 1;

// returns the symbol in the BWT of T$ of a symbol c of S
template <typename sym_t>
inline constexpr sym_t bwt_symbol(sym_t c) { return c == dollar ? sym_t(terminator) : c; }

// frees memory allocated with malloc
struct free_deleter {
    void operator()(void* ptr) const { std::free(ptr); }
};

// an array of values of a trivial type, allocated with malloc (it can be shrunk with realloc)
template <typename val_t>
using malloc_array_t = std::unique_ptr<val_t[], free_deleter>;

// returns an array of n values that are not initialized (on huge pages, see advise_huge_pages)
template <typename val_t>
inline malloc_array_t<val_t> new_malloc_array(uint64_t n)
{
    static_assert(std::is_trivial_v<val_t>);
    val_t* ptr = (val_t*) std::malloc(std::max<uint64_t>(1, n) * sizeof(val_t));
    if (ptr == nullptr) throw std::bad_alloc();
    advise_huge_pages(ptr, n * sizeof(val_t));
    return malloc_array_t<val_t>(ptr);
}

// resizes the array a (see new_malloc_array) to max(1, n) values with realloc; the values that fit are kept
template <typename val_t>
inline void resize_malloc_array(malloc_array_t<val_t>& a, uint64_t n)
{
    val_t* resized = (val_t*) std::realloc(a.release(), std::max<uint64_t>(1, n) * sizeof(val_t));
    if (resized == nullptr) throw std::bad_alloc();
    a.reset(resized);
}

// the fewest and the most bytes per value of a uint_array (32-bit and 48-bit values)
inline constexpr uint8_t min_value_bytes = 4;
inline constexpr uint8_t max_value_bytes = 6;

// returns the number of bytes per value of a uint_array whose values are at most max_value (< 2^48)
inline uint8_t value_bytes(uint64_t max_value)
{
    return std::max<uint8_t>(min_value_bytes, (std::bit_width(max_value) + 7) / 8);
}

/**
 * @brief an array of unsigned values with a fixed number of bytes (min_value_bytes to max_value_bytes) per value,
 *        little endian; a value is read by loading 8 bytes and masking them
 */
class uint_array {
private:
    /* the values: [0..n-1] value i is stored in bytes [i * bytes..(i + 1) * bytes - 1], little endian; followed by
       sizeof(uint64_t) bytes (operator[] loads 8 bytes) */
    malloc_array_t<uint8_t> data;
    uint8_t bytes = 0; // bytes per value
    uint64_t mask = 0; // the low 8 * bytes bits

public:
    uint_array() = default;

    // allocates n values (not initialized) with the given number of bytes each
    uint_array(uint64_t n, uint8_t bytes)
        : uint_array(new_malloc_array<uint8_t>(n * bytes + sizeof(uint64_t)), bytes)
    { }

    // takes the values stored in data (with the given number of bytes each, followed by sizeof(uint64_t) bytes)
    uint_array(malloc_array_t<uint8_t> data, uint8_t bytes)
        : data(std::move(data)), bytes(bytes),
          mask((uint64_t { 1 } << (8 * bytes)) - 1)
    { }

    inline uint64_t operator[](uint64_t i) const
    {
        uint64_t val;
        std::memcpy(&val, data.get() + i * bytes, sizeof(uint64_t));
        return val & mask;
    }

    inline void set(uint64_t i, uint64_t val) { std::memcpy(data.get() + i * bytes, &val, bytes); }

    // fetches the cache line of value i into the cache
    inline void prefetch(uint64_t i) const { __builtin_prefetch(data.get() + i * bytes); }

    // returns the number of bytes per value
    uint8_t width() const { return bytes; }

    // returns the bytes of the values
    uint8_t* raw() { return data.get(); }
};

// number of bytes of the largest blocks of a block_vector
constexpr uint64_t max_block_bytes = uint64_t { 1 } << 18;

/**
 * @brief a sequence of values stored in blocks that are never moved; move_to copies the values into one array and
 *        frees each block after it is copied. The first block has room for min_block_size values, each further block
 *        for twice as many as the one before, up to max_block_size values.
 */
template <typename val_t>
class block_vector {
public:
    static constexpr uint64_t min_block_size = uint64_t { 1 } << 12; // values in the first block
    // values in the largest blocks (max_block_bytes bytes)
    static constexpr uint64_t max_block_size = std::max<uint64_t>(min_block_size, max_block_bytes / sizeof(val_t));

    // reads the values in order without freeing them
    class reader {
    private:
        const block_vector& v;
        uint64_t b = 0; // the block of the next value
        uint64_t i = 0; // the index of the next value in its block

    public:
        // reads from value first on
        explicit reader(const block_vector& v, uint64_t first = 0) : v(v)
        {
            while (b < v.sizes.size() && first >= v.sizes[b]) first -= v.sizes[b++];
            if (b == v.sizes.size() && b != 0) i = v.sizes[--b]; // first = size(): at the end of the last block
            else i = first;
        }

        // returns the next value (there has to be one)
        inline val_t next()
        {
            if (i == v.sizes[b]) {
                b++;
                i = 0;
            }

            return v.blocks[b][i++];
        }
    };

private:
    /* the blocks: [0..blocks.size()-1] blocks[b] stores the values [o_b..o_b + sizes[b] - 1] of the sequence, where
       o_b = sizes[0] + ... + sizes[b - 1]; block b has room for min(min_block_size * 2^b, max_block_size) values */
    std::vector<malloc_array_t<val_t>> blocks;
    // [0..blocks.size()-1] sizes[b] = number of values in block b
    std::vector<uint64_t> sizes;
    uint64_t fill = 0; // number of values in the last block
    uint64_t capacity = 0; // size of the last block

public:
    inline void push_back(val_t val)
    {
        if (fill == capacity) [[unlikely]] {
            capacity = blocks.empty() ? min_block_size : std::min(2 * capacity, max_block_size);
            blocks.emplace_back(new_malloc_array<val_t>(capacity));
            sizes.push_back(0);
            fill = 0;
        }

        blocks.back()[fill++] = val;
        sizes.back() = fill;
    }

    uint64_t size() const
    {
        uint64_t s = 0;
        for (uint64_t b : sizes) s += b;
        return s;
    }

    // frees the unused part of the last block
    void shrink_to_fit()
    {
        if (blocks.empty() || fill == capacity) return;
        resize_malloc_array(blocks.back(), fill);
        capacity = fill;
    }

    // returns the first value (the sequence must not be empty)
    val_t front() const { return blocks.front()[0]; }

    // returns the last value (the sequence must not be empty)
    val_t back() const { return blocks.back()[fill - 1]; }

    /**
     * @brief copies the values to dst (in order), applying fnc to each, and frees the blocks
     * @param dst destination of size() values
     * @param fnc function applied to each value
     * @param threads number of threads
     */
    template <typename dst_t, typename fnc_t>
    void move_to(dst_t* dst, fnc_t fnc, uint16_t threads)
    {
        move_to([&](uint64_t i, const val_t& val) { dst[i] = fnc(val); }, threads);
    }

    /**
     * @brief calls fnc(i, value) for each value i and frees the blocks
     * @param fnc function called for each value
     * @param threads number of threads
     */
    template <typename fnc_t>
    void move_to(fnc_t fnc, uint16_t threads)
    {
        move_blocks([&](uint64_t offset, const val_t* values, uint64_t num) {
            for (uint64_t i = 0; i < num; i++) fnc(offset + i, values[i]);
        }, threads);
    }

    /**
     * @brief calls fnc(offset, values, num) for each block, where values[0..num-1] are the values offset, ...,
     *        offset + num - 1, and frees the blocks; the blocks are processed in parallel and in any order
     * @param fnc function called for each block
     * @param threads number of threads
     */
    template <typename fnc_t>
    void move_blocks(fnc_t fnc, uint16_t threads)
    {
        // offsets[b] = number of values in the blocks before block b
        std::vector<uint64_t> offsets(blocks.size() + 1, 0);
        for (uint64_t b = 0; b < blocks.size(); b++) offsets[b + 1] = offsets[b] + sizes[b];
        std::vector<std::exception_ptr> failures(blocks.size());

        #pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
        for (uint64_t b = 0; b < blocks.size(); b++) {
            try {
                fnc(offsets[b], blocks[b].get(), sizes[b]);
            } catch (...) {
                failures[b] = std::current_exception();
            }

            blocks[b].reset();
        }

        blocks.clear();
        sizes.clear();
        fill = capacity = 0;
        rethrow_first_exception(failures);
    }
};

// number of bytes of the full BWT, SA and LCP array a chunk collects before they are written to their file, and of
// the records the runs file is written in
constexpr uint64_t write_buffer_bytes = uint64_t { 1 } << 20;

/**
 * @brief writes values given with a fixed number of bytes each (little endian) into a uint_array: the value written at
 *        the byte offset o is stored at index o / value_bytes
 */
class array_output : public positioned_output {
private:
    uint_array& values; // the destination
    uint8_t value_bytes; // bytes per written value

public:
    array_output(uint_array& values, uint8_t value_bytes) : values(values), value_bytes(value_bytes) { }

    void write(const void* data, uint64_t size, uint64_t offset) override
    {
        const uint8_t* bytes = (const uint8_t*) data;
        const uint64_t first = offset / value_bytes;

        for (uint64_t i = 0; i < size / value_bytes; i++) {
            uint64_t val = 0;
            std::memcpy(&val, bytes + i * value_bytes, value_bytes);
            values.set(first + i, val);
        }
    }

    void close() override { }
};

/**
 * @brief writes the values given with a fixed number of bytes each (little endian) at the positions of given ranges
 *        into a uint_array: the values of the positions [beg, end) of a range follow those of the ranges before it;
 *        the values at other positions are skipped
 */
class ranges_output : public positioned_output {
private:
    uint_array& values; // the destination
    uint8_t value_bytes; // bytes per written value
    std::vector<std::pair<uint64_t, uint64_t>> ranges; // the ranges [beg, end), sorted and disjoint
    std::vector<uint64_t> offsets; // [0..ranges.size()-1] offsets[i] = the number of positions of the ranges before i

public:
    ranges_output(uint_array& values, uint8_t value_bytes, std::vector<std::pair<uint64_t, uint64_t>> ranges)
        : values(values), value_bytes(value_bytes), ranges(std::move(ranges)), offsets(this->ranges.size(), 0)
    {
        for (uint64_t i = 1; i < this->ranges.size(); i++) {
            offsets[i] = offsets[i - 1] + this->ranges[i - 1].second - this->ranges[i - 1].first;
        }
    }

    void write(const void* data, uint64_t size, uint64_t offset) override
    {
        const uint8_t* bytes = (const uint8_t*) data;
        const uint64_t first = offset / value_bytes, last = first + size / value_bytes;
        // the first range that ends after first
        uint64_t i = std::upper_bound(ranges.begin(), ranges.end(), first,
            [](uint64_t pos, const std::pair<uint64_t, uint64_t>& r) { return pos < r.second; }) - ranges.begin();

        for (; i < ranges.size() && ranges[i].first < last; i++) {
            for (uint64_t pos = std::max(first, ranges[i].first); pos < std::min(last, ranges[i].second); pos++) {
                uint64_t val = 0;
                std::memcpy(&val, bytes + (pos - first) * value_bytes, value_bytes);
                values.set(offsets[i] + pos - ranges[i].first, val);
            }
        }
    }

    void close() override { }
};

/* the level of the recursive prefix-free parsing that the calling thread computes (1: the text, k >= 2: the parse of
   level k - 1, see build_level); it selects the names in the log (see level_name), the number of symbols read at once
   and where the temporary arrays are kept. It is 1 in the threads that a construction starts (read it before). */
inline thread_local uint64_t current_level = 1;

// returns the name of a structure of the current level in the log: name on level 1, name_k on level k >= 2
inline std::string level_name(const std::string& name)
{
    return current_level == 1 ? name : name + "_" + std::to_string(current_level);
}

// the path of a temporary file, which is removed when the object is destroyed
class temporary_file {
private:
    std::string file_path; // the path (empty: no file)

public:
    temporary_file() = default;
    explicit temporary_file(std::string path) : file_path(std::move(path)) { }
    ~temporary_file() { remove(); }

    temporary_file(const temporary_file&) = delete;
    temporary_file& operator=(const temporary_file&) = delete;
    temporary_file(temporary_file&& other) noexcept : file_path(std::exchange(other.file_path, std::string())) { }

    temporary_file& operator=(temporary_file&& other) noexcept
    {
        if (this != &other) {
            remove();
            file_path = std::exchange(other.file_path, std::string());
        }

        return *this;
    }

    const std::string& path() const { return file_path; }

    // removes the file
    void remove()
    {
        if (file_path.empty()) return;
        std::error_code ec;
        std::filesystem::remove(file_path, ec);
        file_path.clear();
    }
};

/**
 * @brief writes bytes to a new file (see direct_ofstream)
 * @param path path of the file
 * @param data the bytes
 * @param size number of bytes
 * @return the file
 */
inline temporary_file write_temporary_file(const std::string& path, const void* data, uint64_t size)
{
    temporary_file file(path);
    direct_ofstream out(path);
    if (!out.is_open()) throw std::runtime_error("pfp_rlbwt: cannot write " + path);
    write_in_chunks(out, (const char*) data, size);
    out.close();
    if (out.fail()) throw std::runtime_error("pfp_rlbwt: cannot write " + path);
    return file;
}

// number of bytes a reader of a temporary file reads at once (see array_file_reader)
constexpr uint64_t temporary_read_bytes = uint64_t { 1 } << 20;

// writes bytes at given offsets into an array in RAM (several threads, in any order)
class memory_output : public positioned_output {
private:
    uint8_t* bytes; // the destination

public:
    explicit memory_output(uint8_t* bytes) : bytes(bytes) { }
    void write(const void* data, uint64_t size, uint64_t offset) override { std::memcpy(bytes + offset, data, size); }
    void close() override { }
};

/**
 * @brief returns whether an array the construction keeps temporarily is kept in RAM instead of in a file: on the levels
 *        after the first of the recursive prefix-free parsing, whose arrays are small (see current_level), and if it
 *        fits into the buffer of a file written with O_DIRECT (then storing it does not lower the peak memory usage)
 * @param bytes the size of the array in bytes
 */
inline bool temporary_arrays_in_memory(uint64_t bytes)
{
    return current_level >= 2 || bytes < direct_filebuf_base::buffer_length;
}

/**
 * @brief an array that the construction keeps temporarily: in a file (see temporary_file) or, on the levels after the
 *        first of the recursive prefix-free parsing, whose arrays are small, in RAM (followed by sizeof(uint64_t)
 *        bytes, as a uint_array)
 */
class temporary_array {
private:
    std::string file_path; // path of the file (empty if the array is stored in RAM)
    temporary_file file; // the file, once it is written (removed with the array)
    malloc_array_t<uint8_t> memory; // the bytes, if the array is stored in RAM

public:
    temporary_array() = default;

    /**
     * @brief an empty array; exactly one of output, write and store fills it
     * @param path path of its file (not created if the array is kept in RAM)
     * @param in_memory true <=> the array is kept in RAM
     */
    temporary_array(const std::string& path, bool in_memory) : file_path(in_memory ? std::string() : path) { }

    // returns whether the array is stored in RAM
    bool is_in_memory() const { return file_path.empty(); }

    // returns the path of the file (empty if the array is stored in RAM)
    const std::string& path() const { return file_path; }

    // returns the bytes of the array (only if it is stored in RAM)
    const uint8_t* data() const { return memory.get(); }

    // returns a destination of size bytes that several threads write to; it has to be closed before the array is read
    std::unique_ptr<positioned_output> output(uint64_t size)
    {
        if (is_in_memory()) {
            memory = new_malloc_array<uint8_t>(size + sizeof(uint64_t));
            return std::make_unique<memory_output>(memory.get());
        }

        file = temporary_file(file_path);
        return std::make_unique<parallel_output_file>(file_path);
    }

    // stores a copy of the bytes data[0..size-1]
    void write(const void* data, uint64_t size)
    {
        if (is_in_memory()) {
            memory = new_malloc_array<uint8_t>(size + sizeof(uint64_t));
            std::memcpy(memory.get(), data, size);
        } else {
            file = write_temporary_file(file_path, data, size);
        }
    }

    /**
     * @brief stores the values of an array allocated with malloc and frees it (in RAM, its memory is taken over)
     * @param a the array a[0..n-1]
     * @param n number of values
     */
    template <typename val_t>
    void store(malloc_array_t<val_t> a, uint64_t n)
    {
        const uint64_t size = n * sizeof(val_t);

        if (!is_in_memory()) {
            file = write_temporary_file(file_path, a.get(), size);
            return;
        }

        memory.reset((uint8_t*) a.release());
        resize_malloc_array(memory, size + sizeof(uint64_t));
    }

    // returns the bytes of an array stored in RAM, followed by sizeof(uint64_t) bytes, and empties the array
    malloc_array_t<uint8_t> release() { return std::move(memory); }

    // removes the file or frees the bytes
    void remove()
    {
        file.remove();
        memory.reset();
    }
};

/**
 * @brief reads the values of a temporary_array (a fixed number of bytes per value, little endian): from its file with
 *        an array_file_reader, or from RAM
 */
class array_reader {
private:
    std::unique_ptr<array_file_reader> file; // the reader of the file (nullptr if the array is in RAM)
    const uint8_t* memory; // the bytes of the array, if it is in RAM
    uint64_t width; // bytes per value
    uint64_t mask; // the low 8 * width bits

public:
    /**
     * @brief opens an array
     * @param a the array
     * @param width bytes per value
     * @param read_bytes number of bytes read from the file at once (at least width)
     */
    array_reader(const temporary_array& a, uint64_t width, uint64_t read_bytes = temporary_read_bytes)
        : memory(a.data()), width(width),
          mask(width >= sizeof(uint64_t) ? ~uint64_t { 0 } : (uint64_t { 1 } << (8 * width)) - 1)
    {
        if (!a.is_in_memory()) file = std::make_unique<array_file_reader>(a.path(), width, read_bytes);
    }

    // returns the bytes of value i
    inline const uint8_t* bytes(uint64_t i) { return file == nullptr ? memory + i * width : file->bytes(i); }

    // returns value i (width <= 8)
    inline uint64_t operator[](uint64_t i)
    {
        uint64_t val;
        std::memcpy(&val, bytes(i), sizeof(uint64_t));
        return val & mask;
    }

    // copies the values [from, to) to dst
    void read(uint64_t from, uint64_t to, void* dst)
    {
        if (file != nullptr) file->read(from, to, dst);
        else std::memcpy(dst, memory + from * width, (to - from) * width);
    }
};

} // namespace pfp_rlbwt

#endif // PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_COMMON_HPP
