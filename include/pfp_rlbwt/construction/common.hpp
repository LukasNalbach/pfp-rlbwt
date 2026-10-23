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

// returns an array of n values that are not initialized
template <typename val_t>
inline malloc_array_t<val_t> new_malloc_array(uint64_t n)
{
    static_assert(std::is_trivial_v<val_t>);
    val_t* ptr = (val_t*) std::malloc(std::max<uint64_t>(1, n) * sizeof(val_t));
    if (ptr == nullptr) throw std::bad_alloc();
    return malloc_array_t<val_t>(ptr);
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

    /**
     * @brief stores the values of an array with as few bytes per value as hold max_value, in the memory of the array
     * @param a the array a[0..n-1]; its memory is taken over and shrunk
     * @param n number of values
     * @param max_value the largest value
     */
    template <typename val_t>
    static uint_array compact(malloc_array_t<val_t> a, uint64_t n, uint64_t max_value)
    {
        const uint8_t bytes = value_bytes(max_value);
        uint8_t* raw = (uint8_t*) a.release();

        // value i moves from offset i * sizeof(val_t) to offset i * bytes <= i * sizeof(val_t), so ascending i is safe
        if (bytes != sizeof(val_t)) {
            for (uint64_t i = 0; i < n; i++) {
                val_t val;
                std::memcpy(&val, raw + i * sizeof(val_t), sizeof(val_t));
                std::memcpy(raw + i * bytes, &val, bytes);
            }
        }

        uint8_t* resized = (uint8_t*) std::realloc(raw, n * bytes + sizeof(uint64_t));
        if (resized == nullptr) throw std::bad_alloc();
        return uint_array(malloc_array_t<uint8_t>(resized), bytes);
    }

    inline uint64_t operator[](uint64_t i) const
    {
        uint64_t val;
        std::memcpy(&val, data.get() + i * bytes, sizeof(uint64_t));
        return val & mask;
    }

    inline void set(uint64_t i, uint64_t val) { std::memcpy(data.get() + i * bytes, &val, bytes); }

    // returns the number of bytes per value
    uint8_t width() const { return bytes; }

    // returns the bytes of the values
    uint8_t* raw() { return data.get(); }
};

/**
 * @brief an array of indices computed by libsais: int32_t values as they are, wider values with as few bytes each
 *        as hold the largest (see uint_array)
 */
template <typename idx_t>
class index_array {
private:
    /* the values: [0..n-1]; for int32_t the array computed by libsais, else that array compacted to the fewest bytes
       per value that hold the largest (see uint_array::compact) */
    std::conditional_t<sizeof(idx_t) == sizeof(int32_t), malloc_array_t<idx_t>, uint_array> values;

public:
    index_array() = default;

    /**
     * @brief takes an array computed by libsais
     * @param a the array a[0..n-1]
     * @param n number of values
     * @param max_value the largest value
     */
    index_array(malloc_array_t<idx_t> a, uint64_t n, uint64_t max_value)
    {
        if constexpr (sizeof(idx_t) == sizeof(int32_t)) values = std::move(a);
        else values = uint_array::compact(std::move(a), n, max_value);
    }

    inline uint64_t operator[](uint64_t i) const { return uint64_t(values[i]); }
};

/**
 * @brief a sequence of values stored in blocks that are never moved; move_to copies the values into one array and
 *        frees each block after it is copied. The first block has room for min_block_size values, each further block
 *        for twice as many as the one before, up to max_block_size values.
 */
template <typename val_t>
class block_vector {
public:
    static constexpr uint64_t min_block_size = uint64_t { 1 } << 12; // values in the first block
    static constexpr uint64_t max_block_size = uint64_t { 1 } << 20; // values in the largest blocks

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
 * @brief a file that several threads write to at given offsets, in any order, with O_DIRECT (if the file system does
 *        not support it, through the page cache). The sectors of a write are written directly; the bytes of a sector
 *        that the write covers only partly are collected until the sector is complete. close writes the remaining
 *        sectors padded with zeros and truncates the file to its length.
 */
class parallel_output_file {
private:
    static constexpr uint64_t sector = direct_io_alignment;

    // a sector that has been written only partly
    struct partial_sector {
        uint8_t* bytes; // the sector (aligned to sector), zero where it has not been written
        uint64_t written; // the number of bytes written
    };

    std::string path; // path of the file
    int fd = -1; // file descriptor
    bool direct = false; // true <=> the file is opened with O_DIRECT
    uint64_t length = 0; // the length of the file: the largest offset after a written byte
    std::mutex lock; // protects partial and length
    std::vector<std::pair<uint64_t, partial_sector>> partial; // the partly written sectors with their indexes

    // writes size bytes from an aligned buffer at an aligned offset
    void write_aligned(const uint8_t* data, uint64_t size, uint64_t offset)
    {
        if (!pwrite_fully(fd, (const char*) data, size, offset))
            throw std::runtime_error("pfp_rlbwt: cannot write " + path);
    }

    // adds the bytes data[0..end-beg-1] at the offsets [beg, end) to the partly written sectors; writes complete ones
    void write_partial(const uint8_t* data, uint64_t beg, uint64_t end)
    {
        std::lock_guard<std::mutex> guard(lock);

        while (beg < end) {
            const uint64_t index = beg / sector;
            const uint64_t num = std::min(end, (index + 1) * sector) - beg;
            auto it = std::find_if(partial.begin(), partial.end(), [&](const auto& e) { return e.first == index; });

            if (it == partial.end()) {
                uint8_t* bytes = (uint8_t*) allocate_aligned(sector);
                std::memset(bytes, 0, sector);
                partial.emplace_back(index, partial_sector { bytes, 0 });
                it = partial.end() - 1;
            }

            std::memcpy(it->second.bytes + (beg - index * sector), data, num);
            it->second.written += num;

            if (it->second.written == sector) {
                write_aligned(it->second.bytes, sector, index * sector);
                platform::free_aligned(it->second.bytes);
                partial.erase(it);
            }

            data += num;
            beg += num;
        }
    }

public:
    // creates (or truncates) the file
    explicit parallel_output_file(const std::string& path) : path(path)
    {
        fd = open_direct_or_buffered(path, platform::file_access::write, direct);
        if (fd < 0) throw std::runtime_error("pfp_rlbwt: cannot write " + path);
        platform::make_sparse(fd);
    }

    ~parallel_output_file()
    {
        for (auto& [index, s] : partial) platform::free_aligned(s.bytes);
        if (fd >= 0) platform::close(fd);
    }

    parallel_output_file(const parallel_output_file&) = delete;
    parallel_output_file& operator=(const parallel_output_file&) = delete;

    /**
     * @brief writes bytes to the file
     * @param data the bytes
     * @param size number of bytes (> 0)
     * @param offset offset in the file
     */
    void write(const void* data, uint64_t size, uint64_t offset)
    {
        const uint8_t* bytes = (const uint8_t*) data;
        const uint64_t end = offset + size;

        {
            std::lock_guard<std::mutex> guard(lock);
            length = std::max(length, end);
        }

        if (!direct) {
            write_aligned(bytes, size, offset);
            return;
        }

        // [first, last) are the sectors the write covers completely
        const uint64_t first = std::min(end, div_ceil(offset, sector) * sector);
        const uint64_t last = std::max(first, (end / sector) * sector);
        if (offset < first) write_partial(bytes, offset, first);
        if (last < end) write_partial(bytes + (last - offset), last, end);
        if (first == last) return;

        // the complete sectors are copied to an aligned buffer of each thread
        static thread_local std::unique_ptr<uint8_t, decltype(&platform::free_aligned)> buffer(
            (uint8_t*) allocate_aligned(write_buffer_bytes), &platform::free_aligned);

        for (uint64_t pos = first; pos < last; pos += write_buffer_bytes) {
            const uint64_t num = std::min(write_buffer_bytes, last - pos);
            std::memcpy(buffer.get(), bytes + (pos - offset), num);
            write_aligned(buffer.get(), num, pos);
        }
    }

    // writes the partly written sectors and truncates the file to its length
    void close()
    {
        for (auto& [index, s] : partial) {
            write_aligned(s.bytes, sector, index * sector);
            platform::free_aligned(s.bytes);
        }

        partial.clear();
        if (!platform::set_size(fd, length)) throw std::runtime_error("pfp_rlbwt: cannot write " + path);
        platform::close(fd);
        fd = -1;
    }
};

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

/**
 * @brief reads the values of an array stored in a file (a fixed number of bytes per value, little endian), with
 *        O_DIRECT (if the file system does not support it, through the page cache, whose read pages are dropped); a
 *        value that is not in the buffer is read together with the values after it, read_bytes bytes at once
 */
class array_file_reader {
private:
    static constexpr uint64_t sector = direct_io_alignment;

    uint64_t width; // bytes per value
    uint64_t mask; // the low 8 * width bits
    uint64_t read_bytes; // number of bytes read at once
    uint64_t num_values = 0; // number of values in the file
    uint8_t* buffer = nullptr; // read_bytes + 2 * sector bytes (rounded up to a multiple of sector), aligned to sector
    int fd = -1; // file descriptor
    bool direct = false; // true <=> the file is opened with O_DIRECT
    uint64_t beg = 0, end = 0; // the values [beg, end) are in the buffer
    const uint8_t* data = nullptr; // the bytes of value beg in the buffer

    // reads the values [i, i + read_bytes / width) into the buffer
    void read_from(uint64_t i)
    {
        beg = i;
        end = std::min(num_values, i + std::max<uint64_t>(1, read_bytes / width));
        const uint64_t first = beg * width, last = end * width; // the bytes [first, last) are needed
        const uint64_t read_beg = direct ? (first / sector) * sector : first;
        const uint64_t read_end = direct ? div_ceil(last, sector) * sector : last;
        const uint64_t num_read = pread_fully(fd, buffer, read_end - read_beg, read_beg);
        if (read_beg + num_read < last) throw std::runtime_error("pfp_rlbwt: cannot read a temporary file");
        if (!direct) platform::drop_cached(fd, read_beg, num_read);
        data = buffer + (first - read_beg);
    }

public:
    /**
     * @brief opens a file of values
     * @param path path of the file
     * @param width bytes per value
     * @param read_bytes number of bytes read at once (at least width)
     */
    array_file_reader(const std::string& path, uint64_t width, uint64_t read_bytes = direct_reader::block)
        : width(width), mask(width >= sizeof(uint64_t) ? ~uint64_t { 0 } : (uint64_t { 1 } << (8 * width)) - 1),
          read_bytes(read_bytes)
    {
        std::error_code ec;
        num_values = std::filesystem::file_size(path, ec) / width;
        if (ec) throw std::runtime_error("pfp_rlbwt: cannot read " + path);
        buffer = (uint8_t*) platform::allocate_aligned_counted(div_ceil(read_bytes, sector) * sector + 2 * sector,
            sector);
        if (buffer == nullptr) throw std::bad_alloc();
        fd = open_direct_or_buffered(path, platform::file_access::read, direct);

        if (fd < 0) {
            platform::free_aligned(buffer);
            throw std::runtime_error("pfp_rlbwt: cannot read " + path);
        }

        if (!direct) advise_uncached(fd);
    }

    ~array_file_reader()
    {
        platform::free_aligned(buffer);
        platform::close(fd);
    }

    array_file_reader(const array_file_reader&) = delete;
    array_file_reader& operator=(const array_file_reader&) = delete;

    // returns the bytes of value i
    inline const uint8_t* bytes(uint64_t i)
    {
        if (i < beg || i >= end) [[unlikely]] read_from(i);
        return data + (i - beg) * width;
    }

    // returns value i (width <= 8; loads 8 bytes, which the buffer holds after each value, and masks them)
    inline uint64_t operator[](uint64_t i)
    {
        uint64_t val;
        std::memcpy(&val, bytes(i), sizeof(uint64_t));
        return val & mask;
    }

    // copies the values [from, to) to dst
    void read(uint64_t from, uint64_t to, void* dst)
    {
        uint8_t* out = (uint8_t*) dst;

        for (uint64_t i = from; i < to;) {
            const uint8_t* src = bytes(i);
            const uint64_t num = std::min(to, end) - i;
            std::memcpy(out, src, num * width);
            out += num * width;
            i += num;
        }
    }
};

} // namespace pfp_rlbwt

#endif // PFP_RLBWT_PFP_RLBWT_CONSTRUCTION_COMMON_HPP
