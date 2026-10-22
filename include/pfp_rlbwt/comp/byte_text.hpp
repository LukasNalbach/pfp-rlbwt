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
 * @file byte_text.hpp
 * @brief a byte text in a file or in memory: the bytes that occur in it, and the text with its bytes mapped to the
 *        symbols of pfp-rlbwt
 */

#ifndef PFP_RLBWT_PFP_RLBWT_COMP_BYTE_TEXT_HPP
#define PFP_RLBWT_PFP_RLBWT_COMP_BYTE_TEXT_HPP

#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>

#include <coding/byte_codes.hpp>
#include <misc/files.hpp>
#include <misc/memory.hpp>
#include <misc/threads.hpp>
#include <pfp_rlbwt/types.hpp>

namespace pfp_rlbwt::comp {

// [0..255] occurs[c] = whether the byte c occurs in the text
using byte_set = std::array<bool, num_byte_values>;

/**
 * @brief a byte text T[0..n-1], stored in a file or in memory
 */
class byte_text {
private:
    std::string path; // the file storing the text (empty if the text is in memory)
    const uint8_t* data = nullptr; // the text in memory (nullptr if it is in a file)
    uint64_t n = 0; // n, the length of the text

public:
    /**
     * @brief a text stored in a file
     * @param path the file
     */
    explicit byte_text(const std::string& path) : path(path)
    {
        std::error_code ec;
        n = std::filesystem::file_size(path, ec);
        if (ec) throw std::runtime_error("cannot read " + path);
    }

    /**
     * @brief a text in memory (not copied; it must stay valid while the object is used)
     * @param data the text
     * @param n its length
     */
    byte_text(const uint8_t* data, uint64_t n) : data(data), n(n) { }

    // returns n, the length of the text
    uint64_t size() const { return n; }

    // returns the file storing the text (empty if the text is in memory)
    const std::string& file() const { return path; }

    // returns the text in memory (nullptr if it is in a file)
    const uint8_t* bytes() const { return data; }

    /**
     * @brief calls fnc(block, len) for consecutive blocks of the bytes T[from..to-1], in ascending order
     * @param from the first position
     * @param to the position after the last one
     * @param fnc the function
     */
    template <typename fnc_t>
    void for_each_block(uint64_t from, uint64_t to, fnc_t fnc) const
    {
        if (from >= to) return;

        if (data != nullptr) {
            fnc(data + from, to - from);
            return;
        }

        bool direct = false;
        const int fd = open_direct_or_buffered(path, platform::file_access::read, direct);
        if (fd < 0) throw std::runtime_error("cannot read " + path);
        platform::close(fd);
        direct_range_reader reader;
        reader.open(path, direct, std::min<uint64_t>(huge_page_bytes, to - from));
        reader.for_each_block(from, to, fnc);
    }

    /**
     * @brief returns the bytes that occur in the text
     * @param num_threads the number of threads
     */
    byte_set occurring_bytes(uint64_t num_threads) const
    {
        byte_set occurs {};
        const uint64_t p = std::max<uint64_t>(1, std::min(num_threads, n));

        // thread t of a team of q threads scans T[t * n / q..(t + 1) * n / q - 1]
        run_on_threads(p, [&](uint64_t t, uint64_t q) {
            byte_set seen {};
            for_each_block((t * n) / q, ((t + 1) * n) / q, [&](const uint8_t* block, uint64_t len) {
                for (uint64_t i = 0; i < len; i++) seen[block[i]] = true;
            });

            #pragma omp critical
            for (uint64_t c = 0; c < num_byte_values; c++) occurs[c] = occurs[c] || seen[c];
        });

        return occurs;
    }
};

/**
 * @brief the text of pfp-rlbwt that consists of the bytes of a byte_text mapped to symbols of type sym_t
 * @tparam sym_t the symbol type (uint8_t or uint16_t)
 */
template <typename sym_t>
class mapped_byte_source : public pfp_rlbwt::text_source<sym_t> {
private:
    const byte_text& text; // the byte text
    std::array<sym_t, num_byte_values> map; // [0..255] map[c] = the symbol of the byte c

    // reads ranges of the mapped text
    class reader : public pfp_rlbwt::text_source<sym_t>::reader {
    private:
        const mapped_byte_source& source; // the mapped text
        direct_reader file; // reads the file storing the text (not opened if the text is in memory)

    public:
        explicit reader(const mapped_byte_source& source) : source(source)
        {
            if (source.text.bytes() == nullptr) file.open_file(source.text.file());
        }

        void read(uint64_t pos, uint64_t len, sym_t* dst) override
        {
            const uint8_t* bytes =
                source.text.bytes() != nullptr ? source.text.bytes() + pos : file.read(pos, pos + len);
            for (uint64_t i = 0; i < len; i++) dst[i] = source.map[bytes[i]];
        }
    };

public:
    /**
     * @brief maps the bytes of a text
     * @param text the text (it must stay valid while the object is used)
     * @param map [0..255] map[c] = the symbol of the byte c (at least pfp_rlbwt::min_symbol for each byte of the text)
     */
    mapped_byte_source(const byte_text& text, const std::array<sym_t, num_byte_values>& map) : text(text), map(map) { }

    uint64_t size() const override { return text.size(); }

    std::unique_ptr<typename pfp_rlbwt::text_source<sym_t>::reader> new_reader() const override
    {
        return std::make_unique<reader>(*this);
    }
};

} // namespace pfp_rlbwt::comp

#endif // PFP_RLBWT_PFP_RLBWT_COMP_BYTE_TEXT_HPP
