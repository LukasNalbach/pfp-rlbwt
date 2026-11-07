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
 * @file text_ranges.hpp
 * @brief the text split into ranges that threads read at the same time
 */

#ifndef PFP_RLBWT_PFP_RLBWT_TEXT_RANGES_HPP
#define PFP_RLBWT_PFP_RLBWT_TEXT_RANGES_HPP

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <memory>

#include <misc/log.hpp>
#include <pfp_rlbwt/types.hpp>

namespace pfp_rlbwt {

/**
 * @brief the text T split into p ranges that p threads read at the same time, each from the start of its range on
 *        (and past its end, as far as it needs); the length of a range is known once its reader has passed its end, so
 *        the length of T may be unknown before it has been read
 * @tparam sym_t type of a symbol (see STRING_UTILS_FOR_EACH_SYMBOL_TYPE)
 */
template <typename sym_t>
class text_ranges {
public:
    // returned by reader::range_length while the reader has not passed the end of its range
    static constexpr uint64_t unknown_length = std::numeric_limits<uint64_t>::max();

    // reads T from the start of a range on
    class reader {
    public:
        virtual ~reader() = default;

        /**
         * @brief returns the next symbols of T, valid until the next call; a block lies in the range or after it
         * @param max_len [1..] the maximum number of symbols
         * @param len (output) the number of symbols (0: T has ended)
         */
        virtual const sym_t* next(uint64_t max_len, uint64_t& len) = 0;

        // returns the number of symbols of the range if the reader has passed its end, else unknown_length
        virtual uint64_t range_length() const = 0;
    };

    virtual ~text_ranges() = default;

    // returns the length of T if it is known before T is read, else an upper bound of it
    virtual uint64_t max_size() const = 0;

    // returns whether max_size() is the length of T
    virtual bool size_known() const = 0;

    // splits T into p ranges (called before open)
    virtual void split(uint64_t p) = 0;

    /**
     * @brief returns a new reader of a range
     * @param i [0..p-1] the range
     * @param progress if not nullptr, counts what is read of the range (max_size() over all ranges)
     */
    virtual std::unique_ptr<reader> open(uint64_t i, phase_progress* progress) = 0;

    // returns whether ranges_read needs the frequencies of the symbols of T (only for symbols of one byte)
    virtual bool needs_frequencies() const { return false; }

    /**
     * @brief called after every range has been read to its end; returns the map that the symbols of T are renamed by
     *        (only for symbols of one byte), or nullptr if they keep their values
     * @param freq if needs_frequencies(), [0..255] freq[c] = the number of occurrences of c in T, else nullptr
     */
    virtual const std::array<uint8_t, 256>* ranges_read(const std::array<uint64_t, 256>* /* freq */) { return nullptr; }
};

/**
 * @brief the ranges of a text_source: range i of p is T[i*n/p .. (i+1)*n/p - 1]
 * @tparam sym_t type of a symbol
 */
template <typename sym_t>
class source_ranges : public text_ranges<sym_t> {
private:
    // reads T from the start of a range on through a reader of the text_source
    class range_reader : public text_ranges<sym_t>::reader {
    private:
        std::unique_ptr<typename text_source<sym_t>::reader> in; // reads the text
        uint64_t n; // the length of T
        uint64_t pos; // the next position to read
        uint64_t beg, end; // the range [beg, end) of T
        phase_progress* progress; // counts the positions of the range read (if not nullptr)

    public:
        range_reader(const text_source<sym_t>& text, uint64_t beg, uint64_t end, phase_progress* progress)
            : in(text.new_reader()), n(text.size()), pos(beg), beg(beg), end(end), progress(progress) { }

        const sym_t* next(uint64_t max_len, uint64_t& len) override
        {
            len = std::min(max_len, (pos < end ? end : n) - pos);
            if (len == 0) return nullptr;
            const sym_t* block = in->view(pos, len);
            if (progress != nullptr && pos < end) progress->advance(len);
            pos += len;
            return block;
        }

        uint64_t range_length() const override { return end - beg; }
    };

    const text_source<sym_t>& text; // the text
    uint64_t p = 1; // the number of ranges

public:
    explicit source_ranges(const text_source<sym_t>& text) : text(text) { }

    uint64_t max_size() const override { return text.size(); }

    bool size_known() const override { return true; }

    void split(uint64_t num_ranges) override { p = num_ranges; }

    std::unique_ptr<typename text_ranges<sym_t>::reader> open(uint64_t i, phase_progress* progress) override
    {
        const uint64_t n = text.size();
        return std::make_unique<range_reader>(text, (i * n) / p, ((i + 1) * n) / p, progress);
    }
};

} // namespace pfp_rlbwt

#endif // PFP_RLBWT_PFP_RLBWT_TEXT_RANGES_HPP
