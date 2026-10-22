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
 * @file bits.hpp
 * @brief reading and writing bit-packed values (the bits of a byte least significant bit first)
 */

#ifndef PFP_RLBWT_PFP_RLBWT_OUTPUT_BITS_HPP
#define PFP_RLBWT_PFP_RLBWT_OUTPUT_BITS_HPP

#include <bit>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

#include <pfp_rlbwt/types.hpp>

namespace pfp_rlbwt {

// number of bits of a byte
inline constexpr uint8_t byte_bits = std::numeric_limits<uint8_t>::digits;

// number of bytes a buffer of bit-packed data has after its last byte, so that read_bits can load 8 bytes anywhere
inline constexpr uint64_t read_padding = sizeof(uint64_t);

// bits per width (of a value or of a part of a block) stored in a block: they hold max_packed_bits
inline constexpr uint8_t width_bits = std::bit_width(max_packed_bits);

// returns 2^bits - 1 (bits <= 64)
inline constexpr uint64_t low_mask(uint8_t bits)
{
    return bits >= std::numeric_limits<uint64_t>::digits ? ~uint64_t { 0 } : (uint64_t { 1 } << bits) - 1;
}

/**
 * @brief returns the value of width bits at a bit offset of bit-packed data
 * @param data the data; read_padding bytes after the byte containing the last bit of the value must be readable
 * @param bit offset (in bits) of the value
 * @param width number of bits of the value (at most max_packed_bits)
 */
inline uint64_t read_bits(const uint8_t* data, uint64_t bit, uint8_t width)
{
    uint64_t word;
    std::memcpy(&word, data + bit / byte_bits, sizeof(uint64_t));
    return (word >> (bit % byte_bits)) & low_mask(width);
}

// returns the value of the low bytes bytes (at most 8) at data, little endian
inline uint64_t read_bytes(const uint8_t* data, uint8_t bytes)
{
    uint64_t value = 0;
    std::memcpy(&value, data, bytes);
    return value;
}

// stores the low bytes bytes (at most 8) of value at data, little endian
inline void write_bytes(uint8_t* data, uint64_t value, uint8_t bytes) { std::memcpy(data, &value, bytes); }

// returns the index of the k-th set bit (k from 0) of a word that has more than k set bits
inline uint8_t select_in_word(uint64_t word, uint64_t k)
{
    for (; k > 0; k--) word &= word - 1;
    return uint8_t(std::countr_zero(word));
}

/**
 * @brief a sequence of bits stored in 64-bit words, least significant bit first; written to a file, the bits of
 *        a byte are therefore also stored least significant bit first
 */
class bit_buffer {
public:
    static constexpr uint8_t word_bits = std::numeric_limits<uint64_t>::digits;

private:
    // the bits: bit j is bit j % 64 of words[j / 64], j in [0..num_bits-1]; the bits after them are 0
    std::vector<uint64_t> words;
    uint64_t num_bits = 0;

public:
    // appends the bits low bits of val (bits <= 64; the other bits of val must be 0)
    inline void put(uint64_t val, uint8_t bits)
    {
        if (bits == 0) return;
        const uint8_t used = num_bits % word_bits;
        if (used == 0) words.push_back(val);
        else {
            words.back() |= val << used;
            if (used + bits > word_bits) words.push_back(val >> (word_bits - used));
        }

        num_bits += bits;
    }

    // appends num 0-bits
    inline void put_zeros(uint64_t num)
    {
        for (; num >= word_bits; num -= word_bits) put(0, word_bits);
        put(0, uint8_t(num));
    }

    // returns the number of bits
    uint64_t size() const { return num_bits; }

    // returns the number of bytes of the bits (the last byte padded with 0-bits)
    uint64_t bytes() const { return (num_bits + byte_bits - 1) / byte_bits; }

    // returns the bits (bytes() bytes)
    const uint8_t* data() const { return (const uint8_t*) words.data(); }

    // pads the bits with 0-bits to a multiple of 8 bits
    void pad_to_byte() { put_zeros((byte_bits - num_bits % byte_bits) % byte_bits); }

    void clear()
    {
        words.clear();
        num_bits = 0;
    }
};

} // namespace pfp_rlbwt

#endif // PFP_RLBWT_PFP_RLBWT_OUTPUT_BITS_HPP
