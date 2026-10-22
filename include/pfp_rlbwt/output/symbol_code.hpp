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
 * @file symbol_code.hpp
 * @brief the canonical Huffman code of the run symbols of the compressed format (see section_symbol_code)
 */

#ifndef PFP_RLBWT_PFP_RLBWT_OUTPUT_SYMBOL_CODE_HPP
#define PFP_RLBWT_PFP_RLBWT_OUTPUT_SYMBOL_CODE_HPP

#include <cstdint>
#include <limits>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <pfp_rlbwt/output/bits.hpp>

namespace pfp_rlbwt {

/**
 * @brief a canonical Huffman code of the run symbols: the codes of length l are the consecutive integers
 *        first_code[l], first_code[l] + 1, ... assigned to the symbols of code length l in increasing order. A code is
 *        stored most significant bit first in the bit order of the stream (least significant bit of a byte first).
 *        Section layout: the largest code length L (64 bits), then for l in [1..L] the number of codes of length l
 *        (64 bits each), then the symbols in the order of their codes (symbol_bytes bytes each).
 */
class symbol_code {
public:
    // the largest code length; a code is read with read_bits
    static constexpr uint8_t max_code_bits = max_packed_bits;

    // codes of at most lookup_bits bits are decoded with one lookup in table
    static constexpr uint8_t lookup_bits = 11;

private:
    std::vector<uint64_t> symbols; // the symbols in the order of their codes
    std::vector<uint64_t> count; // [0..max_length] count[l] = number of codes of length l
    std::vector<uint64_t> first_code; // [0..max_length] first_code[l] = the smallest code of length l
    std::vector<uint64_t> first_index; // [0..max_length] first_index[l] = index in symbols of first_code[l]
    uint8_t max_length = 0; // the largest code length
    uint8_t table_bits = 0; // min(max_length, lookup_bits)
    /* [0..2^table_bits-1] table[x] = (index in symbols << byte_bits) | code length, if the code of length at most
       table_bits is a prefix of the next table_bits bits x of the stream; 0 otherwise */
    std::vector<uint64_t> table;
    // [0..symbols.size()-1] the code of each symbol, bit-reversed (in stream order), and its length
    std::vector<uint64_t> stream_code;
    std::vector<uint8_t> length;
    // the index in symbols of a symbol: in dense_index for symbols that fit into two bytes, else in sparse_index
    std::vector<uint64_t> dense_index;
    std::unordered_map<uint64_t, uint64_t> sparse_index;

    // computes first_code, first_index, table, stream_code and the indexes from symbols and count
    void initialize();

    // returns the index of a symbol in symbols
    inline uint64_t index_of(uint64_t sym) const
    {
        return sym < dense_index.size() ? dense_index[sym] : sparse_index.at(sym);
    }

public:
    symbol_code() = default;

    /**
     * @brief computes the code from the frequencies of the symbols
     * @param frequencies the (symbol, frequency) pairs of the distinct symbols (frequencies >= 1, at least one pair)
     */
    explicit symbol_code(std::vector<std::pair<uint64_t, uint64_t>> frequencies);

    /**
     * @brief reads the code from its section
     * @param data the bytes of the section
     * @param size number of bytes of the section
     * @param symbol_bytes bytes per symbol
     * @return whether the section stores a valid code
     */
    bool load(const uint8_t* data, uint64_t size, uint8_t symbol_bytes);

    /**
     * @brief returns the bytes of the section of the code
     * @param symbol_bytes bytes per symbol
     */
    std::string serialize(uint8_t symbol_bytes) const;

    // returns the number of distinct symbols
    uint64_t num_symbols() const { return symbols.size(); }

    // returns the length of the code of a symbol
    inline uint8_t code_length(uint64_t sym) const { return length[index_of(sym)]; }

    // appends the code of a symbol to bits
    inline void encode(bit_buffer& bits, uint64_t sym) const
    {
        const uint64_t i = index_of(sym);
        bits.put(stream_code[i], length[i]);
    }

    /**
     * @brief decodes the code at a bit offset
     * @param data the bit-packed data (see read_bits)
     * @param bit offset (in bits) of the code; it is advanced by the code length
     * @return the symbol
     */
    inline uint64_t decode(const uint8_t* data, uint64_t& bit) const
    {
        const uint64_t window = read_bits(data, bit, max_length);
        const uint64_t entry = table[window & low_mask(table_bits)];

        if (entry != 0) [[likely]] {
            bit += entry & low_mask(byte_bits);
            return symbols[entry >> byte_bits];
        }

        uint64_t code = 0;

        for (uint8_t l = 1;; l++) {
            code = (code << 1) | ((window >> (l - 1)) & 1);

            if (code - first_code[l] < count[l]) {
                bit += l;
                return symbols[first_index[l] + code - first_code[l]];
            }
        }
    }

    // skips the code at a bit offset (see decode)
    inline void skip(const uint8_t* data, uint64_t& bit) const { decode(data, bit); }
};

} // namespace pfp_rlbwt

#include <pfp_rlbwt/output/symbol_code.tpp>

#endif // PFP_RLBWT_PFP_RLBWT_OUTPUT_SYMBOL_CODE_HPP
