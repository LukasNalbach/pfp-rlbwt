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

#ifndef PFP_RLBWT_PFP_RLBWT_OUTPUT_SYMBOL_CODE_TPP
#define PFP_RLBWT_PFP_RLBWT_OUTPUT_SYMBOL_CODE_TPP

#include <pfp_rlbwt/output/symbol_code.hpp>

#include <algorithm>
#include <functional>
#include <queue>

namespace pfp_rlbwt {

inline symbol_code::symbol_code(std::vector<std::pair<uint64_t, uint64_t>> frequencies)
{
    const uint64_t m = frequencies.size();
    // [0..m-1] the code length of each symbol
    std::vector<uint8_t> lengths(m, 1);

    if (m > 1) {
        std::vector<uint64_t> weight(m);
        for (uint64_t i = 0; i < m; i++) weight[i] = frequencies[i].second;
        /* the nodes of the Huffman tree: [0..m-1] the leaves, [m..2m-2] the inner nodes in the order of their creation;
           parent[v] = the parent of node v, depth[v] = the depth of node v */
        std::vector<uint64_t> parent(2 * m - 1), depth(2 * m - 1);

        // a code longer than max_code_bits is avoided by halving the weights (rounded up) until none is longer
        while (true) {
            using node_t = std::pair<uint64_t, uint64_t>; // (weight, node)
            std::priority_queue<node_t, std::vector<node_t>, std::greater<node_t>> queue;
            for (uint64_t i = 0; i < m; i++) queue.emplace(weight[i], i);
            uint64_t next = m;

            while (queue.size() > 1) {
                const node_t a = queue.top();
                queue.pop();
                const node_t b = queue.top();
                queue.pop();
                parent[a.second] = parent[b.second] = next;
                queue.emplace(a.first + b.first, next++);
            }

            depth[2 * m - 2] = 0;
            for (uint64_t v = 2 * m - 2; v-- > 0;) depth[v] = depth[parent[v]] + 1;
            const uint64_t longest = *std::max_element(depth.begin(), depth.begin() + m);

            if (longest <= max_code_bits) {
                for (uint64_t i = 0; i < m; i++) lengths[i] = uint8_t(depth[i]);
                break;
            }

            for (uint64_t& w : weight) w = (w + 1) / 2;
        }
    }

    // the symbols in the order of their codes: by code length, then by symbol
    std::vector<uint64_t> order(m);
    for (uint64_t i = 0; i < m; i++) order[i] = i;
    std::sort(order.begin(), order.end(), [&](uint64_t a, uint64_t b) {
        return lengths[a] != lengths[b] ? lengths[a] < lengths[b] : frequencies[a].first < frequencies[b].first;
    });

    max_length = *std::max_element(lengths.begin(), lengths.end());
    count.assign(max_length + 1, 0);
    symbols.resize(m);

    for (uint64_t i = 0; i < m; i++) {
        symbols[i] = frequencies[order[i]].first;
        count[lengths[order[i]]]++;
    }

    initialize();
}

inline void symbol_code::initialize()
{
    const uint64_t m = symbols.size();
    first_code.assign(max_length + 1, 0);
    first_index.assign(max_length + 1, 0);
    uint64_t code = 0, index = 0;

    for (uint8_t l = 1; l <= max_length; l++) {
        first_code[l] = code;
        first_index[l] = index;
        code = (code + count[l]) << 1;
        index += count[l];
    }

    stream_code.resize(m);
    length.resize(m);
    table_bits = std::min(max_length, lookup_bits);
    table.assign(uint64_t { 1 } << table_bits, 0);
    uint8_t l = 1;

    for (uint64_t i = 0; i < m; i++) {
        while (i >= first_index[l] + count[l]) l++;
        const uint64_t c = first_code[l] + (i - first_index[l]);
        uint64_t reversed = 0;
        for (uint8_t j = 0; j < l; j++) reversed |= ((c >> j) & 1) << (l - 1 - j);
        stream_code[i] = reversed;
        length[i] = l;

        if (l <= table_bits) {
            for (uint64_t f = 0; f < (uint64_t { 1 } << (table_bits - l)); f++) {
                table[reversed | (f << l)] = (i << byte_bits) | l;
            }
        }
    }

    // the symbols that fit into two bytes are indexed by an array, the others by a hash table
    const uint64_t dense_limit = uint64_t { std::numeric_limits<uint16_t>::max() } + 1;
    uint64_t dense_size = 0;
    for (uint64_t sym : symbols) if (sym < dense_limit) dense_size = std::max(dense_size, sym + 1);
    dense_index.assign(dense_size, 0);
    sparse_index.clear();

    for (uint64_t i = 0; i < m; i++) {
        if (symbols[i] < dense_size) dense_index[symbols[i]] = i;
        else sparse_index.emplace(symbols[i], i);
    }
}

inline bool symbol_code::load(const uint8_t* data, uint64_t size, uint8_t symbol_bytes)
{
    if (size < sizeof(uint64_t)) return false;
    const uint64_t max_len = read_bytes(data, sizeof(uint64_t));
    if (max_len == 0 || max_len > max_code_bits || size < (max_len + 1) * sizeof(uint64_t)) return false;
    max_length = uint8_t(max_len);
    count.assign(max_length + 1, 0);
    uint64_t m = 0;
    // the codes satisfy the Kraft inequality: sum of 2^(max_length - l) over the codes <= 2^max_length
    uint64_t kraft = 0;

    for (uint8_t l = 1; l <= max_length; l++) {
        count[l] = read_bytes(data + l * sizeof(uint64_t), sizeof(uint64_t));
        if (count[l] > (uint64_t { 1 } << l)) return false;
        m += count[l];
        kraft += count[l] << (max_length - l);
    }

    const uint64_t header = (uint64_t { max_length } + 1) * sizeof(uint64_t);
    if (m == 0 || kraft > (uint64_t { 1 } << max_length) || (size - header) / symbol_bytes != m
        || (size - header) % symbol_bytes != 0) {
        return false;
    }

    symbols.resize(m);
    for (uint64_t i = 0; i < m; i++) symbols[i] = read_bytes(data + header + i * symbol_bytes, symbol_bytes);
    initialize();
    return true;
}

inline std::string symbol_code::serialize(uint8_t symbol_bytes) const
{
    std::string bytes((uint64_t { max_length } + 1) * sizeof(uint64_t) + symbols.size() * symbol_bytes, 0);
    uint8_t* data = (uint8_t*) bytes.data();
    write_bytes(data, max_length, sizeof(uint64_t));
    for (uint8_t l = 1; l <= max_length; l++) write_bytes(data + l * sizeof(uint64_t), count[l], sizeof(uint64_t));
    const uint64_t header = (uint64_t { max_length } + 1) * sizeof(uint64_t);
    for (uint64_t i = 0; i < symbols.size(); i++) write_bytes(data + header + i * symbol_bytes, symbols[i], symbol_bytes);
    return bytes;
}

} // namespace pfp_rlbwt

#endif // PFP_RLBWT_PFP_RLBWT_OUTPUT_SYMBOL_CODE_TPP
