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
 * @file fasta.hpp
 * @brief the text of a FASTA file as the input of build (see fasta_source); the transformation, the scan and the
 *        reader of the text are the ones of string-utils (fasta/scan.hpp)
 */

#ifndef PFP_RLBWT_PFP_RLBWT_FASTA_HPP
#define PFP_RLBWT_PFP_RLBWT_FASTA_HPP

#include <algorithm>
#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <fasta/scan.hpp>
#include <misc/files.hpp>
#include <misc/memory.hpp>
#include <pfp_rlbwt/types.hpp>

namespace pfp_rlbwt {

using ::append_fasta_chunk;
using ::fasta_carriage_return;
using ::fasta_chunk;
using ::fasta_header_marker;
using ::fasta_line_feed;
using ::fasta_min_chunk_bytes;
using ::fasta_sample;
using ::fasta_sample_period;
using ::fasta_scan;
using ::fasta_scan_params;
using ::fasta_text_reader;
using ::fasta_transform;
using ::fasta_window;
using ::memory_file;
using ::scan_fasta_lines;
using ::scan_fasta_range;

/**
 * @brief the text of a FASTA file (see fasta_transform) as the input of build, optionally reversed and with its
 *        symbols mapped; a reader starts at the last sample at or before the position it reads (see fasta_sample)
 * @tparam sym_t type of a symbol
 */
template <typename sym_t = uint8_t>
class fasta_source : public text_source<sym_t> {
private:
    // reads the text; consecutive reads of consecutive ranges continue the scan of the file
    template <typename bytes_t>
    class fasta_reader : public text_source<sym_t>::reader {
    public:
        fasta_reader(const fasta_source& src, bytes_t&& file)
            : src(src), text(std::move(file), src.transform, src.samples, src.num_samples) { }

        void read(uint64_t pos, uint64_t len, sym_t* dst) override
        {
            const uint64_t beg = src.reverse ? src.text_size - pos - len : pos;
            uint8_t* bytes; // the symbols of the text before they are reversed and mapped

            if constexpr (std::is_same_v<sym_t, uint8_t>) {
                bytes = dst;
            } else {
                if (buffer.size() < len) resize_no_init(buffer, len);
                bytes = buffer.data();
            }

            text.read(beg, len, bytes);
            if (src.reverse) std::reverse(bytes, bytes + len);

            if (src.mapped) {
                for (uint64_t i = 0; i < len; i++) dst[i] = src.map[bytes[i]];
            } else if constexpr (!std::is_same_v<sym_t, uint8_t>) {
                for (uint64_t i = 0; i < len; i++) dst[i] = sym_t(bytes[i]);
            }
        }

    private:
        const fasta_source& src;
        fasta_text_reader<bytes_t> text; // reads the text
        no_init_vector<uint8_t> buffer; // the bytes of a read, if sym_t is not uint8_t
    };

    std::string path; // path of the FASTA file (empty if it is in memory)
    const uint8_t* data = nullptr; // the FASTA file in memory (nullptr if it is read from path)
    uint64_t data_size = 0; // size of the FASTA file in memory
    uint64_t text_size; // length of the text
    fasta_transform transform;
    std::vector<fasta_sample> computed_samples; // the samples, if they have been computed here
    const fasta_sample* samples; // [0..num_samples-1] the samples, in increasing order of their positions
    uint64_t num_samples;
    bool reverse; // true <=> the text is read reversed
    bool mapped; // true <=> map is applied to each symbol
    std::array<sym_t, 256> map {}; // [0..255] map[c] = the symbol the byte c of the text is read as

    // returns a new reader of the bytes of the FASTA file
    auto open_file() const { return direct_reader(path); }
    memory_file open_memory() const { return memory_file(data, data_size); }

    // sets the samples, the map and computes the samples if none are given
    template <typename open_fnc_t>
    void init(open_fnc_t open, const std::array<sym_t, 256>* map, uint64_t num_threads)
    {
        if (mapped) this->map = *map;

        if (samples == nullptr || num_samples == 0) {
            fasta_scan scan = scan_fasta_lines<fasta_chunk>(open, transform,
                { .num_threads = num_threads }, [](auto&&...) { });
            computed_samples = std::move(scan.samples);
            samples = computed_samples.data();
            num_samples = computed_samples.size();
        }
    }

public:
    fasta_source(const fasta_source&) = delete;
    fasta_source& operator=(const fasta_source&) = delete;

    /**
     * @brief the text of a FASTA file
     * @param path path of the FASTA file
     * @param text_size length of the text
     * @param transform the transformation
     * @param samples the samples of the text in increasing order of their positions (computed here if nullptr)
     * @param num_samples number of samples
     * @param reverse true <=> the text is read reversed
     * @param map if not nullptr, the symbol c of the text is read as (*map)[c]
     * @param num_threads number of threads used to compute the samples
     */
    fasta_source(const std::string& path, uint64_t text_size, const fasta_transform& transform,
        const fasta_sample* samples, uint64_t num_samples, bool reverse = false,
        const std::array<sym_t, 256>* map = nullptr, uint64_t num_threads = 1)
        : path(path), text_size(text_size), transform(transform), samples(samples), num_samples(num_samples),
          reverse(reverse), mapped(map != nullptr)
    {
        init([&] { return open_file(); }, map, num_threads);
    }

    /**
     * @brief the text of a FASTA file in memory (see the constructor above)
     * @param data the FASTA file (not copied; it must stay valid until the construction has finished)
     * @param size the size of the FASTA file
     */
    fasta_source(const uint8_t* data, uint64_t size, uint64_t text_size, const fasta_transform& transform,
        const fasta_sample* samples, uint64_t num_samples, bool reverse = false,
        const std::array<sym_t, 256>* map = nullptr, uint64_t num_threads = 1)
        : data(data), data_size(size), text_size(text_size), transform(transform), samples(samples),
          num_samples(num_samples), reverse(reverse), mapped(map != nullptr)
    {
        init([&] { return open_memory(); }, map, num_threads);
    }

    uint64_t size() const override { return text_size; }

    std::unique_ptr<typename text_source<sym_t>::reader> new_reader() const override
    {
        if (data != nullptr) return std::make_unique<fasta_reader<memory_file>>(*this, open_memory());
        return std::make_unique<fasta_reader<direct_reader>>(*this, open_file());
    }
};

} // namespace pfp_rlbwt

#endif // PFP_RLBWT_PFP_RLBWT_FASTA_HPP
