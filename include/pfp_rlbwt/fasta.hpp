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
 * @brief the text of a FASTA file as the input of pfp_rlbwt::build
 */

#ifndef PFP_RLBWT_PFP_RLBWT_FASTA_HPP
#define PFP_RLBWT_PFP_RLBWT_FASTA_HPP

#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <omp.h>

#include <fasta/scan.hpp>
#include <misc/files.hpp>
#include <misc/memory.hpp>
#include <pfp_rlbwt/text_ranges.hpp>
#include <pfp_rlbwt/types.hpp>
#include <rlbwt/fasta_source.hpp>

namespace pfp_rlbwt {

// the transformation, the scan and the readers of the text of a FASTA file of string-utils (fasta/scan.hpp)
using ::append_fasta_chunk;
using ::fasta_carriage_return;
using ::fasta_chunk;
using ::fasta_first_header;
using ::fasta_header_marker;
using ::fasta_line_feed;
using ::fasta_min_chunk_bytes;
using ::fasta_range_bound;
using ::fasta_range_reader;
using ::fasta_range_starts;
using ::fasta_sample;
using ::fasta_sample_period;
using ::fasta_scan;
using ::fasta_scan_params;
using ::fasta_text_reader;
using ::fasta_transform;
using ::fasta_window;
using ::memory_file;
using ::read_ahead_reader;
using ::scan_fasta_lines;
using ::scan_fasta_range;

// the text of a FASTA file as a text_source of string-utils (rlbwt/fasta_source.hpp)
using ::rlbwt::fasta_source;

/**
 * @brief the text of a FASTA file (see fasta_transform) as the input of build, read once: range i of p starts at the
 *        first offset at or after i * size / p at which the state of the transformation follows from its line (see
 *        fasta_range_start), and each range is scanned while it is parsed (see fasta_range_reader). After the parse,
 *        the ranges are appended to the text (its length and samples, see scan), whose symbol frequencies the parse
 *        counts, and the symbols are renamed by a map computed from them.
 * @tparam chunk_t type of the result of the scan of a range (derived from fasta_chunk)
 */
template <typename chunk_t = fasta_chunk>
class fasta_ranges : public text_ranges<uint8_t> {
public:
    // returns the map that renames the symbols of the text, from the result of the scan (see scan)
    using map_fnc_t = std::function<std::array<uint8_t, 256>(const fasta_scan&)>;
    // called as merge(chunk, shift) for each range in the order of the ranges (see scan_fasta_lines)
    using merge_fnc_t = std::function<void(chunk_t&, uint64_t)>;

private:
    // reads a range of the text (see fasta_range_reader)
    template <typename bytes_t>
    class range_reader : public text_ranges<uint8_t>::reader {
    private:
        fasta_range_reader<bytes_t, chunk_t> in; // reads the text
        const chunk_t& chunk; // the result of the scan of the range
        no_init_vector<uint8_t> buffer; // the symbols of the last block

    public:
        range_reader(bytes_t&& file, const fasta_transform& transform, const fasta_range_bound& begin, uint64_t end,
            chunk_t& chunk, uint64_t sample_period, phase_progress* progress)
            : in(std::move(file), transform, begin, end, chunk, sample_period, progress), chunk(chunk) { }

        const uint8_t* next(uint64_t max_len, uint64_t& len) override
        {
            if (buffer.size() < max_len) resize_no_init(buffer, max_len);
            len = in.read(buffer.data(), max_len);
            return buffer.data();
        }

        uint64_t range_length() const override { return in.range_read() ? chunk.text_size : unknown_length; }
    };

    std::string path; // path of the FASTA file (empty if it is in memory)
    const uint8_t* data = nullptr; // the FASTA file in memory (nullptr if it is read from path)
    uint64_t file_size = 0; // size of the FASTA file
    fasta_transform transform; // the transformation
    map_fnc_t map_of; // computes the map of the symbols (none if empty)
    merge_fnc_t merge; // is called with the results of the ranges (none if empty)
    uint64_t sample_period; // a sample every sample_period positions of a range (0: none)
    std::vector<fasta_range_bound> bounds; // [0..p] the starts of the ranges; bounds[p] is at the end of the file
    std::vector<chunk_t> chunks; // [0..p-1] the results of the scans of the ranges
    // the length, the symbol frequencies and the samples of the text, once the ranges have been read
    fasta_scan text;
    std::array<uint8_t, 256> map {}; // the map of the symbols, once the ranges have been read
    // whether the ranges have been appended to text since the last split (a repeated parse reads them again)
    bool merged = false;

    // returns a new reader of the bytes of the FASTA file
    auto open_file() const { return direct_reader(path); }
    memory_file open_memory() const { return memory_file(data, file_size); }

    // sets the starts of p ranges (see fasta_range_starts)
    template <typename open_fnc_t>
    void set_bounds(open_fnc_t open, uint64_t p)
    {
        std::vector<uint64_t> positions(p);
        for (uint64_t i = 0; i < p; i++) positions[i] = (i * file_size) / p;
        bounds = fasta_range_starts(open, transform, positions, p);
        bounds.push_back(fasta_range_bound { file_size, {} });

        // a range that would start before the range before it is empty
        for (uint64_t i = 1; i < p; i++) {
            if (bounds[i].offset < bounds[i - 1].offset) bounds[i] = bounds[i - 1];
        }
    }

public:
    fasta_ranges(const fasta_ranges&) = delete;
    fasta_ranges& operator=(const fasta_ranges&) = delete;

    /**
     * @brief the text of a FASTA file
     * @param path path of the FASTA file
     * @param transform the transformation
     * @param map_of computes the map of the symbols from the result of the scan (the symbols keep their values if it
     *        is empty)
     * @param merge is called with the result of the scan of each range (see merge_fnc_t; nothing if empty)
     * @param sample_period a sample every sample_period positions of a range (0: none, see fasta_sample)
     */
    fasta_ranges(const std::string& path, const fasta_transform& transform, map_fnc_t map_of = {},
        merge_fnc_t merge = {}, uint64_t sample_period = fasta_sample_period)
        : path(path), file_size(open_file().size()), transform(transform), map_of(std::move(map_of)),
          merge(std::move(merge)), sample_period(sample_period) { }

    /**
     * @brief the text of a FASTA file in memory (see the constructor above)
     * @param data the FASTA file (not copied; it must stay valid until the construction has finished)
     * @param size the size of the FASTA file
     */
    fasta_ranges(const uint8_t* data, uint64_t size, const fasta_transform& transform, map_fnc_t map_of = {},
        merge_fnc_t merge = {}, uint64_t sample_period = fasta_sample_period)
        : data(data), file_size(size), transform(transform), map_of(std::move(map_of)), merge(std::move(merge)),
          sample_period(sample_period) { }

    uint64_t max_size() const override { return file_size; }

    bool size_known() const override { return false; }

    void split(uint64_t p) override
    {
        if (data != nullptr) set_bounds([&] { return open_memory(); }, p);
        else set_bounds([&] { return open_file(); }, p);
        chunks.assign(p, chunk_t());
        merged = false;
    }

    std::unique_ptr<text_ranges<uint8_t>::reader> open(uint64_t i, phase_progress* progress) override
    {
        if (data != nullptr) {
            return std::make_unique<range_reader<memory_file>>(open_memory(), transform, bounds[i],
                bounds[i + 1].offset, chunks[i], sample_period, progress);
        }

        return std::make_unique<range_reader<read_ahead_reader>>(read_ahead_reader(path), transform, bounds[i],
            bounds[i + 1].offset, chunks[i], sample_period, progress);
    }

    bool needs_frequencies() const override { return true; }

    const std::array<uint8_t, 256>* ranges_read(const std::array<uint64_t, 256>* freq) override
    {
        if (!merged) {
            text = fasta_scan();
            const uint8_t start_state = fasta_transform::state_t {}.to_bits();
            if (sample_period != 0) text.samples.push_back(fasta_sample { 0, 0, start_state });

            for (chunk_t& chunk : chunks) {
                const uint64_t shift = append_fasta_chunk(text, chunk, transform.separator);
                if (merge) merge(chunk, shift);
            }

            text.freq = *freq;
            if (map_of) map = map_of(text);
            merged = true;
        }

        return map_of ? &map : nullptr;
    }

    // returns the length, the symbol frequencies and the samples of the text (once build has read it)
    const fasta_scan& scan() const { return text; }
};

} // namespace pfp_rlbwt

#endif // PFP_RLBWT_PFP_RLBWT_FASTA_HPP
