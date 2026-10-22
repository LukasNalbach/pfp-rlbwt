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

#ifndef PFP_RLBWT_PFP_RLBWT_OUTPUT_RUN_DATA_FILE_TPP
#define PFP_RLBWT_PFP_RLBWT_OUTPUT_RUN_DATA_FILE_TPP

#include <pfp_rlbwt/output/run_data_file.hpp>

namespace pfp_rlbwt {

inline bool run_data_file::open(const std::string& file_name)
{
    path = file_name;
    entries.clear();
    bool valid = false;

    try {
        direct_reader file(file_name);
        const uint64_t file_size = file.size();
        uint64_t num_entries = 0;
        if (file_size < sizeof(num_entries)) return false;
        std::memcpy(&num_entries, file.read(file_size - sizeof(num_entries), file_size), sizeof(num_entries));
        if (num_entries > (file_size - sizeof(num_entries)) / sizeof(run_data_entry)) return false;
        // the sections end before the table of contents
        const uint64_t contents_offset = file_size - sizeof(num_entries) - num_entries * sizeof(run_data_entry);
        entries.resize(num_entries);

        if (num_entries > 0) {
            std::memcpy(entries.data(), file.read(contents_offset, file_size - sizeof(num_entries)),
                num_entries * sizeof(run_data_entry));
        }

        valid = true;

        for (uint64_t i = 0; valid && i < num_entries; i++) {
            valid = entries[i].offset <= contents_offset && entries[i].size <= contents_offset - entries[i].offset;
            for (uint64_t j = 0; valid && j < i; j++) valid = entries[j].id != entries[i].id;
        }

        valid = valid && read_contents();
    } catch (const std::exception&) {
        valid = false;
    }

    if (!valid) entries.clear();
    return valid;
}

inline bool run_data_file::read_section(uint64_t id, std::string& bytes) const
{
    const run_data_entry* e = find(id);
    if (e == nullptr) return false;
    bytes.clear();
    if (e->size == 0) return true;

    try {
        direct_reader file(path);
        const uint8_t* data = file.read(e->offset, e->offset + e->size);
        bytes.assign((const char*) data, e->size);
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

inline bool run_data_file::read_padded(uint64_t id, no_init_vector<uint8_t>& bytes) const
{
    const run_data_entry* e = find(id);
    if (e == nullptr) return false;
    resize_no_init(bytes, e->size + read_padding);
    std::memset(bytes.data() + e->size, 0, read_padding);
    if (e->size == 0) return true;
    direct_reader file(path);
    std::memcpy(bytes.data(), file.read(e->offset, e->offset + e->size), e->size);
    return true;
}

inline bool run_data_file::read_contents()
{
    std::string bytes;
    if (!read_section(section_fields, bytes) || bytes.size() != sizeof(run_data_fields)) return false;
    std::memcpy(&values, bytes.data(), sizeof(run_data_fields));
    const run_data_fields& f = values;
    if (f.r == 0 || f.r > f.n || f.runs_per_block == 0 || f.sa_samples > 1 || f.lcp_samples > 1) return false;
    blocks = div_ceil(f.r, f.runs_per_block);
    const bool plain = contains(section_records);
    fmt = plain ? file_format::plain : file_format::compressed;

    if (plain) {
        if (f.position_bytes == 0 || f.position_bytes > sizeof(uint64_t) || f.symbol_bytes == 0
            || f.symbol_bytes > sizeof(uint64_t)) {
            return false;
        }

        const record_layout layout(f);
        return section(section_records).size == f.r * layout.bytes && read_padded(section_record_samples, run_index)
            && run_index.size() == blocks * f.position_bytes + read_padding;
    }

    const uint64_t sa_bytes = div_ceil<uint64_t>(f.r * f.position_bits, byte_bits);
    if (f.position_bits == 0 || f.position_bits > max_packed_bits || f.run_block_offset_bits > max_packed_bits
        || f.lcp_block_offset_bits > max_packed_bits || f.runs_per_symbol_sample == 0
        || f.symbol_bytes == 0 || f.symbol_bytes > sizeof(uint64_t) || !contains(section_run_blocks)) {
        return false;
    }

    if (f.sa_samples != 0 && (!contains(section_sa_run_starts) || !contains(section_sa_run_ends)
        || section(section_sa_run_starts).size != sa_bytes || section(section_sa_run_ends).size != sa_bytes)) {
        return false;
    }

    if (f.lcp_samples != 0 && (!contains(section_lcp_blocks) || !read_padded(section_lcp_block_index, lcp_index)
        || lcp_index.size() != div_ceil<uint64_t>(blocks * f.lcp_block_offset_bits, byte_bits) + read_padding)) {
        return false;
    }

    index_entry_bits = f.run_block_offset_bits + f.position_bits;
    no_init_vector<uint8_t> code_bytes;
    return read_padded(section_run_block_index, run_index)
        && run_index.size() == div_ceil<uint64_t>(blocks * index_entry_bits, byte_bits) + read_padding
        && read_padded(section_symbol_code, code_bytes)
        && code.load(code_bytes.data(), code_bytes.size() - read_padding, uint8_t(f.symbol_bytes));
}

} // namespace pfp_rlbwt

#endif // PFP_RLBWT_PFP_RLBWT_OUTPUT_RUN_DATA_FILE_TPP
