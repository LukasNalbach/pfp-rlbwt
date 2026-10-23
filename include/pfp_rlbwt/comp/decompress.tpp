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
 * @file decompress.tpp
 * @brief decompress: decoding the runs, building M_LF, reverting the pieces of the text between the start pairs and
 *        restoring a FASTA file
 */

#ifndef PFP_RLBWT_PFP_RLBWT_COMP_DECOMPRESS_TPP
#define PFP_RLBWT_PFP_RLBWT_COMP_DECOMPRESS_TPP

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <condition_variable>
#include <cstring>
#include <filesystem>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <utility>
#include <vector>

#include <coding/byte_codes.hpp>
#include <coding/range_coder.hpp>
#include <data_structures/interleaved_bit_aligned_vectors.hpp>
#include <fasta/restore.hpp>
#include <fasta/split.hpp>
#include <misc/bits.hpp>
#include <misc/files.hpp>
#include <misc/log.hpp>
#include <misc/memory.hpp>
#include <misc/threads.hpp>
#include <move_data_structure/move_data_structure.hpp>
#include <pfp_rlbwt/comp/pfp_rlbwt_comp.hpp>
#include <pfp_rlbwt/comp/run_model.hpp>

namespace pfp_rlbwt::comp {

/**
 * @brief reads a file into memory
 * @param path the file
 * @return its bytes
 */
inline no_init_vector<uint8_t> read_file(const std::string& path)
{
    std::error_code ec;
    const uint64_t size = std::filesystem::file_size(path, ec);
    if (ec) throw std::runtime_error("cannot read " + path);
    no_init_vector<uint8_t> bytes;
    resize_no_init(bytes, size);
    if (!read_file_in_blocks(path, size, reinterpret_cast<char*>(bytes.data())))
        throw std::runtime_error("cannot read " + path);
    return bytes;
}

/**
 * @brief writes bytes that are reported from right to left to consecutive positions of a file with O_DIRECT (see
 *        direct_positioned_file), through a buffer that holds a window of the file starting at a multiple of the sector
 *        size; writers of disjoint ranges of the file can be used by different threads
 */
class backward_writer {
private:
    static constexpr uint64_t sector = direct_io_alignment;

    direct_positioned_file* file; // the file
    uint64_t window_bytes; // the size of the window in bytes (a multiple of sector)
    uint64_t window = 0; // the position in the file of buffer[0] (a multiple of sector)
    uint64_t begin = 0; // [begin, end) are the buffered positions of the file
    uint64_t end = 0;
    char* buffer; // [0..window_bytes-1] the byte at position window + i is buffer[i] (aligned to sector)

public:
    /**
     * @brief creates a writer
     * @param file the file
     * @param buffer_bytes the size of the buffer in bytes (rounded up to a multiple of the sector size)
     */
    backward_writer(direct_positioned_file& file, uint64_t buffer_bytes)
        : file(&file), window_bytes(std::max<uint64_t>(1, div_ceil(buffer_bytes, sector)) * sector),
          buffer(static_cast<char*>(allocate_aligned(window_bytes)))
    {
    }

    backward_writer(backward_writer&& other) noexcept
        : file(other.file), window_bytes(other.window_bytes), window(other.window), begin(other.begin),
          end(other.end), buffer(std::exchange(other.buffer, nullptr))
    {
        other.begin = other.end;
    }

    backward_writer(const backward_writer&) = delete;
    backward_writer& operator=(const backward_writer&) = delete;
    backward_writer& operator=(backward_writer&&) = delete;

    ~backward_writer()
    {
        flush();
        platform::free_aligned(buffer);
    }

    // writes the byte c to the position i (consecutive positions from right to left are buffered)
    inline void put(uint64_t i, uint8_t c)
    {
        if (i + 1 != begin || i < window) {
            flush();
            end = i + 1;
            // the window ends at the end of the sector containing i
            const uint64_t window_end = div_ceil(end, sector) * sector;
            window = window_end > window_bytes ? window_end - window_bytes : 0;
        }

        begin = i;
        buffer[i - window] = char(c);
    }

    // writes the buffered bytes to the file
    void flush()
    {
        if (begin == end) return;
        file->write(buffer + (begin - window), end - begin, begin);
        begin = end;
    }
};

/**
 * @brief the decoded content of a compressed file: n, the alphabet, M_LF with the symbols of the runs and the
 *        start pairs
 */
class decompressor {
private:
    using mlf_t = move_data_structure<POS, UNI>;

    uint64_t n = 0; // n, the length of the text
    std::array<uint8_t, num_byte_values + 1> byte_of {}; // [1..sigma] byte_of[c] = the byte with the code c
    mlf_t M_LF; // M_LF of the BWT of T$; row 0 stores the code of each input interval (0 for $, 1..sigma)
    std::vector<uint64_t> start_rows; // [0..m-1] the BWT rows of the start pairs, sorted by text position
    std::vector<uint64_t> start_positions; // [0..m-1] SA[start_rows[j]]

public:
    /**
     * @brief decodes the run-length encoded BWT and the start pairs at the end of a compressed file and builds M_LF
     * @param in reads the compressed file from the start of the run-length encoded BWT
     * @param params parameters
     */
    decompressor(byte_reader in, const decompress_params& params)
    {
        const uint64_t p = std::max<uint64_t>(1, params.num_threads);
        auto time = now();
        log_phase_start(params.log, time, "decoding the runs");
        n = in.varint();
        if (n == 0) return;
        const uint64_t N = n + 1; // the length of T$

        // the alphabet
        uint64_t sigma = 0;

        for (uint64_t c = 0; c < num_byte_values; c += std::numeric_limits<uint8_t>::digits) {
            const uint8_t bits = in.byte();

            for (uint64_t i = 0; i < std::numeric_limits<uint8_t>::digits; i++) {
                if ((bits >> i) & 1) byte_of[++sigma] = uint8_t(c + i);
            }
        }

        const uint64_t r = in.varint();
        const uint64_t terminator_run = in.varint();
        const uint64_t m = in.varint();
        if (sigma == 0 || r < 2 || r > N || terminator_run >= r || m >= n) throw invalid_file();
        std::vector<uint64_t> start_runs;

        if (m > 0) {
            const uint8_t run_width = std::bit_width(r - 1);
            const uint8_t difference_width = in.byte();
            if (difference_width > std::numeric_limits<uint64_t>::digits) throw invalid_file();
            bit_unpacker pairs(in.bytes(packed_bytes(m, run_width + difference_width)));
            start_runs.resize(m);
            start_positions.resize(m);

            for (uint64_t j = 0; j < m; j++) {
                start_runs[j] = pairs.pop(run_width);
                start_positions[j] = pairs.pop(difference_width) + (j == 0 ? 0 : start_positions[j - 1]);

                if (start_runs[j] >= r || start_positions[j] == 0 || start_positions[j] >= n ||
                    (j > 0 && start_positions[j] == start_positions[j - 1])) {
                    throw invalid_file();
                }
            }
        }

        // the coded run symbols and run lengths of each block
        const uint64_t num_blocks = in.varint();
        if (num_blocks == 0 || num_blocks > r) throw invalid_file();
        std::vector<uint64_t> symbol_sizes(num_blocks), length_sizes(num_blocks);

        for (uint64_t b = 0; b < num_blocks; b++) {
            symbol_sizes[b] = in.varint();
            length_sizes[b] = in.varint();
        }

        std::vector<const uint8_t*> block_data(num_blocks); // [0..num_blocks-1] the coded runs of block b

        for (uint64_t b = 0; b < num_blocks; b++) {
            block_data[b] = in.bytes(symbol_sizes[b]);
            in.bytes(length_sizes[b]);
        }

        if (in.remaining() != 0) throw invalid_file();
        auto block_start = [&](uint64_t b) { return uint64_t((uint128_t(b) * r) / num_blocks); };

        /* the pairs (i, LF(i)) of the run starts i, and the codes of the runs: 8 bits per code if they fit into a
           byte (then they order the pairs by LF(i)), else the bit width of sigma */
        const uint8_t code_width =
            sigma < num_byte_values ? std::numeric_limits<uint8_t>::digits : std::bit_width(sigma);
        mds_pairs I(r, N);
        interleaved_bit_aligned_vectors codes({ code_width });
        codes.resize_no_init(r);
        // [b * (sigma + 1) + c] the number of occurrences of the code c in block b; then the number of occurrences of
        // c in the blocks before b plus C[c]
        std::vector<uint64_t> occurrences(num_blocks * (sigma + 1), 0);
        std::vector<uint64_t> block_size(num_blocks); // [0..num_blocks-1] the number of BWT rows of block b
        const uint64_t threads = std::min(p, num_blocks);

        // decodes the blocks; I stores the length of each run
        std::atomic<uint64_t> next_block = 0;

        run_on_threads(threads, [&](uint64_t, uint64_t) {
            for (uint64_t b = next_block++; b < num_blocks; b = next_block++) {
                range_decoder symbols(block_data[b], symbol_sizes[b]);
                range_decoder lengths(block_data[b] + symbol_sizes[b], length_sizes[b]);
                run_model model(sigma, n);
                uint64_t* occ = occurrences.data() + b * (sigma + 1);
                uint64_t rows = 0;

                for (uint64_t x = block_start(b); x < block_start(b + 1); x++) {
                    uint64_t c = 0, len = 1;

                    if (x == terminator_run) {
                        model.skip_terminator();
                    } else {
                        model.decode(symbols, lengths, c, len);
                        c++;
                    }

                    if (c > sigma || len == 0 || len > N - rows) throw invalid_file();
                    I.set<0>(x, len);
                    codes.set_parallel<0, uint64_t>(x, c);
                    occ[c] += len;
                    rows += len;
                }

                block_size[b] = rows;
            }
        });

        // block_size[b] becomes the first BWT row of block b; LF(i) = C[c] + rank(L, c, i - 1) for the start i of a
        // run of c
        uint64_t rows = 0;

        for (uint64_t b = 0; b < num_blocks; b++) {
            const uint64_t size = block_size[b];
            block_size[b] = rows;
            rows += size;
            if (rows > N) throw invalid_file();
        }

        if (rows != N) throw invalid_file();
        uint64_t smaller = 0;

        for (uint64_t c = 0; c <= sigma; c++) {
            for (uint64_t b = 0; b < num_blocks; b++) {
                const uint64_t occ = occurrences[b * (sigma + 1) + c];
                occurrences[b * (sigma + 1) + c] = smaller;
                smaller += occ;
            }
        }

        next_block = 0;

        run_on_threads(threads, [&](uint64_t, uint64_t) {
            for (uint64_t b = next_block++; b < num_blocks; b = next_block++) {
                uint64_t* occ = occurrences.data() + b * (sigma + 1);
                uint64_t i = block_size[b];

                for (uint64_t x = block_start(b); x < block_start(b + 1); x++) {
                    const uint64_t c = codes.get<0, uint64_t>(x);
                    const uint64_t len = I.get<0>(x);
                    I.set<0>(x, i);
                    I.set<1>(x, occ[c]);
                    occ[c] += len;
                    i += len;
                }
            }
        });

        start_rows.resize(m);
        for (uint64_t j = 0; j < m; j++) start_rows[j] = I.get<0>(start_runs[j]);

        log_phase_end(params.log, time);
        log_phase_start(params.log, time, "building M_LF");
        mds_params mlf_params { .num_threads = p, .log = false };
        if (code_width == std::numeric_limits<uint8_t>::digits) mlf_params.q_order_keys = codes.data<uint8_t>();
        const mds_pair_row code_row { .values = &codes };
        mlf_params.pair_rows = &code_row;
        M_LF = mlf_t(std::move(I), N, mlf_params, { uint8_t(std::bit_width(sigma)) });
        log_phase_end(params.log, time);
    }

    // returns n, the length of the text
    uint64_t size() const { return n; }

    // returns the number of pieces the text is reverted in (the number of start pairs plus one)
    uint64_t num_pieces() const { return start_rows.size() + 1; }

    /**
     * @brief returns piece j in [0..m] of the text: T[s_{j-1}..s_j - 1], where s_j is the text position of the j-th
     *        start pair, s_{-1} = 0 and s_m = n
     * @param j the piece
     * @return (s_{j-1}, s_j)
     */
    std::pair<uint64_t, uint64_t> piece(uint64_t j) const
    {
        const uint64_t m = start_rows.size();
        return { j == 0 ? 0 : start_positions[j - 1], j == m ? n : start_positions[j] };
    }

    // the state of a piece that is being reverted
    struct strand_t {
        uint64_t i = 0; // the current BWT row
        uint64_t x = 0; // the input interval of M_LF containing i
        uint64_t t = 0; // the text position of L[i]
        uint64_t first = 0; // the first text position of the piece
    };

    // returns the state at the start of piece j: at the BWT row of the j-th start pair (for j = m at row 0,
    // L[0] = T[n-1]) and the last text position of the piece
    strand_t start_piece(uint64_t j) const
    {
        const auto [first, end] = piece(j);
        const uint64_t i = j == start_rows.size() ? 0 : start_rows[j];
        return { .i = i, .x = M_LF.interval_and_offset(i).first, .t = end - 1, .first = first };
    }

    /**
     * @brief reverts pieces of the text with up to k strands on the calling thread: per turn, each strand reports the
     *        byte at its text position, moves to the next BWT row with M_LF and prefetches the block that its next
     *        move reads; a strand that has reported its piece takes a new one
     * @param k the number of strands
     * @param take function called as take(s, strand, wait); it assigns a piece to strand s (see start_piece) and
     *        returns true, or returns false if there is none; wait is true iff no other strand has a piece
     * @param report function called as report(s, i, c) with c = T[i], from right to left within each piece
     * @param finish function called as finish(s) after strand s has reported its piece
     */
    template <typename take_fnc_t, typename report_fnc_t, typename finish_fnc_t>
    void revert_strands(uint64_t k, take_fnc_t take, report_fnc_t report, finish_fnc_t finish) const
    {
        std::vector<strand_t> strands(k);
        std::vector<uint64_t> active; // the strands that revert a piece
        std::vector<uint64_t> idle(k); // the strands without a piece
        for (uint64_t s = 0; s < k; s++) idle[s] = k - 1 - s;
        active.reserve(k);

        // assigns pieces to idle strands
        auto refill = [&] {
            while (!idle.empty() && take(idle.back(), strands[idle.back()], active.empty())) {
                active.push_back(idle.back());
                idle.pop_back();
            }
        };

        refill();

        while (!active.empty()) {
            for (uint64_t a = 0; a < active.size();) {
                const uint64_t s = active[a];
                strand_t& st = strands[s];
                report(s, st.t, byte_of[M_LF.row<0>(st.x)]);

                if (st.t != st.first) {
                    st.t--;
                    M_LF.move(st.i, st.x);
                    M_LF.prefetch_move(st.x);
                    a++;
                    continue;
                }

                finish(s);
                active[a] = active.back();
                active.pop_back();
                idle.push_back(s);
                refill();
            }
        }
    }

    /**
     * @brief reports the bytes of the text with up to num_threads threads, each with num_strands strands (see
     *        revert_strands): a strand takes one piece after another
     * @param report function called as report(w, i, c) with c = T[i], by several threads at the same time; w in
     *        [0..num_threads * num_strands - 1] is the strand, and each strand reports its pieces one after another,
     *        each from right to left
     * @param num_threads the maximum number of threads
     * @param num_strands the number of strands per thread
     */
    template <typename report_fnc_t>
    void revert(report_fnc_t report, uint64_t num_threads, uint64_t num_strands) const
    {
        if (n == 0) return;
        const uint64_t pieces = num_pieces();
        const uint64_t k = std::max<uint64_t>(1, num_strands);
        const uint64_t threads = std::max<uint64_t>(1, std::min(num_threads, div_ceil(pieces, k)));
        std::atomic<uint64_t> next = 0; // the next piece to take

        run_on_threads(threads, [&](uint64_t thread, uint64_t) {
            revert_strands(k, [&](uint64_t, strand_t& st, bool) {
                if (next.load(std::memory_order_relaxed) >= pieces) return false;
                const uint64_t j = next.fetch_add(1, std::memory_order_relaxed);
                if (j >= pieces) return false;
                st = start_piece(j);
                return true;
            }, [&](uint64_t s, uint64_t i, uint8_t c) { report(thread * k + s, i, c); }, [](uint64_t) {});
        });
    }

    /**
     * @brief reverts the pieces of the text with up to num_threads threads, each with num_strands strands (see
     *        revert_strands), and passes them to consume in the order of their positions; a strand takes the next
     *        piece only while fewer than two pieces per strand have been taken and not consumed
     * @param consume function called as consume(bytes, len) with each piece, by one thread at a time
     * @param num_threads the maximum number of threads
     * @param num_strands the number of strands per thread
     */
    template <typename consume_fnc_t>
    void revert_in_order(consume_fnc_t consume, uint64_t num_threads, uint64_t num_strands) const
    {
        if (n == 0) return;
        const uint64_t pieces = num_pieces();
        const uint64_t k = std::max<uint64_t>(1, num_strands);
        const uint64_t threads = std::max<uint64_t>(1, std::min(num_threads, div_ceil(pieces, k)));
        const uint64_t max_open = 2 * threads * k; // the largest number of pieces taken and not consumed
        std::mutex lock; // guards the variables below
        std::condition_variable progressed; // notified when a piece has been consumed or a thread has failed
        std::vector<no_init_vector<uint8_t>> reverted(pieces); // [0..pieces-1] the reverted pieces not yet consumed
        std::vector<uint8_t> is_reverted(pieces, 0); // [0..pieces-1] whether piece j has been reverted
        uint64_t next = 0; // the next piece to take
        uint64_t consumed = 0; // the number of consumed pieces
        bool consuming = false; // whether a thread consumes pieces
        bool failed = false; // whether a thread has thrown an exception

        run_on_threads(threads, [&](uint64_t, uint64_t) {
            std::vector<uint64_t> piece_of(k); // [0..k-1] the piece of each strand
            std::vector<uint64_t> first_of(k); // [0..k-1] the first text position of the piece of each strand
            std::vector<no_init_vector<uint8_t>> bytes(k); // [0..k-1] the bytes of the piece of each strand

            try {
                revert_strands(k, [&](uint64_t s, strand_t& st, bool wait) {
                    uint64_t j;

                    {
                        std::unique_lock<std::mutex> guard(lock);
                        auto ready = [&] { return failed || next == pieces || next < consumed + max_open; };

                        if (wait) {
                            progressed.wait(guard, ready);
                        } else if (!ready()) {
                            return false;
                        }

                        if (failed || next == pieces) return false;
                        j = next++;
                    }

                    piece_of[s] = j;
                    st = start_piece(j);
                    first_of[s] = st.first;
                    resize_no_init(bytes[s], st.t + 1 - st.first);
                    return true;
                }, [&](uint64_t s, uint64_t i, uint8_t c) {
                    bytes[s][i - first_of[s]] = c;
                }, [&](uint64_t s) {
                    std::unique_lock<std::mutex> guard(lock);
                    reverted[piece_of[s]] = std::move(bytes[s]);
                    bytes[s] = no_init_vector<uint8_t>();
                    is_reverted[piece_of[s]] = 1;
                    if (consuming) return;
                    consuming = true;

                    // consume the reverted pieces that follow the consumed ones
                    while (!failed && consumed < pieces && is_reverted[consumed]) {
                        no_init_vector<uint8_t> ready = std::move(reverted[consumed]);
                        guard.unlock();
                        consume(reinterpret_cast<const char*>(ready.data()), ready.size());
                        ready = no_init_vector<uint8_t>();
                        guard.lock();
                        consumed++;
                        progressed.notify_all();
                    }

                    consuming = false;
                });
            } catch (...) {
                std::lock_guard<std::mutex> guard(lock);
                failed = true;
                progressed.notify_all();
                throw;
            }
        });
    }
};

/**
 * @brief decompresses a compressed file into an output
 * @param file the bytes of the compressed file
 * @param params parameters
 * @param revert_plain is called as revert_plain(text) in the plain mode; it reverts the text into the output and
 *        returns n
 * @param out the output of the FASTA mode
 * @return the length of the output
 */
template <typename revert_plain_t>
inline uint64_t decompress_into(const no_init_vector<uint8_t>& file, const decompress_params& params,
    revert_plain_t revert_plain, fasta_output& out)
{
    const uint64_t p = std::max<uint64_t>(1, params.num_threads);
    byte_reader in(file.data(), file.size());
    const uint8_t mode = in.byte();
    if (mode != text_mode_plain && mode != text_mode_fasta) throw invalid_file();

    if (mode == text_mode_plain) {
        decompressor text(in, params);
        return revert_plain(text);
    }

    const fasta_section section = fasta_section::read(in);
    decompressor text(in, params);
    auto time = now();
    log_phase_start(params.log, time, "reverting the BWT");
    fasta_restorer restorer(section, out);
    text.revert_in_order([&](const char* bytes, uint64_t len) { restorer.append(bytes, len); }, p, params.num_strands);
    restorer.finish();
    const uint64_t size = out.finish();
    log_phase_end(params.log, time);
    return size;
}

inline uint64_t decompress(const std::string& input_file, const std::string& output_file,
    const decompress_params& params)
{
    const uint64_t p = std::max<uint64_t>(1, params.num_threads);
    const no_init_vector<uint8_t> file = read_file(input_file);
    fasta_output out(output_file);

    return decompress_into(file, params, [&](const decompressor& text) {
        const uint64_t n = text.size();
        auto time = now();
        log_phase_start(params.log, time, "reverting the BWT");
        direct_positioned_file output(output_file, n);

        {
            // [0..p * k - 1] the writer of each strand; the strands of a thread share huge_page_bytes of buffers
            const uint64_t k = std::max<uint64_t>(1, params.num_strands);
            std::vector<backward_writer> writers;
            writers.reserve(p * k);
            for (uint64_t w = 0; w < p * k; w++) writers.emplace_back(output, div_ceil(huge_page_bytes, k));
            text.revert([&](uint64_t w, uint64_t i, uint8_t c) { writers[w].put(i, c); }, p, k);
        }

        if (!output.close_file()) throw std::runtime_error("cannot write " + output_file);
        log_phase_end(params.log, time);
        return n;
    }, out);
}

inline std::string decompress(const std::string& input_file, const decompress_params& params)
{
    const uint64_t p = std::max<uint64_t>(1, params.num_threads);
    std::string result;
    fasta_output out(&result);

    decompress_into(read_file(input_file), params, [&](const decompressor& text) {
        resize_no_init(result, text.size());
        auto time = now();
        log_phase_start(params.log, time, "reverting the BWT");
        text.revert([&](uint64_t, uint64_t i, uint8_t c) { result[i] = char(c); }, p, params.num_strands);
        log_phase_end(params.log, time);
        return text.size();
    }, out);

    return result;
}

} // namespace pfp_rlbwt::comp

#endif // PFP_RLBWT_PFP_RLBWT_COMP_DECOMPRESS_TPP
