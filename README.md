# pfp-rlbwt

This repository contains a C++ library that computes the run-length encoded BWT of a text with prefix-free parsing [1], optionally together with the suffix array values at the starts and at the ends of the BWT runs and the LCP values at the starts of the runs, and writes it to a file in one of two formats: plain (one record per run with byte-aligned fields) or compressed (each value readable with O(1) disk reads). It provides streamed and random-access readers of both formats on the disk and two representations in RAM. Optionally, it also writes the BWT with one symbol per position, the suffix array and the LCP array. It also contains a compressor (pfp-rlbwt-comp) that stores a byte text as its run-length encoded BWT.

## Build Instructions

```shell
git clone git@github.com:LukasNalbach/pfp-rlbwt.git
cd pfp-rlbwt
cmake -B build
cmake --build build -j
```

CMake fetches the required submodule [string-utils](https://github.com/LukasNalbach/string-utils) (and, for the compressor, [move-data-structure](https://github.com/LukasNalbach/move-data-structure), see below). The command-line tools are built into `build/cli/` and the tests into `build/tests/`; `PFP_RLBWT_BUILD_CLI` and `PFP_RLBWT_BUILD_TESTS` turn them on or off (both are `ON` if pfp-rlbwt is the top-level project). The general build options (`PFP_RLBWT_MARCH_NATIVE`, `PFP_RLBWT_ENABLE_LTO`, ...) and environment variables (`PFP_RLBWT_TMPDIR`, ...) are described in [string-utils](https://github.com/LukasNalbach/string-utils).

## CLI Tool

**pfp-rlbwt** computes the output of `build` (see below) for a file and writes it to files with the prefix `<input_file>`; `-format compressed` writes the compressed format instead of the plain one, `-sa` and `-lcp` add the SA and LCP samples to the runs file, and it prints its options with `-h`.

```shell
build/cli/pfp-rlbwt -p 8 text.txt              # text.txt.bwt-runs (plain)
build/cli/pfp-rlbwt -p 8 -format compressed -sa -lcp -full-sa text.txt
                     # text.txt.bwt-runs (compressed) with SA and LCP
                     # samples, and text.txt.sa
```

## Usage in C++

```cmake
add_subdirectory(external/pfp-rlbwt EXCLUDE_FROM_ALL)
target_link_libraries(<target> PRIVATE pfp_rlbwt)
```

A text T consists of symbols of 8, 16, 24, 32, 40 or 48 bits (`uint8_t`, `uint16_t`, `pfp_rlbwt::sym24_t`, `uint32_t`, `pfp_rlbwt::sym40_t`, `pfp_rlbwt::sym48_t`); every symbol has to be at least `pfp_rlbwt::min_symbol` (2). `pfp_rlbwt::build` reads it from a `text_source`: a `file_source<sym_t>` (a file with `sizeof(sym_t)` bytes per symbol, little endian, optionally reversed and, for one-byte symbols, with its bytes mapped), a `memory_source<sym_t>` (optionally reversed) or a class derived from `text_source<sym_t>`. It writes the BWT of T$, where $ is the symbol 0, to files with a given prefix; n = |T$| below.

| File | Content |
| --- | --- |
| `<prefix>.bwt-runs` | the runs: its sections (see below), then the table of contents: for each section its identifier, offset and size in bytes (`pfp_rlbwt::run_data_entry`, three 64-bit integers), and the number of sections (a 64-bit integer) |
| `<prefix>.bwt` | the BWT with one symbol per position (`sizeof(sym_t)` bytes each), if `write_bwt` is set |
| `<prefix>.sa` | the suffix array (`byte_width(n)` bytes per value, little endian), if `write_sa` is set |
| `<prefix>.lcp` | the LCP array (`byte_width(n)` bytes per value, little endian), if `write_lcp` is set |

`byte_width(x)` is the fewest whole bytes that hold x. All values are stored little endian, bit-packed data least significant bit of a byte first. Which sections a runs file stores tells its format. Both formats store `section_fields`: the numbers that describe the file (`pfp_rlbwt::run_data_fields`, 64-bit integers): n, the number of runs r, the largest run symbol, the largest run length, whether SA and LCP samples are stored, the number of runs per block, and the widths of the fields of the format.

### The Plain Format (default)

| Section | Content |
| --- | --- |
| `section_records` | one record per run, run i at offset i * (record size): the BWT position of its start (`byte_width(n)` bytes), its symbol (`sizeof(sym_t)` bytes) and, if `write_sa_samples` is set, the SA values at its start and at its end and, if `write_lcp_samples` is set, the LCP value at its start (`byte_width(n)` bytes each) |
| `section_record_samples` | the start of every `runs_per_block`-th run (`byte_width(n)` bytes each), so that the run containing a BWT position is found with one disk read |
| `section_fields` | see above |
| from `first_extra_section` on | the sections of `extra_sections`, each written by its own function |

The records start at offset 0 of the file, and their layout depends only on n, the symbol type and the samples, so other programs can read the file directly (`pfp_rlbwt::record_layout` gives the offsets of the fields).

### The Compressed Format

| Section | Content |
| --- | --- |
| `section_sa_run_starts`, `section_sa_run_ends` | the SA values at the starts and at the ends of the runs, bit-packed with `bit_width(n)` bits each, if `write_sa_samples` is set |
| `section_run_blocks` | the runs in blocks of `runs_per_block` runs: the starts relative to the first start of the block Elias-Fano coded, then the symbols coded with a canonical Huffman code, with the offset of every `runs_per_symbol_sample`-th code (see `pfp_rlbwt::run_block`) |
| `section_lcp_blocks` | the LCP values at the starts of the runs in blocks of `runs_per_block` values, coded with RLE-GEF, the run-length variant of Generalized Elias-Fano [2]: the values minus their minimum are split into low parts (the `split` low bits) and high parts, and the high part is stored once per run of equal high parts; `split` minimizes the size of the block (see `pfp_rlbwt::lcp_block`); if `write_lcp_samples` is set |
| `section_run_block_index` | for each block of runs its offset in `section_run_blocks` and the start of its first run, bit-packed |
| `section_lcp_block_index` | for each LCP block its offset in `section_lcp_blocks`, bit-packed |
| `section_symbol_code` | the canonical Huffman code of the run symbols: the largest code length L, the number of codes of each length 1..L (64 bits each) and the symbols in the order of their codes (`sizeof(sym_t)` bytes each) |
| `section_fields`, from `first_extra_section` on | see above |

A reader keeps the block indexes and the code in RAM (about (2 log2 n) / `runs_per_block` bits per run, plus the code). Each block starts at a byte, so a random access reads one range of the file: the run i, its length and the run containing a BWT position read one block of runs, an SA value reads its bytes, and an LCP value reads one LCP block; each is one disk read (`pread` with `O_DIRECT`).

The coding of the LCP values follows GEF [2] (RLE_GEF; licensed under the Apache License 2.0), but pfp-rlbwt contains none of its code: GEF stores a sequence in partitions of 32000 values whose serialized layout (sdsl-lite vectors and rank directories without an index of the partition offsets) cannot be read with O(1) disk reads per value, so the coding is implemented anew in blocks of `runs_per_block` values that are read with one disk read (the rank within a block is computed with popcounts, so no rank directory is stored). The Apache License 2.0 allows such use from MIT-licensed code; including GEF's code would require keeping its license and notices for those files.

While it runs, `build` stores the suffix array of the dictionary temporarily in `<prefix>.dsa`, the ends of the phrase runs and the symbols before their next phrases in `<prefix>.run-ends` and `<prefix>.run-symbols`, and the inverse suffix array of the parse in `<prefix>.isa`, so that they are not held in RAM while they are not needed; with the compressed format, it writes the records of the plain format to `<prefix>.bwt-runs.records` and the LCP blocks to `<prefix>.bwt-runs.<identifier>` and computes the compressed format from them with `num_threads` threads.

### Readers (`#include <pfp_rlbwt/output/readers.hpp>`)

`run_data_file` opens a runs file: its table of contents, its fields and the parts readers keep in RAM. All readers return a run as `run_record` (start, length, symbol, the SA values at its start and at its end, the LCP value at its start; 0 for values the file does not store) or `run_head` (start and symbol).

| Reader | Format | Access |
| --- | --- | --- |
| `run_stream` | both | the runs [first, end) in order, from any run on, with `O_DIRECT` reads of 1 MiB; one stream per thread (`plain_run_stream`, `compressed_run_stream`) |
| `plain_disk_reader`, `compressed_disk_reader` | plain, compressed | random access on the disk: `run(i)`, `length(i)`, `run_of(j)` (the run containing BWT position j), `sa_start(i)`, `sa_end(i)`, `lcp(i)` with one disk read each; `record(i)` with one disk read (plain) or one per stored part (compressed); one reader per thread |
| `packed_run_data` | both | in RAM, bit-packed: per run its start, its symbol and its samples, interleaved, `bit_width(n)` bits per position and `bit_width(largest symbol)` bits per symbol; loaded with several threads |
| `compressed_run_data` | compressed | in RAM, in the compressed format; loaded with several threads |

| Parameter | Default | Description |
| --- | --- | --- |
| `window` | 10 | length w of a trigger string in symbols |
| `modulus` | 100 | a window is a trigger string with probability 1 / `modulus` |
| `num_threads` | 1 | number of threads |
| `position_bytes` | 0 | bytes per text position in the arrays of the construction, 4 to 6 (0: the fewest that hold n + `window`) |
| `format` | `plain` | the format of `<prefix>.bwt-runs` (`file_format::plain` or `file_format::compressed`) |
| `runs_per_block` | 1024 | runs per block (compressed) or per record sample (plain); a multiple of 8 |
| `runs_per_symbol_sample` | 64 | runs per sample of the symbol codes in a block (compressed), 1 to `runs_per_block` |
| `write_sa_samples`, `write_lcp_samples` | false | store the SA values at the starts and at the ends of the runs, the LCP values at the starts of the runs |
| `write_bwt`, `write_sa`, `write_lcp` | false | write `<prefix>.bwt`, `<prefix>.sa`, `<prefix>.lcp` |
| `extra_sections` | empty | further sections of `<prefix>.bwt-runs`: an identifier (at least `first_extra_section`, distinct) and a function that writes the bytes to a `std::ostream` |
| `log` | false | print the phases and their running times |

```cpp
#include <fstream>
#include <iostream>
#include <string>
#include <pfp_rlbwt/output/readers.hpp>
#include <pfp_rlbwt/pfp_rlbwt.hpp>

int main()
{
    // the text abaababaab in a file, read with 2 threads; additionally write
    // the SA and LCP samples, a section of our own and the suffix array
    std::ofstream("text.txt") << "abaababaab";
    pfp_rlbwt::file_source<uint8_t> text("text.txt");
    pfp_rlbwt::params prm;
    prm.num_threads = 2;
    prm.write_sa_samples = prm.write_lcp_samples = true;
    prm.write_sa = true;
    prm.extra_sections.push_back({ pfp_rlbwt::first_extra_section,
        [](std::ostream& out) { out << "a note"; } });
    pfp_rlbwt::result res = pfp_rlbwt::build(text, "text.txt", prm);
    std::cout << res.n_with_terminator << " " << res.r << std::endl; // 11 5

    // text.txt.bwt-runs in the plain format: records of 5 bytes (start,
    // symbol, SA value at the start and at the end, LCP value)
    pfp_rlbwt::run_data_file data;
    data.open(std::string("text.txt") + pfp_rlbwt::ext_runs);
    std::cout << data.sections().size() << " sections, the note has "
              << data.section(pfp_rlbwt::first_extra_section).size
              << " bytes" << std::endl; // 4 sections, the note has 6 bytes

    // the runs of the BWT of T$ = abaababaab$, with $ = 0
    pfp_rlbwt::run_stream runs(data);
    pfp_rlbwt::run_record run;
    while (!runs.done()) {
        runs.next(run);
        std::cout << (run.symbol == pfp_rlbwt::terminator ? '$'
            : char(run.symbol)) << run.length << " (" << run.sa_start << ", "
            << run.sa_end << ", " << run.lcp << ") ";
    }
    std::cout << std::endl;
    // b3 (10, 2, 0) a1 (8, 8, 1) b1 (5, 5, 2) $1 (0, 0, 5) a5 (3, 4, 3)

    // random access on the disk and in RAM
    pfp_rlbwt::plain_disk_reader disk(data);
    std::cout << disk.run_of(7) << " " << disk.length(4) << std::endl; // 4 5
    pfp_rlbwt::packed_run_data packed(data, 2);
    std::cout << packed.run(1).start << " " << packed.sa_end(0) << std::endl; // 3 2
}
```

The same text in the compressed format, read in RAM in this format:

```cpp
#include <iostream>
#include <string>
#include <pfp_rlbwt/output/readers.hpp>
#include <pfp_rlbwt/pfp_rlbwt.hpp>

int main()
{
    std::string T = "abaababaab";
    pfp_rlbwt::params prm;
    prm.format = pfp_rlbwt::file_format::compressed;
    prm.write_lcp_samples = true;
    pfp_rlbwt::build(pfp_rlbwt::memory_source<uint8_t>(
        (const uint8_t*) T.data(), T.size()), "text", prm);

    pfp_rlbwt::run_data_file data;
    data.open(std::string("text") + pfp_rlbwt::ext_runs);
    pfp_rlbwt::compressed_run_data runs(data);
    for (uint64_t i = 0; i < runs.num_runs(); i++)
        std::cout << runs.run(i).start << ":" << runs.lcp(i) << " ";
    std::cout << std::endl; // 0:0 3:1 4:2 5:5 6:3
}
```

A text in memory, here with 32-bit symbols, and its reverse:

```cpp
#include <iostream>
#include <vector>
#include <pfp_rlbwt/pfp_rlbwt.hpp>

int main()
{
    std::vector<uint32_t> T =
        { 1000, 2000, 1000, 1000, 2000, 1000, 2000, 1000, 1000, 2000 };
    pfp_rlbwt::params prm; // only ints.bwt-runs with the runs
    pfp_rlbwt::result res = pfp_rlbwt::build(
        pfp_rlbwt::memory_source<uint32_t>(T.data(), T.size()), "ints", prm);
    std::cout << res.r << std::endl; // 5

    // the reversed text: the BWT of
    // 2000 1000 1000 2000 1000 2000 1000 1000 2000 1000 $
    pfp_rlbwt::result rev = pfp_rlbwt::build(
        pfp_rlbwt::memory_source<uint32_t>(T.data(), T.size(), true), "ints-rev",
        prm);
    std::cout << rev.r << std::endl; // 7
}
```

### FASTA Files

`pfp_rlbwt::fasta_source<sym_t>` (`#include <pfp_rlbwt/fasta.hpp>`) reads the text of a FASTA file directly from the file (or from the file in memory), so that the text never has to be written. The transformation, the scan and the reader of the text come from [string-utils](https://github.com/LukasNalbach/string-utils) (`fasta/scan.hpp`); `pfp_rlbwt/fasta.hpp` makes them available in the namespace `pfp_rlbwt`. A `fasta_transform` maps the bytes of the file to the text: a header line (starting with `>`) yields one separator, a body byte in the alphabet yields itself (upper-cased if `to_upper`), a run of body bytes that are not in the alphabet yields one separator, and line ends yield nothing. Its options decide whether the header of the first record yields a separator (`separator_before_first_record`), whether the lines before the first header form a record (`lines_before_first_header`) and whether a carriage return is dropped, dropped only before a line feed, or kept (`carriage_return`).

`scan_fasta_lines` scans the file once in parallel and returns the length of the text, the frequencies of its bytes and samples (a text position with the offset of the byte in the file that yields it and the state of the transformation there, one every `fasta_sample_period` positions); a reader of `fasta_source` (`fasta_text_reader`) starts at the last sample before the position it reads. Each thread scans chunks of the file at arbitrary offsets, and a chunk type derived from `fasta_chunk` can collect the header lines and the lengths of the sequence lines (`header_start`, `header_bytes`, `sequence_symbols`, `line_end`), which the caller receives chunk by chunk in the order of the chunks (`merge`).

```cpp
#include <array>
#include <fstream>
#include <iostream>
#include <pfp_rlbwt/fasta.hpp>
#include <pfp_rlbwt/pfp_rlbwt.hpp>

int main()
{
    std::ofstream("genomes.fa") << ">a\nACGTNNAC\ngt\n>b\nACGTAC\n";

    // keeps A, C, G and T (upper-cased); a header or a run of other bytes
    // yields the separator #
    pfp_rlbwt::fasta_transform transform("ACGT", '#');
    pfp_rlbwt::fasta_scan scan = pfp_rlbwt::scan_fasta_lines<pfp_rlbwt::fasta_chunk>(
        [] { return direct_reader("genomes.fa"); }, transform,
        { .num_threads = 2 }, [](auto&&...) { });
    // the text #ACGT#ACGT#ACGTAC
    std::cout << scan.text_size << std::endl; // 17

    // the alphabet # < A < C < G < T is mapped to 2, 3, 4, 5, 6
    std::array<uint8_t, 256> map {};
    uint8_t symbol = pfp_rlbwt::min_symbol;
    for (uint64_t c = 0; c < map.size(); c++)
        if (scan.freq[c] > 0) map[c] = symbol++;

    pfp_rlbwt::fasta_source<uint8_t> text("genomes.fa", scan.text_size,
        transform, scan.samples.data(), scan.samples.size(), false, &map);
    pfp_rlbwt::params prm;
    prm.num_threads = 2;
    pfp_rlbwt::result res = pfp_rlbwt::build(text, "genomes.fa", prm);
    std::cout << res.r << std::endl; // 7
}
```

## Compressor

**pfp-rlbwt-comp** stores a byte text as the run-length encoded BWT [3] of the text, with the run symbols and run lengths compressed by an arithmetic coder with context mixing, and restores the text in parallel with the move data structure [4, 5] for LF. The BWT is computed with `build` (the plain format with SA samples). The compressor needs the submodule [move-data-structure](https://github.com/LukasNalbach/move-data-structure), which CMake fetches only if `PFP_RLBWT_BUILD_COMP` (the library `pfp_rlbwt_comp`) or `PFP_RLBWT_BUILD_COMP_CLI` (the command-line tool) is `ON`; both are `ON` if pfp-rlbwt is the top-level project, and the tests of the compressor (`test-pfp-rlbwt-comp`) are built with the other tests if the library is.

**pfp-rlbwt-comp** compresses a file into `<input_file>.rlbwt` (`-fasta on|off|auto`: FASTA mode, see below; `auto`, the default, uses it for files that start with `>`). With `-d`, it restores the file into `<input_file>` without `.rlbwt`. `-t` sets the number of threads, `-o` the output file; `-h` prints all options.

```shell
# writes text.txt.rlbwt
build/cli/pfp-rlbwt-comp -t 16 text.txt
# restores text.txt into restored.txt
build/cli/pfp-rlbwt-comp -d -t 16 -o restored.txt text.txt.rlbwt
```

### FASTA Mode

With `-fasta on` (`comp::compress_params::fasta`), pfp-rlbwt-comp splits a FASTA file into its sequences, its header lines and its line layout. The sequences of all records, separated by a line feed, form the text whose BWT is stored; the headers and the line layout are coded separately with context models. Any input is restored byte for byte, also files that are not FASTA. `pfp-rlbwt-comp -d` detects the mode from the file.

```shell
# writes genomes.fa.rlbwt
build/cli/pfp-rlbwt-comp -fasta on -t 16 genomes.fa
# restores genomes.fa into restored.fa
build/cli/pfp-rlbwt-comp -d -t 16 -o restored.fa genomes.fa.rlbwt
```

- A line starting with `>` is a header line, every other line is a sequence line of the current record; lines before the first header line form a record without a header.
- The header lines are split into tokens (runs of digits, runs of letters, single other bytes); each token is coded relative to the token at the same position of the previous header: equal, a number with a coded difference, or new.
- The line layout stores per record the line width, the lines that deviate from it and the empty lines at its end; the line end (`\n` or `\r\n`) and a missing line end at the end of the file are stored once.
- Lowercase letters, `N` runs and IUPAC codes stay in the sequences.

The split of a FASTA file into its parts, the header and line layout models, the restoration of the file, the range coder and the context models are the ones of [string-utils](https://github.com/LukasNalbach/string-utils) (`fasta/`, `coding/`).

The compression scans the FASTA file once in parallel: it codes the headers and the line layout and samples the file offsets of the sequence text; the construction of the BWT then reads the sequence text directly from the FASTA file (`pfp_rlbwt::fasta_source`), so the sequences are not written anywhere. The decompression passes the reverted pieces of the sequence text in order to the restoration of the FASTA file, which writes them with the headers and the line ends; it keeps at most two pieces per strand (see below) in memory.

### Usage of the Compressor in C++

```cmake
set(PFP_RLBWT_BUILD_COMP ON)
add_subdirectory(external/pfp-rlbwt EXCLUDE_FROM_ALL)
target_link_libraries(<target> PRIVATE pfp_rlbwt_comp)
```

`comp::compress` reads a text from a file or from memory (`comp::byte_text`) and writes the compressed file; it stores the temporary files of `build` next to the output file. `comp::decompress` restores the text into a file or into a `std::string`. The parameters are the number of threads, the maximum number of start pairs k (`max_start_pairs`, default 4096, see below), the largest number of runs per block (`runs_per_block`, default 2^20, see below), the number of strands per thread of the decompression (`num_strands`, default 16, see below) and whether to log the phases.

```cpp
#include <iostream>
#include <string>
#include <pfp_rlbwt/comp/pfp_rlbwt_comp.hpp>

int main()
{
    // T = abaababaab; the BWT of T$ is bbbab$aaaaa (5 runs)
    std::string T = "abaababaab";
    pfp_rlbwt::comp::byte_text text((const uint8_t*) T.data(), T.size());

    // compress T from memory with 2 threads
    pfp_rlbwt::comp::compress_result res =
        pfp_rlbwt::comp::compress(text, "text.rlbwt", { .num_threads = 2 });
    std::cout << res.n << " " << res.sigma << " " << res.r << " "
              << res.num_start_pairs << std::endl; // 10 2 5 3
    std::cout << res.file_bytes << " bytes" << std::endl;

    // restore T into memory with 4 threads, and into a file
    std::string restored =
        pfp_rlbwt::comp::decompress("text.rlbwt", { .num_threads = 4 });
    std::cout << (restored == T) << std::endl; // 1
    pfp_rlbwt::comp::decompress("text.rlbwt", "text.txt");

    // compress a file, with the phases logged
    pfp_rlbwt::comp::compress("text.txt", "text.txt.rlbwt",
        { .log = true });
}
```

### File Format of the Compressor

A variable-length integer stores 7 bits per byte, the lowest group first, and sets the highest bit of each byte except the last. The file starts with the mode: one byte, 0 for a text, 1 for a FASTA file. In the FASTA mode, the following fields follow:

1. the flags (one byte): bit 0 = the lines end with `\r\n`, bit 1 = the last line has no line end, bit 2 = the file does not start with a header line;
2. the number of records, the size of the coded headers and the size of the coded line layout in bytes (variable-length integers);
3. the coded headers and the coded line layout.

Then follows the BWT of the text T (the sequence text in the FASTA mode). Let n be the length of T, sigma the number of its distinct bytes and r the number of runs of the BWT of T$, where $ is smaller than every byte:

1. n (variable-length integer); for n = 0, the file ends here;
2. the alphabet: 32 bytes, bit c mod 8 of byte c div 8 is set iff the byte c occurs in T; the symbol of the i-th smallest occurring byte is i;
3. r, the index of the run of $ and the number m of start pairs (variable-length integers);
4. if m > 0: the bit width w of the largest difference below, then m start pairs, sorted by text position: the index of a run (bit_width(r - 1) bits) and the difference of the suffix array value s at its start to the one of the previous start pair (to 0 for the first one, w bits), bit-packed lowest bit first and padded to a byte;
5. the number B of blocks of runs (variable-length integer); block b contains the runs with the indices floor(b r / B) to floor((b + 1) r / B) - 1;
6. for each block, the size of its coded run symbols and the size of its coded run lengths in bytes (variable-length integers);
7. for each block, its coded run symbols and its coded run lengths.

The start pairs are chosen as follows: [0, n) is split into k = min(r, `max_start_pairs`) ranges of equal length, and for each range that contains such a value, the run start i with the smallest SA[i] in [1, n - 1] in the range is a start pair. The decompression walks LF from each start pair (and from row 0 for the end of T) to the text position of the previous start pair. Each thread reverts `num_strands` pieces alternately: per turn, a strand reports one byte, moves to the next BWT row with M_LF and prefetches the block of M_LF that its next move reads; a strand that has finished its piece takes the next one. In the plain mode, each strand writes its pieces directly to their positions in the output.

The runs are split into B = ceil(r / `runs_per_block`) blocks (see above), which are coded and decoded in parallel, each with new models. In each block, the run symbols and the run lengths (except for the run of $, of length 1) are coded in two streams with a binary arithmetic coder (32-bit range, carry propagation; the zero bytes at the end of a stream are not stored). A value is coded bit by bit, the highest bit first, along a binary tree of adaptive bit models (the probability of a 0-bit, moved towards each coded bit by 1 / (c + 2) after c updates, and by 1/257 after 255 updates):

- the symbol of a run with three models (contexts: none, the symbol of the previous run, the symbols of the previous two runs), whose predictions are mixed in the logistic domain with one weight vector per tree node;
- the bit width b of the length of a run with three mixed models (contexts: none, the symbol of the run, the bit widths of the lengths of the previous run and of the last run with the same symbol);
- the 8 bits below the leading one of the length with one model per b, and the remaining bits with probability 1/2.

## Tests

```shell
build/tests/test-pfp-rlbwt
build/tests/test-pfp-rlbwt-wide    # the same tests with 40- and 64-bit
                                 # indices for small inputs
build/tests/test-pfp-rlbwt-comp    # the compressor
```

The tests build random inputs in both formats and compare the runs read by every reader with those computed from the suffix array; the tests of the compressor compress random and repetitive texts and FASTA files and restore them with 1, 2, 3 and 8 threads. They check random inputs for at most `PFP_RLBWT_TEST_SECONDS` seconds (default 10) or `PFP_RLBWT_TEST_ITERATIONS` inputs (default 1000) per functionality; `PFP_RLBWT_TEST_SEED=<seed>` replays the inputs of a seed that a failed test printed.

## References

[1] Christina Boucher, Travis Gagie, Alan Kuhnle, Ben Langmead, Giovanni Manzini and Taher Mun. Prefix-free parsing for building big BWTs.
Algorithms for Molecular Biology 14, 2019.

[2] Michelangelo Pucci. Generalized Elias-Fano (GEF). https://github.com/mxpucci/generalized-elias-fano, Apache License 2.0, 2025.

[3] Michael Burrows and David J. Wheeler. A block-sorting lossless data compression algorithm. Technical Report 124, Digital Equipment Corporation, 1994.

[4] Takaaki Nishimoto and Yasuo Tabei. Optimal-time queries on bwt-runs compressed indexes. In 48th International Colloquium on Automata, Languages, and Programming (ICALP), 2021. ([paper](https://arxiv.org/abs/2006.05104))

[5] Nico Bertram, Johannes Fischer and Lukas Nalbach. Move-r: Optimizing the r-index. In Symposium on Experimental Algorithms (SEA), 2024. ([paper](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.SEA.2024.1))

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
