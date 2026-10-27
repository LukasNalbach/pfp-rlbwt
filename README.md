# pfp-rlbwt

This repository contains a C++ library that computes the run-length encoded BWT of a text with prefix-free parsing
[1, 2], and a compressor based on it.

- **pfp_rlbwt** (`pfp_rlbwt::build`) computes the runs of the BWT of a text in parallel. The text has symbols of 1 to
  8 bytes and is read from a file, from memory or from a FASTA file. Optionally, it also writes the suffix array values
  at the starts and ends of the runs, the LCP values at the starts of the runs, the BWT, the suffix array and the LCP
  array.
- The runs are written to a runs file in a plain or a compressed format. The readers of both formats (sequential, with
  random access on the disk and in RAM) are part of [string-utils](https://github.com/LukasNalbach/string-utils)
  (`rlbwt/`).
- **pfp-rlbwt-comp** (`pfp_rlbwt_comp`) compresses a text or a FASTA file as its run-length encoded BWT and restores it
  in parallel.

## Contents

- [Build](#build)
- [Command-line tools](#command-line-tools)
- [C++ usage](#c-usage)
- [Tests](#tests)
- [References](#references)
- [License](#license)

## Build

Requirements: a C++20 compiler (GCC 13 or later, Clang 18 or later), CMake 3.21 or later and OpenMP.

```shell
git clone https://github.com/LukasNalbach/pfp-rlbwt.git
cd pfp-rlbwt
cmake -B build
cmake --build build -j
```

`cmake` downloads the submodules [string-utils](https://github.com/LukasNalbach/string-utils) and
[move-data-structure](https://github.com/LukasNalbach/move-data-structure). The command-line tools are written to
`build/cli/` and the tests to `build/tests/`.

| Option | Description |
| --- | --- |
| `PFP_RLBWT_BUILD_CLI` | build `pfp-rlbwt` |
| `PFP_RLBWT_BUILD_COMP`, `PFP_RLBWT_BUILD_COMP_CLI` | build the compressor library `pfp_rlbwt_comp` and `pfp-rlbwt-comp` |
| `PFP_RLBWT_BUILD_TESTS` | build the tests |

These options are `ON` if pfp-rlbwt is the top-level project. The general options and environment variables
(`PFP_RLBWT_MARCH_NATIVE`, `PFP_RLBWT_TMPDIR`, ...) are described in string-utils.

## Command-line tools

`<tool> -h` prints the options of a tool.

- **pfp-rlbwt**: writes the runs of the BWT of a file to `<input_file>.bwt-runs`; `-format compressed` selects the
  compressed format, `-sa` and `-lcp` add the SA and LCP samples, `-full-sa` also writes the suffix array.
  `pfp-rlbwt -p 8 -format compressed -sa -lcp text.txt`

  The file may contain every byte value. If a byte below 2 (`pfp_rlbwt::min_symbol`) occurs, the distinct bytes are
  mapped to 2, 3, ... in increasing order, to 16-bit symbols if more than 254 distinct bytes occur. The BWT then
  consists of the mapped bytes, and the runs file stores the map as the section `section_byte_map`: the symbol of the
  smallest byte (64 bits), then the distinct bytes in increasing order (one byte each); `rlbwt::byte_map::load` reads
  it. With `-fasta on` (or `auto`, the default, for files that start like FASTA files), the text is the one of the
  FASTA file: each header line and each run of symbols outside the alphabet of `-A` (default `ACGT`, after
  upper-casing) is replaced with one separator symbol (`0x7F`), the lines before the first header are skipped and
  carriage returns are dropped.
- **pfp-rlbwt-comp**: compresses a file into `<input_file>.rlbwt` and restores it with `-d`. A FASTA file is stored as
  its sequences, its header lines and its line layout (`-fasta on|off|auto`).
  `pfp-rlbwt-comp -p 16 text.txt` and `pfp-rlbwt-comp -d -p 16 -o restored.txt text.txt.rlbwt`

## C++ usage

Add pfp-rlbwt to a CMake project and link the target `pfp_rlbwt` (or `pfp_rlbwt_comp` for the compressor, with
`PFP_RLBWT_BUILD_COMP` set to `ON`):

```cmake
add_subdirectory(pfp-rlbwt/ EXCLUDE_FROM_ALL)
target_link_libraries(<target> PRIVATE pfp_rlbwt)
```

The BWT runs of a text in a file, with the SA and LCP samples:

```cpp
#include <fstream>
#include <iostream>
#include <string>
#include <pfp_rlbwt/pfp_rlbwt.hpp>
#include <rlbwt/readers.hpp>

int main()
{
    // the text abaababaab in a file, read with 2 threads
    std::ofstream("text.txt") << "abaababaab";
    pfp_rlbwt::file_source<uint8_t> text("text.txt");
    pfp_rlbwt::params prm;
    prm.num_threads = 2;
    prm.write_sa_samples = prm.write_lcp_samples = true;
    pfp_rlbwt::result res = pfp_rlbwt::build(text, "text.txt", prm);
    std::cout << res.n_with_terminator << " " << res.r << std::endl; // 11 5

    // the runs of the BWT of T$ = abaababaab$, with $ = 0
    pfp_rlbwt::run_data_file data;
    data.open(std::string("text.txt") + pfp_rlbwt::ext_runs);
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

The compressor:

```cpp
#include <iostream>
#include <string>
#include <pfp_rlbwt/comp/pfp_rlbwt_comp.hpp>

int main()
{
    std::string T = "abaababaab";
    pfp_rlbwt::comp::byte_text text((const uint8_t*) T.data(), T.size());

    // compress T from memory with 2 threads, restore it with 4 threads
    pfp_rlbwt::comp::compress(text, "text.rlbwt", { .num_threads = 2 });
    std::string restored =
        pfp_rlbwt::comp::decompress("text.rlbwt", { .num_threads = 4 });
    std::cout << (restored == T) << std::endl; // 1
}
```

## Tests

```shell
build/tests/test-pfp-rlbwt
build/tests/test-pfp-rlbwt-wide
build/tests/test-pfp-rlbwt-comp
```

The tests check random inputs for at most `PFP_RLBWT_TEST_SECONDS` seconds (default 10) or `PFP_RLBWT_TEST_ITERATIONS`
inputs (default 1000) per function. A failed test prints its seed; `PFP_RLBWT_TEST_SEED=<seed>` repeats its inputs.

## References

[1] Christina Boucher, Travis Gagie, Alan Kuhnle, Ben Langmead, Giovanni Manzini and Taher Mun. Prefix-free parsing for
building big BWTs. Algorithms for Molecular Biology 14, 2019.

[2] Marco Oliva, Travis Gagie and Christina Boucher. Recursive prefix-free parsing for building big BWTs. In Data
Compression Conference (DCC), 2023.

## License

This project is licensed under the MIT License, see [LICENSE](LICENSE). The submodules in [external/](external/) have
their own licenses.
