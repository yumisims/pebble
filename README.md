# Pebble

**Author:** Yumi Sims

Pebble is a command-line tool and C library for smoothing genomic sequencing coverage tracks. It reads per-base coverage from BED or BedGraph files, applies a sliding trimmed mean, and writes the result as BedGraph.

The algorithm follows the coverage-smoothing approach used in [stepStone](https://github.com/wtsi-hpag/stepStone) (`-denoise 1`): within each sliding window, extreme values are trimmed before averaging. This suppresses spikes from mapping artifacts, repeats, and other outliers while preserving broad coverage trends.

## How it works

With the default settings (`-W 1000`, `-S 100`, `--trim-low 0.40`, `--trim-high 0.40`):

1. Take a **1000 bp window** of per-base coverage (first window: bases 0–999).
2. **Sort** the 1000 values and discard the lowest 400 and highest 400.
3. Compute the **mean of the middle 200** values. This is the smoothed coverage for the window start.
4. Write that value as a BedGraph interval `[window_start, window_start + step)`.
5. **Slide the window forward by 100 bp** (next window: bases 100–1099, output interval 100–200) and repeat.

Example output coordinates for the first three steps on a contig starting at 0:

| Step | Window (0-based bases) | BedGraph interval | Smoothed value |
|------|------------------------|-------------------|----------------|
| 0 | 0–999 | 0–100 | mean of middle 200 in window |
| 1 | 100–1099 | 100–200 | mean of middle 200 in window |
| 2 | 200–1199 | 200–300 | mean of middle 200 in window |

Contigs shorter than the window size are skipped. Smoothed values are **rounded to the nearest integer** in the output BedGraph.

## Build

Requires a C11 compiler.

**Makefile:**

```bash
make
make test
```

**CMake:**

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build
```

The CLI binary is written to `build/pebble`.

On Linux, CMake **statically links** executables by default so they run reliably on HPC systems where conda or other environments override `LD_LIBRARY_PATH`. To build dynamically instead:

```bash
cmake -S . -B build -DPEBBLE_STATIC=OFF
cmake --build build
```

If you see a segmentation fault on a dynamically linked build inside a conda environment, either rebuild with static linking or run with a cleared library path:

```bash
LD_LIBRARY_PATH= ./build/pebble --demo
```

## Usage

```bash
# Built-in demo with mock coverage data
./build/pebble --demo

# Smooth a BedGraph file
./build/pebble -i examples/mock_scaffold.bedgraph -o smoothed.bedgraph

# Smooth genome2cov-style 4-column BED (chrom, start, end, coverage)
./build/pebble -i genome2cov.bed -o smoothed.bedgraph

# Standard 6-column BED (score in column 5)
./build/pebble -i regions.bed -f bed -o smoothed.bedgraph

# Single contig
./build/pebble -i genome2cov.bed -c HAP1_SCAFFOLD_1 -o scaffold.bedgraph

# Short contig (smaller window than default 1000)
./build/pebble -i examples/genome2cov_style.bed -W 50 -S 10 -o smoothed.bedgraph
```

### Options

| Flag | Description | Default |
|------|-------------|---------|
| `-i input` | Input BED or BedGraph file | — |
| `-o output` | Output BedGraph file | stdout |
| `-c chrom` | Process only the named chromosome or contig | all |
| `-f bed\|bedgraph` | Input format | inferred from extension |
| `-W window` | Sliding window size | 1000 |
| `-S step` | Step size | 100 |
| `--trim-low frac` | Fraction trimmed from the low end | 0.40 |
| `--trim-high frac` | Fraction trimmed from the high end | 0.40 |
| `--demo` | Run the built-in mock scaffold demo | — |

If no input file is given, Pebble runs in demo mode.

## Input and output

**Input formats**

| Format | Columns | Notes |
|--------|---------|-------|
| BedGraph | `chrom start end value` | Inferred from `.bedgraph` extension |
| genome2cov BED | `chrom start end coverage` | 4-column; inferred from `.bed` extension |
| BED | `chrom start end name score …` | Score taken from column 5; use `-f bed` |

Pebble builds a dense per-base coverage array for each contig from the input intervals, then applies the sliding trimmed mean.

**Output** — 4-column BedGraph: `chrom start end smoothed_coverage`, where `start`/`end` span one step (`[window_start, window_start + step)`) and `smoothed_coverage` is a rounded integer.

## Library API

The core logic lives in `pebble_core`. The public headers are:

- `include/pebble/pebble.h` — `pebble_process()` and configuration
- `include/pebble/pebble_io.h` — BED/BedGraph read and write helpers

```c
#include "pebble/pebble.h"
#include "pebble/pebble_io.h"

pebble_config_t config = {
    .window_size = 1000,
    .step_size   = 100,
    .trim_low    = 0.40,
    .trim_high   = 0.40,
};

size_t n = pebble_output_count(coverage_len, &config);
double *out = malloc(n * sizeof(double));
size_t out_len;

pebble_process(coverage, coverage_len, &config, out, n, &out_len);
pebble_write_bedgraph_file("out.bedgraph", "chr1", 0, &config, out, out_len);
```

Link against `pebble_core` and include the `include/` directory.

## Examples

The `examples/` directory contains mock scaffold and genome2cov-style test data:

```bash
make example   # runs pebble on examples/mock_scaffold.bedgraph
./build/pebble -i examples/genome2cov_style.bed -W 50 -S 10 -o build/genome2cov.smoothed.bedgraph
```

## Reference

Pebble implements the sliding-window trimmed-mean smoothing used in stepStone's `-denoise 1` mode. If you use Pebble in your work, please cite stepStone:

> Zemin Ning. stepStone: a pipeline for identification of chromothripsis breakpoints and cancer rearrangements. GitHub: https://github.com/wtsi-hpag/stepStone

stepStone's coverage plot command documents `-denoise 1` as averaging data points after filtering high and low values within the window size. Pebble's default parameters (`-W 1000`, `-S 100`, `--trim-low 0.40`, `--trim-high 0.40`) match that trimming behaviour (middle 20% mean of a 1000 bp window). See `tests/test_pebble.c` for validation.

For questions about stepStone, contact Zemin Ning (zn1@sanger.ac.uk).
