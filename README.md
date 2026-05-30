# Pebble

**Author:** Yumi Sims

Pebble is a command-line tool and C library for smoothing genomic sequencing coverage tracks. It reads per-base coverage from BED or BedGraph files, applies a sliding trimmed mean, and writes the result as BedGraph.

The algorithm follows the coverage-smoothing approach used in [stepStone](https://github.com/wtsi-hpag/stepStone) (`-denoise 1`): within each sliding window, extreme values are trimmed before averaging. This suppresses spikes from mapping artifacts, repeats, and other outliers while preserving broad coverage trends.

## How it works

For each position along a chromosome, Pebble:

1. Collects coverage values in a sliding window (default: 1000 bases).
2. Sorts the window and trims the lowest and highest fractions (default: 40% from each end).
3. Computes the mean of the remaining values.
4. Advances by a step size (default: 100 bases) and repeats.

This produces a smoothed track that is less sensitive to local noise than a plain sliding average.

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

## Usage

```bash
# Built-in demo with mock coverage data
./build/pebble --demo

# Smooth a BedGraph file
./build/pebble -i examples/mock_scaffold.bedgraph -o smoothed.bedgraph

# Smooth a BED file, single chromosome
./build/pebble -i regions.bed -f bed -c chr1 -o smoothed.bedgraph
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

**Input** — BedGraph (4-column: chrom, start, end, value) or BED (coverage encoded as interval width). Format is inferred from the file extension (`.bedgraph` or `.bed`) unless overridden with `-f`.

**Output** — BedGraph with smoothed values at each step interval.

## Library API

The core logic lives in `libpebble` (`pebble_core`). The public headers are:

- `include/pebble/pebble.h` — `pebble_process()` and configuration
- `include/pebble/pebble_io.h` — BED/BedGraph read and write helpers

```c
#include "pebble/pebble.h"

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
```

Link against `pebble_core` and include the `include/` directory.

## Examples

The `examples/` directory contains mock scaffold data:

```bash
make example   # runs pebble on examples/mock_scaffold.bedgraph
```

## Reference

Pebble implements the sliding-window trimmed-mean smoothing used in stepStone's `-denoise 1` mode. If you use Pebble in your work, please cite stepStone:

> Zemin Ning. stepStone: a pipeline for identification of chromothripsis breakpoints and cancer rearrangements. GitHub: https://github.com/wtsi-hpag/stepStone

stepStone's coverage plot command documents `-denoise 1` as averaging data points after filtering high and low values within the window size. Pebble's default parameters (`-W 1000`, `-S 100`, `--trim-low 0.40`, `--trim-high 0.40`) match that behavior and are validated against stepStone's reference output (see `tests/test_pebble.c`).

For questions about stepStone, contact Zemin Ning (zn1@sanger.ac.uk).
