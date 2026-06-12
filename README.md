# Pebble

**Authors:** Yumi Sims, Danil Zilov, Jo Wood, Zemin Ning

Pebble is a command-line tool and C library for smoothing genomic sequencing coverage tracks. It reads per-base coverage from BED or BedGraph files, applies a sliding trimmed mean, and writes the result as BedGraph or BigWig.

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

Contigs shorter than the window size are skipped. The fourth column of every BedGraph row is a **rounded integer** in both smoothed and normalised output.

### Coverage normalisation

After smoothing all contigs, Pebble computes a **length-weighted genome average** of the smoothed track, then writes a second **normalised** BedGraph:

1. **Genome average** — weighted mean of smoothed coverage across all intervals (weighted by interval length).
2. **Normalise factor** — `genome_average / 2`.
3. **Normalised coverage** — `smoothed_coverage / normalise_factor`.

The whole-genome mean of the normalised track is therefore **2.0**. This matches the diploid baseline used in stepStone-style CNV plots: autosomes at typical depth map to **~2**, not to 1.

**Why divide by 2?** In a diploid genome, each autosome is present in two copies. Scaling so the genome-wide average equals **2** turns coverage into copy-number-like units:

| Relative depth | Normalised value | Interpretation |
|----------------|------------------|----------------|
| Typical autosomal | ~2 | Diploid |
| Half autosomal depth | ~1 | One copy |
| Absent / unmapped | ~0 | No copy |

This makes broad deviations from 2 easy to see on a plot without converting raw read depth mentally.

**Sex chromosome identification.** Because the normalise factor is driven mostly by autosomes (diploid, ~2 after scaling), sex chromosomes stand out when their true copy number differs:

| Karyotype | Autosomes | X | Y |
|-----------|-----------|---|---|
| XY (male) | ~2 | ~1 | ~1 |
| XX (female) | ~2 | ~2 | ~0 |

On a normalised track, an X contig sitting near **1** while autosomes sit near **2** suggests XY; X near **2** with Y near **0** suggests XX. This is a visual screen, not a formal karyotype call — PAR regions, aneuploidy, or misassembly can affect the pattern — but it is the same rationale stepStone uses for denoised coverage plots.

Pebble logs the genome average and normalise factor to stderr, for example:

```
pebble: genome average=54.24, normalise factor=27.1180 (avg/2), target mean=2.0
```

**Output files.** Pebble always writes **both** BedGraph tracks:

| File | Contents |
|------|----------|
| `<basename>.smoothed.bedgraph` | Raw smoothed coverage |
| `<basename>.normalised.bedgraph` | Normalised coverage |

The basename comes from `-o` (extension stripped) or from the input filename when `-o` is omitted.

Examples:

- `-i cov.bed` → `cov.smoothed.bedgraph` + `cov.normalised.bedgraph`
- `-o build/out.bedgraph` → `build/out.smoothed.bedgraph` + `build/out.normalised.bedgraph`
- `-o build/sample.smoothed.bedgraph` → `build/sample.smoothed.bedgraph` + `build/sample.normalised.bedgraph`

BigWig (with `--sizes`) is built from the **normalised** BedGraph. Demo mode (`--demo`) prints both smoothed and normalised columns to stdout.

## Build

Requires a C11 compiler.

Clone with submodules (needed for BigWig builds with `make PEBBLE_BIGWIG=1`):

```bash
git clone --recurse-submodules <repo-url>
# or after a plain clone:
git submodule update --init third_party/libBigWig
```

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

BigWig output is **optional** and off by default. Enable it when building if you need BigWig with `--sizes examples/chrom.sizes`:

```bash
cmake -S . -B build -DPEBBLE_BIGWIG=ON
cmake --build build
# or: make PEBBLE_BIGWIG=1
```

On Linux with static linking (the default), BigWig builds need the **static** zlib archive (`libz.a`, provided by `zlib1g-dev` on Debian/Ubuntu). If `libz.a` is missing, CMake will report an error; the Makefile falls back to a dynamic zlib link without `-static`.

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

Example data in `examples/`:

| File | Description |
|------|-------------|
| `example_cov.bed` | genome2cov-style 4-column BED (`chrom start end coverage`) for `scaffold_1` (~10k intervals) |
| `chrom.sizes` | Chromosome/scaffold lengths (88 sequences) for BigWig output |

```bash
# Built-in demo with mock coverage data
./build/pebble --demo

# Smooth example coverage (genome2cov-style BED)
./build/pebble -i examples/example_cov.bed -o build/example_cov.smoothed.bedgraph
# -> build/example_cov.smoothed.bedgraph + build/example_cov.normalised.bedgraph

# Single contig (example_cov.bed contains scaffold_1 only)
./build/pebble -i examples/example_cov.bed -c scaffold_1 -o build/scaffold_1.smoothed.bedgraph

# Normalised BedGraph + BigWig (requires -DPEBBLE_BIGWIG=ON build)
./build/pebble -i examples/example_cov.bed \
  -o build/example_cov.smoothed.bedgraph \
  --sizes examples/chrom.sizes
# -> build/example_cov.smoothed.bedgraph + build/example_cov.normalised.bedgraph + build/example_cov.bw

# BedGraph to BigWig without smoothing (pass-through; input must be .bedgraph)
./build/pebble -i raw.bedgraph \
  -o build/raw.bw \
  --sizes examples/chrom.sizes \
  --no-smooth
```

Other input types:

```bash
# BedGraph input
./build/pebble -i track.bedgraph -o smoothed.bedgraph

# Standard 6-column BED (score in column 5)
./build/pebble -i regions.bed -f bed -o smoothed.bedgraph

# Custom window/step (e.g. short contigs)
./build/pebble -i examples/genome2cov_style.bed -W 50 -S 10 -o smoothed.bedgraph
```

### Coverage plots (stepStone style)

Pebble smoothing matches [stepStone](https://github.com/wtsi-hpag/stepStone) `plot -denoise 1` (1000 bp window, mean of middle 20% after trimming). `--stepstone` uses stepStone line style (black, linewidth 3, ymin 1, ~9000 display points); Y max scales from data (stepStone’s default `-hight 180` is for shallow cancer WGS — use `--ymax 180` only when depth stays below that):

```bash
# Smooth + plot one contig (plots the normalised track)
scripts/pebble_plot.sh examples/example_cov.bed scaffold_1 mysample ./plots
# -> ./plots/mysample.smoothed.bedgraph, ./plots/mysample.normalised.bedgraph, ./plots/mysample_scaffold_1.png

# Plot an existing normalised BedGraph
python3 scripts/plot_bedgraph.py build/example_cov.normalised.bedgraph \
  -c scaffold_1 -s mysample -o plot.png --stepstone
```

For deep coverage (e.g. peaks ~800), rely on auto ymax or set explicitly: `--ymax 900`.

### Options

| Flag | Description | Default |
|------|-------------|---------|
| `-i input` | Input BED or BedGraph file | — |
| `-o output` | Output basename (writes `<basename>.smoothed.bedgraph` and `<basename>.normalised.bedgraph`; default: input basename) | input basename |
| `--sizes file` | Chromosome sizes (`name<TAB>length`); also writes a `.bw` alongside BedGraph | — |
| `--no-smooth` | Convert coverage BED/BedGraph to BigWig without smoothing (requires `--sizes`, `-o`) | off |
| `--format-out bedgraph\|bigwig` | With `--sizes`, treat `-o` as BigWig basename when set to `bigwig` | BedGraph always; BigWig when `--sizes` given |
| `-c chrom` | Process only the named chromosome or contig | all |
| `-f bed\|bedgraph` | Input format | inferred from extension |
| `-W window` | Sliding window size | 1000 |
| `-S step` | Step size | 100 |
| `--trim-low frac` | Fraction trimmed from the low end | 0.40 |
| `--trim-high frac` | Fraction trimmed from the high end | 0.40 |
| `--demo` | Run the built-in mock scaffold demo | — |

If no input file is given, Pebble runs in demo mode.

BigWig output uses [libBigWig](https://github.com/dpryan79/libBigWig). Enable it at build time with `-DPEBBLE_BIGWIG=ON` (CMake) or `make PEBBLE_BIGWIG=1`. **BedGraph is always written** (smoothed and normalised when `-o` is a file). When `--sizes` is given, Pebble also writes a `.bw` file from the **normalised** BedGraph (same intervals as UCSC `bedGraphToBigWig`). Every contig in the BigWig must appear in the chromosome sizes file.

You can also convert an existing BedGraph with UCSC tools if you prefer:

```bash
bedGraphToBigWig smoothed.bedgraph chrom.sizes smoothed.bw
```

## Input and output

**Input formats**

| Format | Columns | Notes |
|--------|---------|-------|
| BedGraph | `chrom start end value` | Inferred from `.bedgraph` extension |
| genome2cov BED | `chrom start end coverage` | 4-column; inferred from `.bed` extension |
| BED | `chrom start end name score …` | Score taken from column 5; use `-f bed` |

Pebble builds a dense per-base coverage array for each contig from the input intervals (anchored at coordinate 0), fills gaps between intervals with the previous coverage value, and extends from base 0 up to the first interval with that interval’s value. Then it applies the sliding trimmed mean.

**Output**

| Format | When | Notes |
|--------|------|-------|
| BedGraph (smoothed) | always | `<basename>.smoothed.bedgraph` |
| BedGraph (normalised) | always | `<basename>.normalised.bedgraph`; genome-wide mean ≈ 2 |
| BigWig | with `--sizes` (from normalised BedGraph), or `--no-smooth` (direct BedGraph conversion) | Requires `-o` and `-DPEBBLE_BIGWIG=ON` build |

BedGraph intervals use `[window_start, window_start + step)` and coverage values are rounded integers. BigWig carries the same intervals and values in binary form; the chromosome sizes file supplies the sequence lengths for the BigWig header (UCSC `chrom.sizes` format: one `name<TAB>length` per line, `#` comments allowed).

**Chromosome sizes file** — tab- or space-separated, one sequence per line (see `examples/chrom.sizes`):

```
scaffold_1	418431624
scaffold_10	135097871
scaffold_176	7000
```

## Library API

The core logic lives in `pebble_core`. The public headers are:

- `include/pebble/pebble.h` — `pebble_process()` and configuration
- `include/pebble/pebble_io.h` — BED/BedGraph read and write helpers
- `include/pebble/pebble_bigwig.h` — BigWig output (when built with `PEBBLE_BIGWIG`)

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

Bundled data (`examples/`):

```
example_cov.bed   # genome2cov BED for scaffold_1
chrom.sizes       # scaffold lengths for BigWig
```

Quick runs:

```bash
make example   # mock scaffold demo (small; see examples/mock_scaffold.bedgraph)

# Real example coverage
./build/pebble -i examples/example_cov.bed -o build/example_cov.smoothed.bedgraph

# Normalised BedGraph + BigWig (after cmake -DPEBBLE_BIGWIG=ON or make PEBBLE_BIGWIG=1)
./build/pebble -i examples/example_cov.bed \
  -o build/example_cov.smoothed.bedgraph \
  --sizes examples/chrom.sizes
# -> build/example_cov.smoothed.bedgraph, build/example_cov.normalised.bedgraph, build/example_cov.smoothed.bw
```

Unit tests also use smaller fixtures (`examples/mock_scaffold.bedgraph`, `examples/genome2cov_style.bed`).

## Reference

Pebble implements the sliding-window trimmed-mean smoothing used in stepStone's `-denoise 1` mode. If you use Pebble in your work, please cite stepStone:

> Zemin Ning. stepStone: a pipeline for identification of chromothripsis breakpoints and cancer rearrangements. GitHub: https://github.com/wtsi-hpag/stepStone

stepStone's coverage plot command documents `-denoise 1` as averaging data points after filtering high and low values within the window size. Pebble's default parameters (`-W 1000`, `-S 100`, `--trim-low 0.40`, `--trim-high 0.40`) match that trimming behaviour (middle 20% mean of a 1000 bp window). See `tests/test_pebble.c` for validation.

For questions about stepStone, contact Zemin Ning (zn1@sanger.ac.uk).
