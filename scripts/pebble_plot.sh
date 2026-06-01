#!/usr/bin/env bash
# Smooth genome2cov-style BED and plot one contig.
#
# Usage:
#   scripts/pebble_plot.sh input.bed scaffold_1 sample_name output_dir
#
# Writes:
#   output_dir/sample_name.smoothed.bedgraph
#   output_dir/sample_name_scaffold_1.svg  (or .png if matplotlib is installed)

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PEBBLE="${PEBBLE_BIN:-$ROOT/build/pebble}"

if [[ $# -lt 4 ]]; then
  echo "Usage: $0 input.bed contig sample_name output_dir" >&2
  echo "  Example: $0 genome2cov.bed scaffold_1 mysample ./plots" >&2
  exit 1
fi

INPUT="$1"
CONTIG="$2"
SAMPLE="$3"
OUTDIR="$4"

SMOOTHED="$OUTDIR/${SAMPLE}.smoothed.bedgraph"

mkdir -p "$OUTDIR"

if [[ ! -x "$PEBBLE" ]]; then
  echo "pebble_plot: binary not found at $PEBBLE (set PEBBLE_BIN or build first)" >&2
  exit 1
fi

if [[ ! -f "$INPUT" ]]; then
  echo "pebble_plot: input not found: $INPUT" >&2
  exit 1
fi

echo "==> pebble: smoothing $CONTIG from $INPUT"
"$PEBBLE" -i "$INPUT" -c "$CONTIG" -o "$SMOOTHED"

LINES=$(wc -l < "$SMOOTHED" | tr -d ' ')
if [[ "$LINES" == "0" ]]; then
  echo "pebble_plot: pebble produced empty output for $CONTIG" >&2
  echo "  Contig may be shorter than window (default 1000 bp)." >&2
  echo "  Try: $PEBBLE -i $INPUT -c $CONTIG -W 500 -S 50 -o $SMOOTHED" >&2
  exit 1
fi

echo "==> plot: $LINES bedgraph rows"
python3 "$ROOT/scripts/plot_bedgraph.py" \
  "$SMOOTHED" -c "$CONTIG" -s "$SAMPLE" \
  -o "$OUTDIR/${SAMPLE}_${CONTIG}.svg" --format svg --ymax auto

echo "Done:"
echo "  $SMOOTHED"
echo "  $OUTDIR/${SAMPLE}_${CONTIG}.svg"
