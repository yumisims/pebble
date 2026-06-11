#!/usr/bin/env bash
# Smooth genome2cov-style BED (or plot an existing BedGraph) — stepStone-style plot.
#
# Usage:
#   scripts/pebble_plot.sh input.bed scaffold_1 sample_name output_dir
#   scripts/pebble_plot.sh pebble_1k.bedgraph scaffold_1 sample_name output_dir
#
# Smoothing matches stepStone -denoise 1 (W=1000, step=100, trimmed mean).
# Plotting matches stepStone defaults (-hight 180, black line, ~9000 display points).
#
# Writes:
#   output_dir/sample_name.smoothed.bedgraph
#   output_dir/sample_name.normalised.bedgraph
#   output_dir/sample_name_scaffold_1.png  (or .svg without matplotlib)

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
NORMALISED="$OUTDIR/${SAMPLE}.normalised.bedgraph"

mkdir -p "$OUTDIR"

if [[ ! -f "$INPUT" ]]; then
  echo "pebble_plot: input not found: $INPUT" >&2
  exit 1
fi

if [[ "$INPUT" == *.bedgraph ]]; then
  echo "==> input is BedGraph; plotting $CONTIG directly (skip pebble smoothing)"
  SMOOTHED="$INPUT"
else
  if [[ ! -x "$PEBBLE" ]]; then
    echo "pebble_plot: binary not found at $PEBBLE (set PEBBLE_BIN or build first)" >&2
    exit 1
  fi

  echo "==> pebble: smoothing $CONTIG from $INPUT (stepStone denoise 1: -W 1000 -S 100)"
  "$PEBBLE" -i "$INPUT" -c "$CONTIG" -W 1000 -S 100 -o "$SMOOTHED"

  LINES=$(wc -l < "$SMOOTHED" | tr -d ' ')
  if [[ "$LINES" == "0" ]]; then
    echo "pebble_plot: pebble produced empty output for $CONTIG" >&2
    echo "  Contig may be shorter than window (default 1000 bp)." >&2
    echo "  Try: $PEBBLE -i $INPUT -c $CONTIG -W 500 -S 50 -o $SMOOTHED" >&2
    exit 1
  fi
fi

LINES=$(wc -l < "$SMOOTHED" | tr -d ' ')
if [[ "$LINES" == "0" ]]; then
  echo "pebble_plot: BedGraph is empty: $SMOOTHED" >&2
  exit 1
fi

PLOT_OUT="$OUTDIR/${SAMPLE}_${CONTIG}.png"
echo "==> plot: normalised track (stepStone style; ymax auto from data)"
python3 "$ROOT/scripts/plot_bedgraph.py" \
  "$NORMALISED" -c "$CONTIG" -s "$SAMPLE" \
  -o "$PLOT_OUT" --stepstone

if [[ ! -f "$PLOT_OUT" ]]; then
  PLOT_OUT="$OUTDIR/${SAMPLE}_${CONTIG}.svg"
fi

echo "Done:"
echo "  $SMOOTHED"
echo "  $NORMALISED"
echo "  $PLOT_OUT"
