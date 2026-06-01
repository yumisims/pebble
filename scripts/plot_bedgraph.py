#!/usr/bin/env python3
"""
Plot Pebble BedGraph smoothed coverage tracks.

Writes PNG when matplotlib is available, otherwise SVG (stdlib only).

Examples:
  ./build/pebble -i genome2cov.bed -c scaffold_1 -o smoothed.bedgraph
  python3 scripts/plot_bedgraph.py smoothed.bedgraph -c scaffold_1 -o plot.png
  python3 scripts/plot_bedgraph.py smoothed.bedgraph -c scaffold_1 -o plot.svg --format svg
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable

HAS_MPL = False
try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    HAS_MPL = True
except ImportError:
    pass


def read_bedgraph(path: Path) -> dict[str, list[tuple[int, int, float]]]:
    tracks: dict[str, list[tuple[int, int, float]]] = defaultdict(list)

    with path.open("r", encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue

            parts = line.split()
            if len(parts) < 4:
                raise ValueError(f"{path}:{line_no}: expected 4 columns, got {len(parts)}")

            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            value = float(parts[3])
            tracks[chrom].append((start, end, value))

    for chrom in tracks:
        tracks[chrom].sort(key=lambda row: row[0])

    return tracks


def bedgraph_to_xy(intervals: list[tuple[int, int, float]]) -> tuple[list[int], list[float], int]:
    if not intervals:
        return [], [], 0

    xs = [intervals[0][0]]
    ys = [intervals[0][2]]
    xmax = intervals[0][1]
    for start, end, value in intervals[1:]:
        xs.append(start)
        ys.append(value)
        if end > xmax:
            xmax = end

    return xs, ys, xmax


def resolve_ymax(values: Iterable[float], ymax_arg: str) -> float:
    vals = list(values)
    if not vals:
        return 180.0
    if ymax_arg == "auto":
        peak = max(vals)
        if peak <= 0.0:
            return 180.0
        return max(10.0, peak * 1.15)
    return float(ymax_arg)


def plot_track_png(
    chrom: str,
    xs: list[int],
    ys: list[float],
    xmax: int,
    *,
    sample: str,
    ymax: float,
    width: float,
    height: float,
    output: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(width, height))
    ax.plot(xs, ys, color="black", linewidth=1.2)
    ax.set_xlabel("Chromosome coordinate")
    ax.set_ylabel("Base coverage")
    ax.set_title(f"{sample} — {chrom}")
    ax.set_ylim(0, ymax)
    ax.set_xlim(xs[0], xmax)
    ax.grid(True, alpha=0.25, linewidth=0.5)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=150)
    plt.close(fig)


def plot_track_svg(
    chrom: str,
    xs: list[int],
    ys: list[float],
    xmax: int,
    *,
    sample: str,
    ymax: float,
    width: int,
    height: int,
    output: Path,
) -> None:
    margin_left = 70
    margin_right = 20
    margin_top = 40
    margin_bottom = 50
    plot_w = width - margin_left - margin_right
    plot_h = height - margin_top - margin_bottom

    xmin = xs[0]
    xspan = max(1, xmax - xmin)

    def x_px(x: int) -> float:
        return margin_left + (x - xmin) * plot_w / xspan

    def y_px(y: float) -> float:
        clipped = min(max(y, 0.0), ymax)
        return margin_top + plot_h - (clipped * plot_h / ymax)

    points = " ".join(f"{x_px(x):.1f},{y_px(y):.1f}" for x, y in zip(xs, ys))

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as out:
        out.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        out.write(
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">\n'
        )
        out.write(f'  <rect width="100%" height="100%" fill="white"/>\n')
        out.write(
            f'  <text x="{margin_left}" y="24" font-size="14" '
            f'font-family="sans-serif">{sample} — {chrom}</text>\n'
        )
        out.write(
            f'  <line x1="{margin_left}" y1="{margin_top + plot_h}" '
            f'x2="{margin_left + plot_w}" y2="{margin_top + plot_h}" stroke="#333"/>\n'
        )
        out.write(
            f'  <line x1="{margin_left}" y1="{margin_top}" '
            f'x2="{margin_left}" y2="{margin_top + plot_h}" stroke="#333"/>\n'
        )
        out.write(
            f'  <polyline fill="none" stroke="black" stroke-width="1.5" points="{points}"/>\n'
        )
        out.write(
            f'  <text x="{margin_left + plot_w / 2:.0f}" y="{height - 12}" '
            f'text-anchor="middle" font-size="12" font-family="sans-serif">'
            f"Chromosome coordinate</text>\n"
        )
        out.write(
            f'  <text x="16" y="{margin_top + plot_h / 2:.0f}" '
            f'transform="rotate(-90 16 {margin_top + plot_h / 2:.0f})" '
            f'text-anchor="middle" font-size="12" font-family="sans-serif">'
            f"Base coverage</text>\n"
        )
        out.write(
            f'  <text x="{margin_left - 8}" y="{margin_top + plot_h + 4}" '
            f'text-anchor="end" font-size="10" font-family="sans-serif">0</text>\n'
        )
        out.write(
            f'  <text x="{margin_left - 8}" y="{margin_top + 4}" '
            f'text-anchor="end" font-size="10" font-family="sans-serif">{ymax:.0f}</text>\n'
        )
        out.write("</svg>\n")


def plot_track(
    chrom: str,
    intervals: list[tuple[int, int, float]],
    *,
    sample: str,
    ymax: float,
    width: float,
    height: float,
    output: Path,
    fmt: str,
) -> bool:
    xs, ys, xmax = bedgraph_to_xy(intervals)
    if not xs:
        print(f"pebble-plot: skipping empty track {chrom}", file=sys.stderr)
        return False

    use_svg = fmt == "svg" or (fmt == "auto" and not HAS_MPL)
    if use_svg:
        if output.suffix.lower() not in {".svg", ""}:
            output = output.with_suffix(".svg")
        plot_track_svg(
            chrom,
            xs,
            ys,
            xmax,
            sample=sample,
            ymax=ymax,
            width=int(width * 80),
            height=int(height * 80),
            output=output,
        )
    else:
        if output.suffix.lower() not in {".png", ".pdf", ".jpg", ".jpeg", ""}:
            output = output.with_suffix(".png")
        plot_track_png(
            chrom,
            xs,
            ys,
            xmax,
            sample=sample,
            ymax=ymax,
            width=width,
            height=height,
            output=output,
        )

    print(f"wrote {output} ({len(xs)} points, ymax={ymax:.1f})")
    return True


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot BedGraph smoothed coverage (Pebble output)."
    )
    parser.add_argument("bedgraph", type=Path, help="Smoothed BedGraph from pebble")
    parser.add_argument(
        "-c",
        "--chrom",
        action="append",
        dest="chroms",
        help="Contig to plot (repeatable). Default: all contigs",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output image file, or directory for multiple contigs",
    )
    parser.add_argument(
        "-s",
        "--sample",
        default="sample",
        help="Sample name for plot title (default: sample)",
    )
    parser.add_argument(
        "--ymax",
        default="auto",
        help="Y-axis max, or 'auto' from data (default: auto). Use 180 for stepStone style.",
    )
    parser.add_argument(
        "--format",
        choices=("auto", "png", "svg"),
        default="auto",
        help="Output format (default: auto → png if matplotlib installed, else svg)",
    )
    parser.add_argument(
        "--width",
        type=float,
        default=12.0,
        help="Figure width in inches for PNG (default: 12)",
    )
    parser.add_argument(
        "--height",
        type=float,
        default=4.0,
        help="Figure height in inches for PNG (default: 4)",
    )
    parser.add_argument(
        "--list-chroms",
        action="store_true",
        help="List contig names in the BedGraph and exit",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if not args.bedgraph.is_file():
        print(f"pebble-plot: file not found: {args.bedgraph}", file=sys.stderr)
        return 1

    if args.bedgraph.stat().st_size == 0:
        print(
            f"pebble-plot: input is empty: {args.bedgraph}\n"
            "  Run pebble first, e.g.:\n"
            "  ./build/pebble -i your.bed -c scaffold_1 -o smoothed.bedgraph",
            file=sys.stderr,
        )
        return 1

    try:
        tracks = read_bedgraph(args.bedgraph)
    except ValueError as exc:
        print(f"pebble-plot: {exc}", file=sys.stderr)
        return 1

    if not tracks:
        print(f"pebble-plot: no data rows in {args.bedgraph}", file=sys.stderr)
        return 1

    if args.list_chroms:
        for chrom in sorted(tracks):
            print(chrom)
        return 0

    if args.format == "png" and not HAS_MPL:
        print(
            "pebble-plot: matplotlib not installed; use --format svg or "
            "'pip install matplotlib'",
            file=sys.stderr,
        )
        return 1

    chroms = args.chroms if args.chroms else sorted(tracks.keys())
    missing = [chrom for chrom in chroms if chrom not in tracks]
    if missing:
        available = ", ".join(sorted(tracks))
        print(
            f"pebble-plot: unknown contig(s): {', '.join(missing)}\n"
            f"  Available: {available}",
            file=sys.stderr,
        )
        return 1

    written = 0
    if len(chroms) == 1:
        chrom = chroms[0]
        intervals = tracks[chrom]
        ymax = resolve_ymax((v for _s, _e, v in intervals), args.ymax)
        out = args.output
        if out.suffix == "":
            ext = ".svg" if args.format == "svg" or (args.format == "auto" and not HAS_MPL) else ".png"
            out = out / f"{args.sample}_{chrom}{ext}"
        if plot_track(
            chrom,
            intervals,
            sample=args.sample,
            ymax=ymax,
            width=args.width,
            height=args.height,
            output=out,
            fmt=args.format,
        ):
            written += 1
    else:
        out_dir = args.output
        if out_dir.suffix:
            print(
                "pebble-plot: for multiple contigs, -o must be a directory",
                file=sys.stderr,
            )
            return 1
        out_dir.mkdir(parents=True, exist_ok=True)
        for chrom in chroms:
            intervals = tracks[chrom]
            ymax = resolve_ymax((v for _s, _e, v in intervals), args.ymax)
            safe = chrom.replace("/", "_")
            ext = ".svg" if args.format == "svg" or (args.format == "auto" and not HAS_MPL) else ".png"
            if plot_track(
                chrom,
                intervals,
                sample=args.sample,
                ymax=ymax,
                width=args.width,
                height=args.height,
                output=out_dir / f"{args.sample}_{safe}{ext}",
                fmt=args.format,
            ):
                written += 1

    if written == 0:
        print("pebble-plot: no plots written", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
