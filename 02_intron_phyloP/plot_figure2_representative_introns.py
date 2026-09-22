#!/usr/bin/env python3
"""
plot_figure2_representative_introns.py

Generates Figure 2: per-site phyloP (LRT, CONACC) evolutionary-rate profiles
for four representative Cereus introns (manuscript Results 3.2 / Figure 2
legend): two "regime shift" loci (a near-neutral segment followed by a
clade-restricted accelerated block), one whole-locus acceleration example,
and one near-neutral locus shown for contrast.

Reads the per-site CSVs produced by phylopv4_linux_filtered.R
(csv_filtered/resultados_detalhados/<og>_LRT_scores.csv, one per intron,
columns include coord/score/clade). Clade colors are Paul Tol's "bright"
qualitative palette (colorblind-safe for deuteranopia/protanopia/
tritanopia). Outgroup is drawn first (bottom z-order) since it is not the
object of the study and would otherwise sit on top of, and hide, ingroup
points; its legend entry keeps the natural clade order regardless. Plot
points are sized separately from the legend markers (fixed-size Line2D
proxies) so shrinking points to fight overplotting doesn't shrink the key.

The TIFF output is sized/resolved per the Biological Journal of the Linnean
Society figure guidelines (max published size 168 x 225 mm; TIFF, LZW
lossless compression to keep the file small; the journal disallows lossy
formats like JPEG, not lossless compression).

Usage
-----
    python plot_figure2_representative_introns.py \
        --detail-dir ../csv_filtered/resultados_detalhados \
        --out-prefix Figure2
"""

import argparse

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Paul Tol's "bright" qualitative palette (colorblind-safe, 7 categories).
CLADE_COLORS = {
    'A1': '#4477AA',  # blue
    'A2': '#66CCEE',  # cyan
    'B':  '#228833',  # green
    'C':  '#CCBB44',  # yellow
    'D':  '#EE6677',  # red
    'E':  '#AA3377',  # purple
    'Outgroup': '#BBBBBB',  # grey, not the object of the study
}
LEGEND_ORDER = ['A1', 'A2', 'B', 'C', 'D', 'E', 'Outgroup']
DRAW_ORDER = ['Outgroup', 'A1', 'A2', 'B', 'C', 'D', 'E']

POINT_SIZE = 1.8       # scatter marker size on the plot (7 clades overplot heavily)
LEGEND_MARKERSIZE = 7  # legend marker size in points, independent of POINT_SIZE
SIG_THRESHOLD = 1.3    # phyloP score threshold, P < 0.05

PANELS = [
    ('OG0081328_allintron', 'OG0081328\nRegime shift — clade C'),
    ('OG0078211_allintron', 'OG0078211\nRegime shift — clade D'),
    ('OG0068401_allintron', 'OG0068401\nWhole-locus acceleration — clades A2 + D'),
    ('OG0094833_allintron', 'OG0094833\nNear-neutral profile'),
]


def load_scores(detail_dir: str, og: str) -> pd.DataFrame:
    return pd.read_csv(f'{detail_dir}/{og}_LRT_scores.csv')


def draw(ax, detail_dir: str, og: str, title: str) -> None:
    df = load_scores(detail_dir, og)
    for zorder, clade in enumerate(DRAW_ORDER):
        sub = df[df['clade'] == clade].sort_values('coord')
        ax.scatter(sub['coord'], sub['score'], s=POINT_SIZE, alpha=0.65,
                   linewidths=0, zorder=zorder, color=CLADE_COLORS[clade])
    ax.axhline(0, color='black', lw=0.8)
    ax.axhline(SIG_THRESHOLD, color='red', ls='--', lw=0.8)
    ax.axhline(-SIG_THRESHOLD, color='red', ls='--', lw=0.8)
    ax.set_title(title, fontsize=9, fontweight='bold')
    ax.set_xlabel('Alignment position (bp)', fontsize=8)
    ax.set_ylabel('PhyloP score (LRT)', fontsize=8)
    ax.tick_params(labelsize=7)


def build_figure(detail_dir: str, width_mm: float, height_mm: float):
    # BJLS max published figure size is 168 x 225 mm; the legend is reserved
    # inside that canvas via tight_layout's rect, not added afterwards with
    # bbox_inches='tight' on save, which would silently grow the file past it.
    fig, axes = plt.subplots(2, 2, figsize=(width_mm / 25.4, height_mm / 25.4))
    for ax, (og, title) in zip(axes.flat, PANELS):
        draw(ax, detail_dir, og, title)

    legend_handles = [
        Line2D([0], [0], marker='o', linestyle='', color=CLADE_COLORS[c],
               markersize=LEGEND_MARKERSIZE)
        for c in LEGEND_ORDER
    ]
    fig.tight_layout(rect=[0, 0, 0.83, 0.95])
    # 'center left' anchors the legend's left edge at bbox_to_anchor's x; the
    # body extends rightward from there, so this x must leave room for its
    # own width or it renders clipped past the canvas edge.
    fig.legend(legend_handles, LEGEND_ORDER, title='Clade', loc='center left',
               bbox_to_anchor=(0.845, 0.5), fontsize=7, title_fontsize=8,
               frameon=False, handletextpad=0.3, labelspacing=0.4)
    fig.suptitle('Evolutionary rate profiles of representative Cereus introns',
                 fontsize=10, fontweight='bold')
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--detail-dir', default='../csv_filtered/resultados_detalhados',
                         help='Directory with <og>_LRT_scores.csv per intron.')
    parser.add_argument('--out-prefix', default='Figure2')
    parser.add_argument('--width-mm', type=float, default=165)
    parser.add_argument('--height-mm', type=float, default=135)
    args = parser.parse_args()

    fig = build_figure(args.detail_dir, args.width_mm, args.height_mm)

    fig.savefig(f'{args.out_prefix}_preview.png', dpi=150)
    # Submission TIFF: no bbox_inches='tight' (would grow past the 168 x 225 mm
    # limit); LZW is lossless, unlike the JPEG/PPT/DOC formats BJLS disallows.
    fig.savefig(f'{args.out_prefix}.tif', dpi=600, format='tiff',
                pil_kwargs={'compression': 'tiff_lzw'})
    print(f'saved {args.out_prefix}_preview.png / {args.out_prefix}.tif')


if __name__ == '__main__':
    main()
