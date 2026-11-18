#!/usr/bin/env python3
"""
plot_figure2_representative_introns.py

Generates Figure 2: per-site phyloP (LRT, CONACC) evolutionary-rate profiles
for four representative Cereus introns (manuscript Results 3.2 / Figure 2
legend): two "regime shift" loci (a near-neutral segment followed by a
clade-restricted accelerated block), one whole-locus acceleration example,
and one near-neutral locus shown for contrast.

Reads the per-site CSVs produced by phylopv4_linux_filtered.R
(csv_filtered/resultados_detalhados/<og>_LRT_scores.csv, one per intron).
Clade colors follow Paul Tol's "bright" colorblind-safe palette; the
outgroup is drawn first, below the ingroup points. Figure size follows the
Biological Journal of the Linnean Society limit of 168 x 225 mm.

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
    ('OG0068401_allintron', 'OG0068401\nWhole-locus acceleration —\nclades A2 + D'),
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
    fig, axes = plt.subplots(2, 2, figsize=(width_mm / 25.4, height_mm / 25.4))
    for ax, (og, title) in zip(axes.flat, PANELS):
        draw(ax, detail_dir, og, title)

    legend_handles = [
        Line2D([0], [0], marker='o', linestyle='', color=CLADE_COLORS[c],
               markersize=LEGEND_MARKERSIZE)
        for c in LEGEND_ORDER
    ]
    fig.tight_layout(rect=[0, 0, 0.83, 0.99])
    fig.legend(legend_handles, LEGEND_ORDER, title='Clade', loc='center left',
               bbox_to_anchor=(0.845, 0.5), fontsize=7, title_fontsize=8,
               frameon=False, handletextpad=0.3, labelspacing=0.4)

    renderer = fig.canvas.get_renderer()
    to_fig = fig.transFigure.inverted()
    for row, letters in zip(axes, ('AB', 'CD')):
        boxes = [to_fig.transform_bbox(ax.get_tightbbox(renderer)) for ax in row]
        top = max(b.y1 for b in boxes)
        for b, letter in zip(boxes, letters):
            fig.text(b.x0, top, letter, fontsize=11, fontweight='bold', ha='left', va='top')
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
    fig.savefig(f'{args.out_prefix}.tif', dpi=600, format='tiff',
                pil_kwargs={'compression': 'tiff_lzw'})
    fig.savefig(f'{args.out_prefix}.png', dpi=600)
    fig.savefig(f'{args.out_prefix}.pdf')
    print(f'saved {args.out_prefix}_preview.png / {args.out_prefix}.tif / .png / .pdf')


if __name__ == '__main__':
    main()
