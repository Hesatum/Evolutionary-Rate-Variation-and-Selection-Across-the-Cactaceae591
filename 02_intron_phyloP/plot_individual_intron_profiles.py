#!/usr/bin/env python3
"""
plot_individual_intron_profiles.py

Generates one per-site phyloP profile plot per intron (all 7 clades in a
single panel, method=LRT, mode=CONACC), for every orthogroup listed in an
exclude-list CSV — by default, the 17 introns with sustained clade-specific
acceleration (03_species_tree/exclude_lists/introns_accelerated_17.csv).

Same drawing conventions as plot_figure2_representative_introns.py (Paul
Tol's "bright" colorblind-safe palette, Outgroup drawn first/bottom z-order
since it isn't the object of the study, |score| > 1.3 significance lines),
factored out here so a single-panel version doesn't drift from the Figure 2
four-panel version.

Usage
-----
    python plot_individual_intron_profiles.py \
        --detail-dir ../csv_filtered/resultados_detalhados \
        --exclude-list ../03_species_tree/exclude_lists/introns_accelerated_17.csv \
        --out-dir individual_profiles
"""

import argparse
import csv
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Paul Tol's "bright" qualitative palette (colorblind-safe, 7 categories) —
# same as plot_figure2_representative_introns.py.
CLADE_COLORS = {
    'A1': '#4477AA',  # blue
    'A2': '#66CCEE',  # cyan
    'B':  '#228833',  # green
    'C':  '#CCBB44',  # yellow
    'D':  '#EE6677',  # red
    'E':  '#AA3377',  # purple
    'Outgroup': '#BBBBBB',  # grey — not the object of the study
}
LEGEND_ORDER = ['A1', 'A2', 'B', 'C', 'D', 'E', 'Outgroup']
DRAW_ORDER = ['Outgroup', 'A1', 'A2', 'B', 'C', 'D', 'E']

POINT_SIZE = 2.5
LEGEND_MARKERSIZE = 7
SIG_THRESHOLD = 1.3  # phyloP score threshold, P < 0.05


def read_og_list(path: Path) -> list[str]:
    with open(path, encoding='utf-8') as f:
        rows = list(csv.reader(f))
    header = [c.strip().lower() for c in rows[0]] if rows else []
    data_rows = rows[1:] if 'og' in header else rows
    og_col = header.index('og') if 'og' in header else 0
    return [r[og_col].strip() for r in data_rows if r and r[og_col].strip()]


def load_scores(detail_dir: str, og: str) -> pd.DataFrame:
    return pd.read_csv(f'{detail_dir}/{og}_allintron_LRT_scores.csv')


def plot_one(detail_dir: str, og: str, out_dir: Path) -> None:
    df = load_scores(detail_dir, og)
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    for zorder, clade in enumerate(DRAW_ORDER):
        sub = df[df['clade'] == clade].sort_values('coord')
        ax.scatter(sub['coord'], sub['score'], s=POINT_SIZE, alpha=0.65,
                   linewidths=0, zorder=zorder, color=CLADE_COLORS[clade])
    ax.axhline(0, color='black', lw=0.8)
    ax.axhline(SIG_THRESHOLD, color='red', ls='--', lw=0.8)
    ax.axhline(-SIG_THRESHOLD, color='red', ls='--', lw=0.8)

    n = len(df) // 7  # sites per clade (same alignment length for all 7)
    accel_any = df.groupby('coord')['score'].min().lt(-SIG_THRESHOLD).mean() * 100
    ax.set_title(f'{og}\nN={n} sites | {accel_any:.1f}% accelerated in >=1 clade',
                 fontsize=10, fontweight='bold')
    ax.set_xlabel('Alignment position (bp)')
    ax.set_ylabel('PhyloP score (LRT)')

    legend_handles = [
        Line2D([0], [0], marker='o', linestyle='', color=CLADE_COLORS[c],
               markersize=LEGEND_MARKERSIZE)
        for c in LEGEND_ORDER
    ]
    ax.legend(legend_handles, LEGEND_ORDER, title='Clade', fontsize=7,
              title_fontsize=8, frameon=False, loc='upper right', ncol=2)

    plt.tight_layout()
    fig.savefig(out_dir / f'{og}_profile_LRT.png', dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--detail-dir', default='../csv_filtered/resultados_detalhados',
                         help='Directory with <og>_allintron_LRT_scores.csv per intron.')
    parser.add_argument('--exclude-list', type=Path,
                         default=Path('../03_species_tree/exclude_lists/introns_accelerated_17.csv'),
                         help='CSV with one orthogroup ID per row (og column); default is the '
                              '17 accelerated introns.')
    parser.add_argument('--out-dir', type=Path, default=Path('individual_profiles'))
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    ogs = read_og_list(args.exclude_list)
    print(f'{len(ogs)} orthogroups in {args.exclude_list}')

    for og in ogs:
        plot_one(args.detail_dir, og, args.out_dir)
        print(f'  {og} -> {args.out_dir}/{og}_profile_LRT.png')

    print(f'Done: {len(ogs)} profiles in {args.out_dir}')


if __name__ == '__main__':
    main()
