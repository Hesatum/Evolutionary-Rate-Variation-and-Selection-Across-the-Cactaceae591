"""
Figure 2 — Evolutionary-rate profiles of representative Cereus introns.

Reads the per-site phyloP LRT output produced by `phylopv4_linux_filtered.R`
(one CSV per intron under csv_filtered/resultados_detalhados/) and reproduces
the four-panel summary figure used in the manuscript: two "regime shift" loci
(a near-neutral segment followed by a clade-restricted accelerated block),
one whole-locus acceleration example, and one near-neutral locus shown for
contrast.

Clade colors follow Paul Tol's "bright" qualitative palette — 7 hues, checked
for deuteranopia/protanopia/tritanopia distinguishability, chosen over the
Okabe-Ito set used elsewhere in the pipeline for a more saturated, print-ready
look. The Outgroup (grey, not the object of the study) is drawn first so it
sits in the bottom z-order layer, underneath every ingroup clade.

Output is written both as a quick-look PNG and as a publication-ready TIFF
sized and resolved per the Biological Journal of the Linnean Society figure
guidelines (max size 168 x 225 mm; >=600 dpi for combination
line/graph figures; uncompressed TIFF).

Usage:
    python plot_figure2_representative_introns.py
"""

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DETAIL_DIR = "../csv_filtered/resultados_detalhados"  # adjust if run from elsewhere
OUT_PNG = "Figure2_representative_introns_preview.png"
OUT_TIFF = "Figure2.tif"

# Paul Tol's "bright" qualitative palette (colorblind-safe, 7 categories).
CLADE_COLORS = {
    "A1": "#4477AA",  # blue
    "A2": "#66CCEE",  # cyan
    "B":  "#228833",  # green
    "C":  "#CCBB44",  # yellow
    "D":  "#EE6677",  # red
    "E":  "#AA3377",  # purple
    "Outgroup": "#BBBBBB",  # grey — not the object of the study, drawn first (bottom layer)
}
LEGEND_ORDER = ["A1", "A2", "B", "C", "D", "E", "Outgroup"]
DRAW_ORDER = ["Outgroup", "A1", "A2", "B", "C", "D", "E"]

POINT_SIZE = 1.8       # scatter marker size in the plot (small — 7 clades overlap a lot)
LEGEND_MARKERSIZE = 7  # legend marker size in points, independent of POINT_SIZE

SIG_THRESHOLD = 1.3  # phyloP score threshold, P < 0.05

PANELS = [
    ("OG0081328_allintron", "OG0081328\nRegime shift — clade C"),
    ("OG0078211_allintron", "OG0078211\nRegime shift — clade D"),
    ("OG0068401_allintron", "OG0068401\nWhole-locus acceleration — clades A2 + D"),
    ("OG0094833_allintron", "OG0094833\nNear-neutral profile"),
]


def load_scores(og: str) -> pd.DataFrame:
    return pd.read_csv(f"{DETAIL_DIR}/{og}_LRT_scores.csv")


def draw(ax, og: str, title: str) -> None:
    df = load_scores(og)
    # Draw in DRAW_ORDER (Outgroup first -> bottom z-layer) so the ingroup
    # clades, which are the actual focus of the study, are never hidden
    # underneath Outgroup points.
    for zorder, clade in enumerate(DRAW_ORDER):
        sub = df[df["clade"] == clade].sort_values("coord")
        ax.scatter(
            sub["coord"], sub["score"],
            s=POINT_SIZE, alpha=0.65, linewidths=0, zorder=zorder,
            color=CLADE_COLORS[clade], label=clade,
        )
    ax.axhline(0, color="black", lw=0.8)
    ax.axhline(SIG_THRESHOLD, color="red", ls="--", lw=0.8)
    ax.axhline(-SIG_THRESHOLD, color="red", ls="--", lw=0.8)
    ax.set_title(title, fontsize=9, fontweight="bold")
    ax.set_xlabel("Alignment position (bp)", fontsize=8)
    ax.set_ylabel("PhyloP score (LRT)", fontsize=8)
    ax.tick_params(labelsize=7)


def build_figure(width_mm: float = 165, height_mm: float = 135):
    # BJLS max published figure size is 168 x 225 mm — the legend must be
    # reserved INSIDE that canvas (not added via bbox_inches="tight" on
    # save, which would silently grow the file past the limit).
    fig, axes = plt.subplots(
        2, 2, figsize=(width_mm / 25.4, height_mm / 25.4)
    )
    for ax, (og, title) in zip(axes.flat, PANELS):
        draw(ax, og, title)

    # Legend uses its own fixed-size proxy markers (Line2D), independent of
    # POINT_SIZE, in LEGEND_ORDER regardless of DRAW_ORDER — so shrinking
    # the on-plot points to fight overplotting never shrinks the legend.
    from matplotlib.lines import Line2D
    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="", color=CLADE_COLORS[c],
               markersize=LEGEND_MARKERSIZE)
        for c in LEGEND_ORDER
    ]

    # "center left" anchors the LEGEND'S left edge at bbox_to_anchor's x —
    # the legend body extends further right from there, so this x must
    # leave enough room for its own width or it renders clipped off-canvas.
    fig.tight_layout(rect=[0, 0, 0.83, 0.95])
    fig.legend(
        legend_handles, LEGEND_ORDER, title="Clade", loc="center left",
        bbox_to_anchor=(0.845, 0.5), fontsize=7,
        title_fontsize=8, frameon=False, handletextpad=0.3, labelspacing=0.4,
    )
    fig.suptitle(
        "Evolutionary rate profiles of representative Cereus introns",
        fontsize=10, fontweight="bold",
    )
    return fig


if __name__ == "__main__":
    fig = build_figure()

    # Quick-look PNG for checking on screen / in the repo.
    fig.savefig(OUT_PNG, dpi=150)
    print(f"Preview salvo em {OUT_PNG}")

    # Submission-ready TIFF: BJLS requires >=300 dpi (photographs) / higher
    # for combination figures with line art and labels, max published size
    # 168 x 225 mm, TIFF format. Saved at the figure's exact declared size
    # (no bbox_inches="tight", which would silently grow past the limit)
    # with lossless LZW compression to keep the file a manageable size.
    fig.savefig(
        OUT_TIFF, dpi=600, format="tiff",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    print(f"TIFF para submissao salvo em {OUT_TIFF}")
