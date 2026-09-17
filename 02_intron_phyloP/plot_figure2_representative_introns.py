"""
Figure 2 — Evolutionary-rate profiles of representative Cereus introns.

Reads the per-site phyloP LRT output produced by `phylopv4_linux_filtered.R`
(one CSV per intron under csv_filtered/resultados_detalhados/) and reproduces
the four-panel summary figure used in the manuscript: two "regime shift" loci
(a near-neutral segment followed by a clade-restricted accelerated block),
one whole-locus acceleration example, and one near-neutral locus shown for
contrast.

Clade colors follow the Okabe-Ito colorblind-safe palette, matching the one
defined in phylopv4_linux_filtered.R (`clade_colors`), so this figure is
visually consistent with every other clade-colored plot in the pipeline.

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

# Okabe-Ito colorblind-safe palette — same as clade_colors in phylopv4_linux_filtered.R
CLADE_COLORS = {
    "A1": "#E69F00",
    "A2": "#56B4E9",
    "B": "#009E73",
    "C": "#F0E442",
    "D": "#0072B2",
    "E": "#D55E00",
    "Outgroup": "#999999",
}
CLADE_ORDER = ["A1", "A2", "B", "C", "D", "E", "Outgroup"]

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
    for clade in CLADE_ORDER:
        sub = df[df["clade"] == clade].sort_values("coord")
        ax.scatter(
            sub["coord"], sub["score"],
            s=4, alpha=0.65, linewidths=0,
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

    # "center left" anchors the LEGEND'S left edge at bbox_to_anchor's x —
    # the legend body extends further right from there, so this x must
    # leave enough room for its own width or it renders clipped off-canvas.
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.tight_layout(rect=[0, 0, 0.83, 0.95])
    fig.legend(
        handles, labels, title="Clade", loc="center left",
        bbox_to_anchor=(0.845, 0.5), markerscale=2.5, fontsize=7,
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
