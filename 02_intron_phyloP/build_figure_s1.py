#!/usr/bin/env python3
"""
build_figure_s1.py

Composes the per-locus PNGs from plot_individual_intron_profiles.py into a
single multi-panel Figure S1, one panel per accelerated intron (17 panels,
labeled A-Q), matching Biological Journal of the Linnean Society's
Supporting Information convention of one numbered file per supplementary
figure (https://academic.oup.com/biolinnean/pages/General_Instructions:
"Supplementary figures ... should be numbered 'Figure S1' ...", files
recommended under 2 MB).

Usage
-----
    python build_figure_s1.py --input-dir individual_profiles/ \
        --out-prefix Figure_S1_accelerated_intron_profiles \
        --cols 3
"""

import argparse
import glob
import os

from PIL import Image, ImageDraw, ImageFont

TITLE = ("Figure S1. Per-site phyloP evolutionary-rate profiles for the "
          "{n} introns with sustained clade-specific acceleration.")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-dir", required=True,
                         help="Directory with one PNG per locus (plot_individual_intron_profiles.py output).")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--cols", type=int, default=3)
    parser.add_argument("--panel-width", type=int, default=850)
    parser.add_argument("--panel-height", type=int, default=420)
    args = parser.parse_args()

    files = sorted(glob.glob(os.path.join(args.input_dir, "*.png")))
    n = len(files)
    if n == 0:
        raise SystemExit(f"No PNGs found in {args.input_dir}")
    if n > 26:
        raise SystemExit("More than 26 panels: extend the A-Z labeling scheme first.")

    cols = args.cols
    rows = -(-n // cols)
    label_h = 36
    margin = 20
    gap = 16
    cell_w, cell_h = args.panel_width, args.panel_height + label_h

    grid_w = margin * 2 + cols * cell_w + (cols - 1) * gap
    grid_h = margin * 2 + rows * cell_h + (rows - 1) * gap + 60

    canvas = Image.new("RGB", (grid_w, grid_h), "white")
    draw = ImageDraw.Draw(canvas)
    try:
        font_title = ImageFont.truetype("arial.ttf", 30)
        font_label = ImageFont.truetype("arialbd.ttf", 24)
    except Exception:
        font_title = ImageFont.load_default()
        font_label = ImageFont.load_default()

    draw.text((margin, 15), TITLE.format(n=n), fill="black", font=font_title)

    letters = [chr(ord("A") + i) for i in range(n)]
    for i, f in enumerate(files):
        row, col = divmod(i, cols)
        x = margin + col * (cell_w + gap)
        y = 60 + margin + row * (cell_h + gap)
        im = Image.open(f).convert("RGB").resize((args.panel_width, args.panel_height), Image.LANCZOS)
        draw.text((x, y), letters[i], fill="black", font=font_label)
        canvas.paste(im, (x, y + label_h))

    png_path = f"{args.out_prefix}.png"
    pdf_path = f"{args.out_prefix}.pdf"
    canvas.save(png_path, optimize=True)
    canvas.save(pdf_path, "PDF", resolution=150.0)

    print(f"{n} panels, {cols}x{rows} grid")
    print(f"Wrote {png_path} ({os.path.getsize(png_path) / 1e6:.2f} MB)")
    print(f"Wrote {pdf_path} ({os.path.getsize(pdf_path) / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
