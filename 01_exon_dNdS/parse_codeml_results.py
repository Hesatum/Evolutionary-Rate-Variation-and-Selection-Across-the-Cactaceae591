#!/usr/bin/env python3
"""
parse_codeml_results.py

Parses a batch of codeml output files (as produced by run_codeml_batch_default.py
/ EasyPAML), computes likelihood-ratio tests (LRT) between nested models,
extracts Bayes Empirical Bayes (BEB) sites, and applies the Benjamini-Hochberg
false-discovery-rate (FDR) correction within each family of tests, matching
manuscript Methods 2.3:

  "Both models were validated against their respective null hypothesis using
  likelihood ratio tests with two degrees of freedom ... For branch-level
  analyses ... likelihood-ratio tests with six degrees of freedom ... the
  likelihood-ratio p-values were further adjusted for multiple comparisons
  using the Benjamini-Hochberg false discovery rate (FDR) procedure ... Sites
  with high posterior probabilities (> 0.95) were considered candidates for
  positive selection."

Site-parsing regex and lnL/omega extraction follow the same patterns used in
EasyPAML's SitesParser/CodemlBatchAnalysis (src/backend/sites_parser.py,
codeml_backend.py); the LRT/FDR aggregation step below is a clean, documented
reimplementation, since the batch-analysis GUI itself does not compute FDR.

Expected input layout (as written by run_codeml_batch_default.py):
    <results_dir>/<locus>/<model>/<locus>_<model>.txt

Usage
-----
    python parse_codeml_results.py --results-dir codeml_results/ \
        --null M1a --alt M2a --df 2 --out site_M1aM2a.tsv

    python parse_codeml_results.py --results-dir codeml_results/ \
        --null M0 --alt Branch --df 6 --out branch_M0Branch.tsv
"""

import argparse
import re
from pathlib import Path

from scipy import stats

LNL_RE = re.compile(r'lnL\([^)]*\):\s*(-?\d+\.\d+)')

# "   159 R      0.990**       8.976 +- 1.573"
BEB_SITE_RE = re.compile(r'\s*(\d+)\s+([A-Z])\s+([\d.]+)(\*{0,2})\s+([\d.]+)\s*\+-\s*([\d.]+)')


def extract_lnl(text: str) -> float | None:
    m = LNL_RE.search(text)
    return float(m.group(1)) if m else None


def extract_beb_sites(text: str, method: str = "BEB") -> list[dict]:
    """Extract 'Positively selected sites' from a codeml BEB/NEB block."""
    pattern = (
        rf"{method}\b.*?analysis"
        r".*?Positively selected sites"
        r".*?\(amino acids refer to[^)]*\)"
        r".*?Pr\(w>1\)[^\n]*\n"
        r"\s*(.*?)"
        r"(?:\n\s*\n|Time used:|$)"
    )
    match = re.search(pattern, text, re.DOTALL | re.IGNORECASE)
    if not match:
        return []
    sites = []
    for line in match.group(1).splitlines():
        m = BEB_SITE_RE.search(line)
        if not m:
            continue
        pos, aa, pr, stars, mean, se = m.groups()
        pr = float(pr)
        sites.append({
            'position': int(pos), 'amino_acid': aa, 'pr_w_gt_1': pr,
            'post_mean': float(mean), 'post_se': float(se),
            'significant_95': pr >= 0.95, 'significant_99': pr >= 0.99,
        })
    return sites


def benjamini_hochberg(pvalues: list[float]) -> list[float]:
    """Benjamini-Hochberg FDR correction (Benjamini & Hochberg 1995).

    Returns q-values in the same order as `pvalues`.
    """
    n = len(pvalues)
    order = sorted(range(n), key=lambda i: pvalues[i])
    q = [0.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        k = n - rank + 1
        val = min(prev, pvalues[i] * n / k)
        q[i] = val
        prev = val
    return q


def find_result_file(results_dir: Path, locus: str, model: str) -> Path | None:
    candidate = results_dir / locus / model / f"{locus}_{model}.txt"
    return candidate if candidate.exists() else None


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results-dir", required=True, type=Path,
                         help="Root directory of codeml runs (<locus>/<model>/<locus>_<model>.txt).")
    parser.add_argument("--null", required=True, help="Null model name (e.g. M1a, M7, M0).")
    parser.add_argument("--alt", required=True, help="Alternative model name (e.g. M2a, M8, Branch).")
    parser.add_argument("--df", required=True, type=int,
                         help="Degrees of freedom for the LRT (2 for site models, 6 for the branch model used here).")
    parser.add_argument("--out", required=True, type=Path, help="Output TSV path.")
    parser.add_argument("--fdr-alpha", type=float, default=0.05)
    args = parser.parse_args()

    loci = sorted(p.name for p in args.results_dir.iterdir() if p.is_dir())

    rows = []
    for locus in loci:
        null_f = find_result_file(args.results_dir, locus, args.null)
        alt_f = find_result_file(args.results_dir, locus, args.alt)
        if not (null_f and alt_f):
            continue

        null_text = null_f.read_text(encoding='utf-8', errors='ignore')
        alt_text = alt_f.read_text(encoding='utf-8', errors='ignore')
        lnl_null = extract_lnl(null_text)
        lnl_alt = extract_lnl(alt_text)
        if lnl_null is None or lnl_alt is None:
            continue

        lrt = 2 * (lnl_alt - lnl_null)
        p_raw = stats.chi2.sf(lrt, args.df) if lrt >= 0 else float('nan')
        beb_sites = extract_beb_sites(alt_text)
        beb_summary = "; ".join(
            f"{s['position']}{s['amino_acid']}:{s['pr_w_gt_1']:.3f}"
            f"{'**' if s['significant_99'] else '*' if s['significant_95'] else ''}"
            for s in beb_sites
        )

        rows.append({
            'locus': locus, 'lnL_null': lnl_null, 'lnL_alt': lnl_alt,
            'lrt_2dl': lrt, 'p_raw': p_raw, 'beb_sites': beb_summary or '-',
        })

    # LRT statistics that are negative for nested models are impossible and
    # signal a codeml optimisation failure (see manuscript Figure 3 legend);
    # they are kept in the table but excluded from the FDR family.
    valid = [r for r in rows if r['lrt_2dl'] >= 0]
    qvals = benjamini_hochberg([r['p_raw'] for r in valid])
    for r, q in zip(valid, qvals):
        r['fdr_q'] = q
        r['significant_fdr'] = q < args.fdr_alpha
    for r in rows:
        r.setdefault('fdr_q', float('nan'))
        r.setdefault('significant_fdr', False)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, 'w', encoding='utf-8') as fh:
        header = ['locus', 'lnL_null', 'lnL_alt', 'lrt_2dl', 'p_raw', 'fdr_q',
                   'significant_fdr', 'beb_sites']
        fh.write('\t'.join(header) + '\n')
        for r in rows:
            fh.write('\t'.join(str(r[h]) for h in header) + '\n')

    n_raw = sum(1 for r in rows if not (r['p_raw'] != r['p_raw']) and r['p_raw'] < 0.05)
    n_fdr = sum(1 for r in rows if r['significant_fdr'])
    print(f"{args.null} vs {args.alt}: n={len(rows)}, raw p<0.05: {n_raw}, "
          f"FDR<{args.fdr_alpha}: {n_fdr}")
    print(f"Written: {args.out}")


if __name__ == '__main__':
    main()
