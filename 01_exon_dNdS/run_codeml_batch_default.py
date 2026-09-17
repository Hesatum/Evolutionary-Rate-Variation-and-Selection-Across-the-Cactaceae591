#!/usr/bin/env python3
"""
run_codeml_batch_default.py

Batch execution of PAML's codeml across a directory of per-locus codon
alignments, using the default control-file (.ctl) parameters of EasyPAML
(https://github.com/<user>/EasyPAML) for the site models (M0, M1a, M2a, M7,
M8) and the branch model (Branch) used in this study.

This script is a minimal, standalone distillation of the batch/default-ctl
logic in EasyPAML's `codeml_backend.py` (CodemlBatchAnalysis.MODEL_CONFIGS
and generate_ctl_content). It is published here only to document how the
codeml runs reported in the manuscript were parameterized and batched; the
full EasyPAML application (GUI, LRT/BEB parsing, warm-start heuristics,
report generation, etc.) is maintained and distributed separately.

Requirements
------------
- PAML's `codeml` executable available on PATH (or pass --codeml-bin).
- One fixed tree file (Newick), shared across all loci, matching the
  IQ-TREE topology used as input for all codeml analyses (see manuscript
  Methods 2.3).
- One codon alignment per locus in the input directory (PHYLIP or FASTA,
  in-frame, gap-free per manuscript Methods 2.2), named "<locus>.<ext>".

Usage
-----
    python run_codeml_batch_default.py \\
        --input-dir alignments/ \\
        --tree species_tree.nwk \\
        --output-dir codeml_results/ \\
        --models M0 M1a M2a M7 M8 Branch \\
        --threads 4

Each locus/model combination is run in its own working directory
(<output_dir>/<locus>/<model>/), with a freshly generated .ctl file, so
that runs can be executed in parallel without clobbering codeml's
temporary files (rst, rub, 2NG.dN, etc.).
"""

import argparse
import concurrent.futures
import shutil
import subprocess
from pathlib import Path

# ---------------------------------------------------------------------------
# Default model parameters, exactly as used by EasyPAML for this study.
# Only the models actually reported in the manuscript are included:
#   - Site models (2.3): M0 (baseline), M1a vs M2a, M7 vs M8
#   - Branch model (2.3): M0 (null) vs Branch (free omega among clades)
# ---------------------------------------------------------------------------
MODEL_CONFIGS = {
    "M0": {"model": 0, "NSsites": 0, "fix_omega": 0, "omega": 0.5, "CodonFreq": 7},
    "M1a": {"model": 0, "NSsites": 1, "fix_omega": 0, "omega": 0.5, "CodonFreq": 7},
    "M2a": {"model": 0, "NSsites": 2, "fix_omega": 0, "omega": 0.5, "CodonFreq": 7},
    "M7": {"model": 0, "NSsites": 7, "fix_omega": 0, "omega": 0.5, "CodonFreq": 7},
    "M8": {"model": 0, "NSsites": 8, "fix_omega": 0, "omega": 0.5, "CodonFreq": 7},
    "Branch": {"model": 2, "NSsites": 0, "fix_omega": 0, "omega": 0.5, "CodonFreq": 7},
}

ALIGNMENT_EXTENSIONS = (".fas", ".fasta", ".phy", ".phylip")


def generate_ctl_content(seqfile: str, treefile: str, outfile: str,
                          model_config: dict, cleandata: int = 1) -> str:
    """Reproduce EasyPAML's default .ctl template for a given model."""
    return f"""      seqfile = {seqfile}
     treefile = {treefile}
      outfile = {outfile}

        noisy = 1              * 0-9: output detail
      verbose = 1              * More or less detailed report in outfile
      seqtype = 1              * Codon data
        ndata = 1              * Number of data sets or loci
        icode = 0              * Universal genetic code
    cleandata = {cleandata}              * Remove sites with ambiguity data?

        model = {model_config['model']}         * Models for omega varying across lineages
      NSsites = {model_config['NSsites']}          * Models for omega varying across sites
    CodonFreq = {model_config['CodonFreq']}        * Codon frequencies (F3x4)
      estFreq = 0              * Use observed frequencies
        clock = 0              * No molecular clock
    fix_omega = {model_config['fix_omega']}         * Estimate omega
        omega = {model_config['omega']}        * Initial omega value
"""


def find_alignments(input_dir: Path):
    files = []
    for ext in ALIGNMENT_EXTENSIONS:
        files.extend(sorted(input_dir.glob(f"*{ext}")))
    return files


def run_one(locus_file: Path, tree_file: Path, model_name: str,
            output_dir: Path, codeml_bin: str) -> tuple[str, str, bool, str]:
    """Generate the .ctl and run codeml for a single locus/model pair."""
    locus = locus_file.stem
    run_dir = output_dir / locus / model_name
    run_dir.mkdir(parents=True, exist_ok=True)

    seq_dst = run_dir / locus_file.name
    tree_dst = run_dir / tree_file.name
    shutil.copy(locus_file, seq_dst)
    shutil.copy(tree_file, tree_dst)

    ctl_content = generate_ctl_content(
        seqfile=seq_dst.name,
        treefile=tree_dst.name,
        outfile=f"{locus}_{model_name}.txt",
        model_config=MODEL_CONFIGS[model_name],
    )
    ctl_path = run_dir / f"{locus}_{model_name}.ctl"
    ctl_path.write_text(ctl_content)

    try:
        result = subprocess.run(
            [codeml_bin, ctl_path.name],
            cwd=run_dir,
            capture_output=True,
            text=True,
            timeout=3600,
        )
        ok = result.returncode == 0
        msg = "" if ok else result.stderr[-500:]
    except subprocess.TimeoutExpired:
        ok, msg = False, "timeout"
    except FileNotFoundError:
        ok, msg = False, f"codeml binary not found: {codeml_bin}"

    return locus, model_name, ok, msg


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-dir", required=True, type=Path,
                         help="Directory with one codon alignment per locus.")
    parser.add_argument("--tree", required=True, type=Path,
                         help="Fixed Newick tree used for all loci (IQ-TREE topology).")
    parser.add_argument("--output-dir", required=True, type=Path,
                         help="Directory where per-locus/model results are written.")
    parser.add_argument("--models", nargs="+", default=list(MODEL_CONFIGS.keys()),
                         choices=list(MODEL_CONFIGS.keys()),
                         help="Models to run for every locus (default: all).")
    parser.add_argument("--codeml-bin", default="codeml",
                         help="Path to the codeml executable (default: 'codeml' on PATH).")
    parser.add_argument("--threads", type=int, default=1,
                         help="Number of loci/model runs to execute in parallel.")
    args = parser.parse_args()

    loci = find_alignments(args.input_dir)
    if not loci:
        raise SystemExit(f"No alignment files found in {args.input_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    jobs = [(f, args.tree, m) for f in loci for m in args.models]

    print(f"{len(loci)} loci x {len(args.models)} models = {len(jobs)} codeml runs")

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as pool:
        futures = [pool.submit(run_one, f, t, m, args.output_dir, args.codeml_bin)
                   for f, t, m in jobs]
        for fut in concurrent.futures.as_completed(futures):
            locus, model_name, ok, msg = fut.result()
            status = "OK" if ok else f"FAILED ({msg})"
            print(f"[{locus}] {model_name}: {status}")


if __name__ == "__main__":
    main()
