# Cactaceae591 — Selection and Evolutionary-Rate Analysis Scripts

Analysis code for:

> Vieira da Silva M, Moraes EM, Romeiro Brito M, Franco FF. *Evolutionary Rate
> Variation and Selection Across the Cactaceae591 Probe Set.* (in preparation)

This repository accompanies the manuscript and provides the scripts used to
detect signatures of selection in exonic regions (dN/dS, PAML/codeml) and
evolutionary-rate acceleration in intronic regions (phyloP/RPhast) of the
Cactaceae591 target-capture panel, using *Cereus* (Cactaceae) as the study
genus. Raw target-capture reads are publicly available on NCBI SRA (see
manuscript Table 1 / Data Availability); this repository does not redistribute
raw sequencing data.

## Repository structure

```
01_exon_dNdS/                    Exon selection analysis (dN/dS, codon models)
  reading_frame_validation.py      6-frame translation + internal stop-codon screen (2.2)
  frame_validation_blastx.py       homology-based frame check via NCBI blastx (2.2)
  run_codeml_batch_default.py      batches codeml across loci with default .ctl settings (2.3)
  parse_codeml_results.py          LRT + Benjamini-Hochberg FDR + BEB sites (2.3)
  figure3_lrt_density.py           Figure 3 (LRT density plots)
  codeml_configs/                  example .ctl files for M0, M1a, M2a, M7, M8, Branch

02_intron_phyloP/                Intron evolutionary-rate analysis (RPhast phyloP)
  phylopv4_linux_filtered.R        main pipeline: QC, neutral model, per-site phyloP
  phylopv4_linux_filtered_byfeature.R  whole-locus (per-feature) phyloP variant;
                                      exploratory only — NOT used for the published results
  plot_figure2_representative_introns.py   builds Figure 2 (4 representative intron
                                      profiles) from phylopv4_linux_filtered.R's per-site
                                      output; also exports a submission-ready TIFF

03_species_tree/                 Species-tree estimation and sensitivity re-runs
  build_reduced_supermatrix.py     removes chosen orthogroups from the full supermatrix + partitions
  run_iqtree.sh                    IQ-TREE 2 command (ModelFinder Plus, 1000 UFBoot)
  plot_tree.R                      draws a tree in the manuscript's Figure 1/S2 style, single or mirrored
  taxa_codes.txt                   the 18 final taxa (Table 1), in figure order
  exclude_lists/                   which orthogroups define each sensitivity tree (see below)
  trees/exon_selection_removed/    tree 1 output: .contree/.treefile, build log, figure (done)

data/                            Small derived tables only (e.g. Table S1/S2 sources).
                                  Do NOT put raw reads or full alignments here — link to
                                  SRA/BioProject instead (see manuscript Table 1).
  08_statistics_by_gene.csv        phyloP results summarized per intron x clade (70 x 7 =
                                      490 rows: n_sites, mean_score, %conserved/%accelerated).
                                      Derived from phylopv4_linux_filtered.R's per-site output
                                      (06_all_phylop_results.csv, 208 MB, kept out of git —
                                      this is the small summary, not the full per-site table
                                      submitted to the journal as Table S2).
```

## codeml and IQ-TREE: models used in this study

This study uses a subset of what PAML/codeml and IQ-TREE can do. To keep the
repository self-explanatory for readers who don't work with these tools day
to day:

**Site models (per-codon selection, tested on the whole tree at once)**
- `M1a` (nearly neutral) vs. `M2a` (adds a class of codons with ω > 1): a
  nested comparison, LRT with 2 degrees of freedom.
- `M7` (β-distributed ω, all values < 1) vs. `M8` (same, plus a class with
  ω > 1): the same comparison with a continuous null instead of a discrete
  one, more sensitive to weak signals of positive selection than M2a.
- Both pairs are run for every exon; a locus counts as under positive
  selection only when both M2a and M8 reject their null after
  Benjamini-Hochberg correction (manuscript Table 2). Sites are called
  positively selected from the Bayes Empirical Bayes (BEB) posterior,
  threshold Pr(ω > 1) > 0.95.

**Branch model (lineage-specific selection)**
- `M0` (one ω for the whole tree) vs. a free-ratio `Branch` model where every
  major *Cereus* clade (A1, A2, B, C, D, E, plus the outgroup) is labeled
  separately in the tree and gets its own ω. This is not a single
  foreground-vs-background branch test: all seven clades are marked at once,
  which is why the LRT has 6 degrees of freedom (7 ω categories minus 1).
- Branch-site models are supported by codeml and by EasyPAML but are **not**
  used in this manuscript, so they're outside the scope of this repository.

**IQ-TREE 2 (species tree used as fixed input to every codeml run)**
- ModelFinder Plus picks the nucleotide substitution model for the
  concatenated supercontig alignment by BIC, instead of assuming one a
  priori.
- 1000 ultrafast bootstrap (UFBoot2) replicates give branch support.
- The resulting topology is fixed and reused, unmodified, as the tree input
  for every codeml run above: codeml does not re-estimate topology, only
  branch lengths and ω under each model.

## phyloP: model and workflow used in this study

The intron side uses RPhast's `phyloP()` quite differently from how codeml is
used above — one neutral model is fit once, then every clade is tested
against it, rather than fitting a fresh model per hypothesis. To keep this
self-explanatory:

**Neutral model (fit once, shared by every test)**
- A single background substitution model (GTR, tried at HKY85 → REV → JC69 in
  that order of preference, first one that fits) is fit on a *concatenated*
  alignment built from a subset of the introns (the most diverse and the
  longest ones), not on each intron individually. Every subsequent phyloP
  call reuses this one fitted model as its null.

**Per-site LRT (`method = "LRT"`, `mode = "CONACC"`)**
- For each of the 7 clades (A1, A2, B, C, D, E, Outgroup) in turn, phyloP
  rescales the single branch leading to that clade's ancestor and compares
  the likelihood against the fixed neutral model; a likelihood-ratio test
  run independently *per alignment column*, not once per whole intron.
- The signed score is `-log10(P)`, positive for conservation (slower than
  neutral) and negative for acceleration (faster than neutral). `|score| >
  1.3` corresponds to raw `P < 0.05`; this is **not** Benjamini-Hochberg
  corrected at the site level. The manuscript's actual safeguard against
  false positives from testing thousands of correlated sites is requiring
  **sustained blocks** of accelerated sites (Results 3.2), not a per-site FDR
  — a single accelerated site in isolation is not treated as evidence of
  anything.
- `phylopv4_linux_filtered_byfeature.R` explores an alternative: testing each
  whole intron as one feature (one LRT per intron per clade, via phyloP's
  `features=` argument) instead of per site, so a locus-level
  Benjamini-Hochberg correction could be applied without the site-level
  multiple-testing problem. It was a useful diagnostic; the per-site test
  has essentially no power to detect conservation column-by-column even
  where the whole-locus test finds it clearly — but it uses a materially
  different statistical unit (490 whole-locus tests vs. ~500k per-site
  tests) and gives very different-looking numbers. **It was not used for any
  number reported in the manuscript**; it's kept here only so the reasoning
  behind the per-site + sustained-block approach is reproducible, not just
  asserted.

## Sensitivity trees (Figure S2)

`03_species_tree/build_reduced_supermatrix.py` takes the full supermatrix
and partition file, drops every partition belonging to a chosen list of
orthogroups, and writes a reduced supermatrix ready for `run_iqtree.sh`.
Three trees are relevant to this study, run the same way with three
different `--exclude-list` files:

1. **Exons under positive selection removed** (21 orthogroups, Table 2):
   `exclude_lists/exons_under_selection_21.csv`.
2. **Introns with sustained clade-specific acceleration removed** (Table
   S2): `exclude_lists/introns_accelerated_17.csv`.
3. **Both removed together**, the tree actually reported in the manuscript
   (Figure S2): run with a list that concatenates 1 and 2.

The intron list has no automated formula behind it — "sustained block vs.
scattered sites" (manuscript section 3.2) was classified by looking at each
of the 70 retained introns' per-site phyloP profile. Two automated proxies
were tried first (longest contiguous accelerated run; overall % of sites
accelerated) and both failed to recover even the loci already named in the
manuscript text — e.g. OG0074859 has only 1.4% accelerated sites and never
ranks in the top 20 by either measure, yet visually shows a short, clean
accelerated block. This confirms the classification is a qualitative call,
not a threshold on either statistic.

`exclude_lists/introns_accelerated_17.csv` is the human-verified list (M.
Vieira da Silva, reviewed against per-locus plots): **17 orthogroups**, not
the 10 stated in the current manuscript draft — OG0074859 was reconsidered
and dropped on review, and the borderline cases initially set aside as
ambiguous were, on reflection, included rather than excluded. **The
manuscript text (Results 3.2, Discussion, and the 46-partition /
522-locus arithmetic tied to it) has not been updated to match yet** — this
is the single open item before the sensitivity analysis and the prose agree
with each other.

Tree 1 has been run end to end against the real 568-partition supermatrix
(`trees/exon_selection_removed/`): 534 of 568 partitions kept, 541,816
sites, IQ-TREE 2.2.2.6, seed 141309. Every node reaches 100% ultrafast
bootstrap support, including the two clades that sat at 92% and 98% in the
full tree, and the topology is identical to the full tree (unrooted
Robinson-Foulds distance = 0, checked with `ape::dist.topo`). Trees 2 and 3
are not run yet: `supermatrix_full.fasta` and `partitions_full.txt` (the
568-partition inputs `build_reduced_supermatrix.py` needs) aren't on this
machine — only the already-reduced 534-partition output of tree 1 is.

Every tree, whichever exclusion list produced it, is plotted the same way:

    Rscript plot_tree.R single <tree.contree> <output_prefix> "<panel title>"
    Rscript plot_tree.R mirror <full_tree.contree> <reduced_tree.contree> <prefix> "Full" "Reduced"

`plot_tree.R` roots on the *Cipocereus* outgroup, orders tips to match the
manuscript's Figure 1, and labels clades A1/A2/B/C/D/E/Outgroup; it works on
any Newick/`.contree` file with the same 18 tips, not just the three trees
above.

## Status

All three folders are populated and checked against the manuscript's Methods
section and figure legends. Two open items remain, both in "Sensitivity
trees" above: the intron exclusion list is finalized at 17 orthogroups
(not the 10 currently stated in the manuscript — text update pending), and
trees 2/3 can't be regenerated on this machine until `supermatrix_full.fasta`
/ `partitions_full.txt` (the 568-partition inputs) are located or copied in.

### A note on `parse_codeml_results.py`

EasyPAML's batch-analysis backend extracts lnL, omega, and BEB sites, and it
computes raw LRT p-values, but it does not correct for multiple testing. The
Benjamini-Hochberg step reported in the manuscript (Methods 2.3, Table 2) was
run by hand during manuscript revision, not from a saved script. The version
here reimplements that same BH procedure as a documented step, so the
published numbers (25 to 21 significant exons under M1a vs. M2a, 31 to 28
under M7 vs. M8, 2 to 1 under the branch model) can be reproduced from the
raw codeml output.

### A note on codeml batch runs

Every codeml run reported in the manuscript was launched through EasyPAML, a
batch GUI for PAML/codeml. Rather than including that whole application here,
`run_codeml_batch_default.py` pulls out the batch-execution logic and the
exact default `.ctl` parameters EasyPAML uses for the models this study
reports (M0, M1a, M2a, M7, M8, Branch), so the codeml settings behind the
manuscript's results can be inspected and rerun without the GUI.

## Software versions used

fastp v0.20.1, HybPiper v2.1.3, MAFFT, trimAL v1.4.rev22, AliView 1.27,
MACSE v2, IQ-TREE 2.0.7 + ModelFinder Plus, PAML 4.9j (codeml), R + RPhast
(phyloP), Python 3.x.

## How this maps to the manuscript

| Manuscript section | Folder |
|---|---|
| 2.2 Data Filtering and Processing (reading-frame validation, Blastx check) | `01_exon_dNdS/` |
| 2.3 Signatures of selection — site/branch models, BH correction, BEB sites | `01_exon_dNdS/` |
| Figure 3 (LRT density) | `01_exon_dNdS/figure3_lrt_density.py` |
| 2.3 Species-tree estimation (IQ-TREE, sensitivity re-estimation) | `03_species_tree/` |
| 2.3 Evolutionary rates in introns (phyloP/RPhast) — see "phyloP: model and workflow" above | `02_intron_phyloP/` |
| Figure 2 (representative intron profiles) | `02_intron_phyloP/plot_figure2_representative_introns.py` |
| Table S1 (codeml results), Table S2 (phyloP scores) | `data/` |

## License

Code released under the MIT License (see `LICENSE`). This covers the analysis
scripts only — sequence data are governed by their original NCBI SRA terms.

## Citation

If you use these scripts, please cite the manuscript above. A versioned,
citable release of this repository is archived on Zenodo: **DOI: (add after
first release)**.
