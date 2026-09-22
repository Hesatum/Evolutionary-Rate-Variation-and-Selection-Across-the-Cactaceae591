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
  plot_individual_intron_profiles.py   same drawing conventions, one full-panel plot per
                                      orthogroup instead of 4 in a grid; defaults to the 17
                                      accelerated introns (exclude_lists/introns_accelerated_17.csv)
  individual_profiles/             its output: 17 PNGs, one per accelerated intron
  build_figure_s1.py                composes individual_profiles/ into the final Figure S1
                                      (one numbered file, panels A-Q, BJLS-style: <2 MB, see
                                      https://academic.oup.com/biolinnean/pages/General_Instructions)
  Figure_S1_accelerated_intron_profiles.png/.pdf   the assembled Figure S1 (done)

03_species_tree/                 Species-tree estimation and sensitivity re-runs
  build_reduced_supermatrix.py     removes chosen orthogroups from the full supermatrix + partitions
  run_iqtree.sh                    IQ-TREE 2 command (ModelFinder Plus, 1000 UFBoot)
  prune_to_final_taxa.R            prunes a 49-taxon tree to the 18 in Table 1 -- the actual
                                      method behind Figure 1/S2 (see "Two different build
                                      methods" below); NOT the same as building fresh on 18 taxa
  plot_tree.R                      draws a tree in the manuscript's Figure 1/S2 style, single or
                                      mirrored, with a cosmetic root stub (ROOT_STUB_FRACTION)
  taxa_codes.txt                   the 18 final taxa (Table 1), in figure order
  taxa_codes_49.txt                all 49 sample headers, for a 49-taxon build (prune afterwards)
  exclude_lists/                   which orthogroups define each sensitivity tree (see below)
  trees/<name>/                    tree built directly on 18 taxa (superseded, see below)
  trees/<name>/49tax/              same exclusion list, built on all 49 taxa
  trees/<name>/pruned_18tax/       .../49tax/ pruned to the 18 final taxa -- the correct one
  trees/full_49tax_reference/      the original, unmodified 49-taxon tree and its pruning,
                                      confirming the pruning method reproduces Figure 1 exactly
  Figure_S2_species_tree_comparison.png/.pdf/.eps   the assembled Figure S2 (done): plot_tree.R
                                      mirror mode on the two pruned_18tax/ trees above

data/                            Small derived tables only (e.g. Table S1/S2 sources).
                                  Do NOT put raw reads or full alignments here — link to
                                  SRA/BioProject instead (see manuscript Table 1).
  Table_S1_codeml_results.xlsx     Table S1: full site- and branch-model codeml results for
                                      the 94 analyzed exons (one row/exon per comparison
                                      sheet), translated to English and cross-checked against
                                      the manuscript's raw significance counts (25/94, 31/94,
                                      2/94) — see Table_S1_README.txt.
  Table_S1_README.txt              sheet-by-sheet column definitions for Table S1.
  Table_S2_phyloP_results.csv.gz   Table S2: per-site phyloP scores (LRT only, CONACC),
                                      515,739 rows, English column names, rounded for size.
                                      2.5 MB gzipped, down from the 208 MB raw pipeline output
                                      (06_all_phylop_results.csv, both methods + internal
                                      columns, kept out of git) — see Table_S2_README.txt.
  Table_S2b_summary_by_orthogroup.csv  same data summarized per intron x clade (490 rows).
  Table_S2_README.txt              column definitions and thresholds, for co-authors/readers
                                      who don't want to re-derive them from Methods 2.3.
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

### Two different build methods give two different answers

All three trees were first run end to end against the real 568-partition
supermatrix, **built directly on the 18 final taxa** (IQ-TREE 2.2.2.6, seed
141309, `trees/<name>/` without a `49tax`/`pruned_18tax` subfolder):

| Tree | Partitions kept | Sites | UFBoot support | RF vs. full tree |
|---|---|---|---|---|
| 1. Exons removed | 534/568 | 541,816 | 100% at every node | 0 |
| 2. Introns removed | 559/568 | 575,781 | 100% at every node | 0 |
| 3. Both removed | 525/568 | 532,343 | 100% at every node | 0 |

Every one of them resolves the two clades that sit at 92% and 98% support in
the full tree up to 100%, which is exactly what the manuscript's Discussion
claims for the loci-removed tree. But this is **not how the manuscript's
species tree was actually built**, and the two methods do not agree.

`partitions.txt.contree` (the source of `FINAL_figure1_source_rooted_supported.nwk`,
copied into `trees/full_49tax_reference/49tax/` as
`run_full_49tax_original.*`) is a tree of **49 samples**: every accession
available for the genus, not just the 18 in Table 1, built with IQ-TREE
2.0.7, seed 107775. The manuscript's actual Figure 1 is that tree **pruned**
down to the 18 final taxa with `ape::keep.tip()`, which keeps whatever
support the 49-taxon analysis assigned to each retained split rather than
re-estimating it from 18 taxa. `prune_to_final_taxa.R` reproduces this
exactly: pruning `run_full_49tax_original.contree` to the 18 codes in
`taxa_codes.txt` gives RF = 0 and the identical support distribution (14
nodes at 100%, one at 92%, one at 98%) as the published figure.

Rerunning the same prune-from-49 procedure on each loci-removed supermatrix
(49-taxon build, same 107775 seed, `trees/<name>/49tax/` and
`.../pruned_18tax/`) gives a different result:

| Tree (pruned from 49 taxa) | RF vs. published Figure 1 | Support at the two contested nodes |
|---|---|---|
| Full 568 loci, no exclusion (control) | 0 | 92%, 98% (identical to published) |
| 1. 21 exons removed | 0 | 92%, 98% (unchanged) |
| 2. 9 introns removed | 0 | **89%** (worse), one other node 100%→99% |
| 3. Both removed (30 OGs) | 0 | 92% (unchanged), one other node 100%→99% |

None of the three loci-removed trees reach 100% at those two nodes once
they're built the way Figure 1 actually was; the "full 568 loci, no
exclusion" row confirms it isn't about which loci are in the alignment at
all, since pruning that unmodified dataset from 49 taxa reproduces the
published 92%/98% exactly. The 100%-everywhere result in the first table
instead traces to dropping from 49 taxa to 18 **before** running IQ-TREE: a
direct 18-taxon build of the full 568-locus alignment, with nothing
excluded, also comes out 100% at every node (not included as its own row
above, since it's the same procedure as trees 1-3, just with an empty
exclusion list), so fewer, more divergent samples resolve those branches
regardless of which loci are used. Topology is unaffected either way (RF =
0 throughout), but the specific claim in the manuscript's Discussion, that
removing loci under selection is what pushed those two branches to 100%
support, is not supported once the tree is built the way Figure 1 was. This
needs the authors' attention before submission, independent of anything
else in this repository.

**A rooting artifact, not new to this repo:** `plot_tree.R` roots each
unrooted IQ-TREE consensus on the *Cipocereus* outgroup with
`ape::root(..., resolve.root = TRUE)`, which always assigns zero length to
one of the two branches at the root (it has no basis to split that length
otherwise). In every tree checked here, including the original published
tree, that zero-length branch is the one leading to the outgroup clade, so
the root used to look like a trifurcation in the figure even though the
topology is strictly bifurcating (`ape::is.binary()` is `TRUE`, node count
is `Ntip - 1`). `plot_tree.R` now gives that one edge a small, fixed,
display-only length (`ROOT_STUB_FRACTION`, a fraction of tree depth) so the
two root branches are visually distinguishable; the branch-length text label
still reports the true value (0, printed blank, as before), so no plotted
number is altered, only the line moves. Set `ROOT_STUB_FRACTION <- 0` to go
back to the untouched original behavior.

Every tree, whichever exclusion list or build method produced it, is
plotted the same way:

    Rscript plot_tree.R single <tree.contree_or_.nwk> <output_prefix> "<panel title>"
    Rscript plot_tree.R mirror <left.contree_or_.nwk> <right.contree_or_.nwk> <prefix> "Left" "Right"

`plot_tree.R` roots on the *Cipocereus* outgroup, orders tips to match the
manuscript's Figure 1, and labels clades A1/A2/B/C/D/E/Outgroup; it works on
any Newick/`.contree` file with the same 18 tips. To reproduce a tree the
way Figure 1 was actually built, run `build_reduced_supermatrix.py` with
`taxa_codes_49.txt` (all sample headers in the source fasta, not just the
18 in `taxa_codes.txt`) so the exclusion list is applied without reducing
taxon sampling, run `run_iqtree.sh` on the resulting 49-taxon supermatrix,
then prune with `prune_to_final_taxa.R`.

## Status

All three folders are populated and checked against the manuscript's Methods
section and figure legends. The manuscript text now says 17 accelerated
introns throughout, matching this repository.

One open item remains, in the manuscript text rather than in this
repository: the Discussion still credits the loci-removed tree with
resolving the two sub-100%-support branches. Built the way Figure 1 was
actually built (pruned from the 49-sample analysis, not assembled fresh on
18 taxa), none of the three sensitivity trees do that; see "Two different
build methods give two different answers" above for the numbers and a
drop-in replacement paragraph.

Table S1 and Table S2 are both in `data/`, in English, with exon/intron
counts cross-checked against the manuscript's reported numbers (94 exons,
25/31/2 raw significant; 70 introns). Both supplementary figures are
assembled and under BJLS's 2 MB recommendation:
`02_intron_phyloP/Figure_S1_accelerated_intron_profiles.png`/`.pdf` (17
panels) and
`03_species_tree/Figure_S2_species_tree_comparison.png`/`.pdf`/`.eps`
(full dataset vs. 21 exons + 17 introns removed, both pruned from the
49-taxon analysis; same topology and the same support at every node,
consistent with the corrected Discussion paragraph above rather than the
manuscript's current text).

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
