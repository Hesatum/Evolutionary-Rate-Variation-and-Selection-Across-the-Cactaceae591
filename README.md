# Cactaceae591: Selection and Evolutionary-Rate Analysis Scripts

Analysis code for:

> Vieira da Silva M, Moraes EM, Romeiro Brito M, Franco FF. *Evolutionary Rate
> Variation and Selection across Cactaceae591 nuclear loci in Cereus
> (Cactaceae).* (in preparation)

This repository provides the scripts used to detect signatures of selection
in exonic regions (dN/dS, PAML/codeml) and evolutionary-rate acceleration in
intronic regions (phyloP/RPhast) of the Cactaceae591 target-capture panel,
using *Cereus* (Cactaceae) as the study genus. Raw target-capture reads are
available on NCBI SRA (accession numbers in Table 1 of the manuscript); this
repository does not redistribute raw sequencing data.

## Repository structure

```
01_exon_dNdS/                    Exon selection analysis (dN/dS, codon models)
  reading_frame_validation.py      6-frame translation + internal stop-codon screen (2.2)
  frame_validation_blastx.py       homology-based frame check via NCBI blastx (2.2)
  run_codeml_batch_default.py      batches codeml across loci with default .ctl settings (2.3)
  parse_codeml_results.py          LRT + Benjamini-Hochberg FDR + BEB sites (2.3)
  figure3_lrt_density.py           Figure 3 and Figure S2 (LRT density plots)
  Figure_S2_M7_M8_LRT_density.png/.pdf   Figure S2
  codeml_configs/                  example .ctl files for M0, M1a, M2a, M7, M8, Branch

02_intron_phyloP/                Intron evolutionary-rate analysis (RPhast phyloP)
  phylopv4_linux_filtered.R        main pipeline: QC, neutral model, per-site phyloP
  phylopv4_linux_filtered_byfeature.R  whole-locus (per-feature) phyloP variant;
                                      not used for the published results
  plot_figure2_representative_introns.py   Figure 2 (four representative intron profiles)
  plot_individual_intron_profiles.py   one plot per intron; defaults to the 17
                                      accelerated introns (exclude_lists/introns_accelerated_17.csv)
  individual_profiles/             17 PNGs, one per accelerated intron
  build_figure_s1.py               composes individual_profiles/ into Figure S1 (panels A-Q)
  Figure_S1_accelerated_intron_profiles.png/.pdf   Figure S1

03_species_tree/                 Species-tree estimation and sensitivity re-runs
  build_reduced_supermatrix.py     removes chosen orthogroups from the full supermatrix + partitions
  run_iqtree.sh                    IQ-TREE 2 command (ModelFinder Plus, 1000 UFBoot)
  plot_tree.R                      draws a tree in the Figure 1/S2 style, single or mirrored
  taxa_codes.txt                   sample codes of the 18 final taxa (Table 1)
  exclude_lists/                   orthogroups removed in each sensitivity tree (see below)
  trees/full_18tax/                full dataset (no exclusion); Figure 1
  trees/<name>/                    each sensitivity tree
  Figure_S3_species_tree_comparison.png/.pdf/.eps   Figure S3: trees/full_18tax/ vs.
                                      trees/exons_and_introns_removed/

data/                            Supplementary tables
  Table_S1_codeml_results.xlsx     Table S1: site- and branch-model codeml results for the
                                      94 analyzed exons (one sheet per comparison)
  Table_S1_README.txt              column definitions for Table S1
  Table_S2_summary_by_orthogroup.xlsx  Table S2: phyloP scores (LRT, CONACC) for the 70
                                      analyzed introns, summarized per intron x clade (490 rows)
  Table_S2_per_site_phyloP_scores.csv.gz  per-site phyloP scores behind Table S2 (515,739 rows)
  Table_S2_README.txt              column definitions and thresholds for Table S2
```

## codeml and IQ-TREE models

**Site models (per-codon selection, tested on the whole tree)**
- `M1a` (nearly neutral) vs. `M2a` (adds a class of codons with ω > 1),
  LRT with 2 degrees of freedom.
- `M7` (β-distributed ω < 1) vs. `M8` (β plus a class with ω > 1), LRT with
  2 degrees of freedom; more sensitive to weak signals than M1a vs. M2a.
- Both comparisons are run for every exon and p-values are corrected with
  the Benjamini-Hochberg procedure. Exons significant under both comparisons
  are listed in Table 2 of the manuscript. Positively selected sites are
  identified by Bayes Empirical Bayes (BEB), Pr(ω > 1) > 0.95.

**Branch model (lineage-specific selection)**
- `M0` (one ω for the whole tree) vs. a `Branch` model in which each major
  *Cereus* clade (A1, A2, B, C, D, E) and the outgroup has its own ω
  (7 ω categories, LRT with 6 degrees of freedom).
- Branch-site models are not used in this study.

**IQ-TREE 2 (species tree)**
- The supermatrix is partitioned by alignment (568 partitions); ModelFinder
  Plus selects the substitution model of each partition (BIC).
- 1000 ultrafast bootstrap (UFBoot2) replicates give branch support.
- The ML topology is the fixed input tree for every codeml run; codeml
  re-estimates branch lengths and ω, not the topology.

## phyloP model and workflow

**Neutral model**
- A single neutral model is fit with phyloFit on a concatenated alignment of
  a subset of introns (selected by nucleotide diversity and alignment
  length). The script tries HKY85, then REV, then JC69, and keeps the first
  that converges. Every phyloP test uses this model as its null.

**Per-site LRT (`method = "LRT"`, `mode = "CONACC"`)**
- For each clade (A1, A2, B, C, D, E, Outgroup), phyloP tests the branch
  leading to the clade's ancestor against the neutral model, once per
  alignment column.
- The score is the signed `-log10(P)`: positive for conservation, negative
  for acceleration. `|score| > 1.3` corresponds to `P < 0.05`, without
  site-level multiple-testing correction. An intron is classified as
  accelerated only when accelerated sites form sustained blocks (Results 3.2).
- `phylopv4_linux_filtered_byfeature.R` tests each whole intron as a single
  feature (one LRT per intron per clade, phyloP `features=` argument). It
  was not used for any result reported in the manuscript.

## Sensitivity trees (Figure S3)

`03_species_tree/build_reduced_supermatrix.py` removes every partition of
the orthogroups in an exclusion list from the full supermatrix and partition
file. Three trees were estimated, one per list:

1. **Exons under positive selection removed** (21 orthogroups, Table 2):
   `exclude_lists/exons_under_selection_21.csv`.
2. **Introns with sustained clade-specific acceleration removed** (17
   orthogroups): `exclude_lists/introns_accelerated_17.csv`.
3. **Both removed** (Figure S3):
   `exclude_lists/exons_and_introns_combined_38.csv`.

`exclude_lists/empty_none.csv` is an empty list that reproduces the full
dataset through the same script.

The intron list was defined by visual inspection of the per-site phyloP
profile of each of the 70 introns (sustained block vs. scattered sites).

### Results

All trees were estimated on the 18 final taxa from the 568-partition
supermatrix (IQ-TREE 2.2.2.6, seed 141309):

| Tree | Partitions kept | Sites | UFBoot support | RF vs. full tree |
|---|---|---|---|---|
| Full dataset (`trees/full_18tax/`, Figure 1) | 568/568 | 585,254 | 100% at every node | - |
| 1. Exons removed | 534/568 | 541,816 | 100% at every node | 0 |
| 2. Introns removed | 559/568 | 575,781 | 100% at every node | 0 |
| 3. Both removed | 525/568 | 532,343 | 100% at every node | 0 |

### Plotting

    Rscript plot_tree.R single <tree.contree_or_.nwk> <output_prefix> "<panel title>"
    Rscript plot_tree.R mirror <left.contree_or_.nwk> <right.contree_or_.nwk> <prefix> "Left" "Right"

`plot_tree.R` roots the tree on the *Cipocereus* outgroup
(`ape::root(..., resolve.root = TRUE)`), orders tips as in Figure 1, and
labels clades A1/A2/B/C/D/E/Outgroup. The support of the root bipartition
(*Cipocereus* | *Cereus*) is shown on both of its nodes. Rooting assigns zero
length to the outgroup branch at the root; that edge is drawn with a short
display-only length (`ROOT_STUB_FRACTION`) and its branch-length label is
left blank.

## Multiple-testing correction

EasyPAML extracts lnL, ω and BEB sites and computes raw LRT p-values, but
does not correct for multiple testing. `parse_codeml_results.py` applies the
Benjamini-Hochberg correction reported in the manuscript (significant exons
before and after correction: 25 and 21 for M1a vs. M2a, 31 and 28 for M7 vs.
M8, 2 and 1 for the branch model).

## codeml batch runs

All codeml runs were launched through EasyPAML
(https://github.com/LAGEVOL/EasyPAML). `run_codeml_batch_default.py`
reproduces its batch execution and default `.ctl` parameters for the models
used here (M0, M1a, M2a, M7, M8, Branch).

## Software versions

fastp v0.20.1, HybPiper v2.1.3, MAFFT, trimAl v1.4.rev22, AliView 1.27,
MACSE v2, IQ-TREE 2.2.2.6 + ModelFinder Plus, PAML 4.9j (codeml), R 4.4.1
(ape 5.8.1, RPhast), Python 3.

## How this maps to the manuscript

| Manuscript section | Location |
|---|---|
| 2.2 Data filtering and processing (reading-frame validation, Blastx check) | `01_exon_dNdS/` |
| 2.3 Site and branch models, BH correction, BEB sites | `01_exon_dNdS/` |
| 2.3 Species-tree estimation and sensitivity re-estimation | `03_species_tree/` |
| 2.3 Evolutionary rates in introns (phyloP/RPhast) | `02_intron_phyloP/` |
| Figure 1 | `03_species_tree/trees/full_18tax/` |
| Figure 2 | `02_intron_phyloP/plot_figure2_representative_introns.py` |
| Figure 3, Figure S2 | `01_exon_dNdS/figure3_lrt_density.py` |
| Figure S1 | `02_intron_phyloP/` |
| Figure S3 | `03_species_tree/` |
| Table S1, Table S2 | `data/` |

## License

Code released under the MIT License (see `LICENSE`). This covers the analysis
scripts only; sequence data are governed by their original NCBI SRA terms.

## Citation

If you use these scripts, please cite the manuscript above. A versioned
release of this repository is archived on Zenodo:
[10.5281/zenodo.22888231](https://doi.org/10.5281/zenodo.22888231).
