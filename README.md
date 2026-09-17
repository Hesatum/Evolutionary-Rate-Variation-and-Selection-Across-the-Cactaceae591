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
  (add) phylopv4_linux_filtered.R    main pipeline: QC, neutral model, per-site phyloP

03_species_tree/                 Species-tree estimation
  (add) IQ-TREE commands/partition files (ModelFinder Plus, 1000 UFBoot)
  (add) sensitivity re-estimation after excluding fast-evolving loci

data/                            Small derived tables only (e.g. Table S1/S2 sources).
                                  Do NOT put raw reads or full alignments here — link to
                                  SRA/BioProject instead (see manuscript Table 1).
```

## Status

`01_exon_dNdS/` is done and checked against the manuscript's Methods section
and figure legends. `02_intron_phyloP/` and `03_species_tree/` are still
placeholders. The phyloP scripts are on a second machine and haven't been
copied in yet.

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
| 2.3 Evolutionary rates in introns (phyloP/RPhast) | `02_intron_phyloP/` |
| Table S1 (codeml results), Table S2 (phyloP scores) | `data/` |

## License

Code released under the MIT License (see `LICENSE`). This covers the analysis
scripts only — sequence data are governed by their original NCBI SRA terms.

## Citation

If you use these scripts, please cite the manuscript above. A versioned,
citable release of this repository is archived on Zenodo: **DOI: (add after
first release)**.
