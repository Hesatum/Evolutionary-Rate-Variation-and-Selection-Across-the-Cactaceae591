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
01_exon_dNdS/       Exon selection analysis (dN/dS, codon models)
  codeml_configs/      PAML/codeml control files (M0, M1a, M2a, M7, M8, branch model)
  (add) reading_frame_validation.py   6-frame translation + internal stop-codon screen
  (add) parse_codeml_results.py       parses codeml output, LRTs, BH correction, BEB sites

02_intron_phyloP/   Intron evolutionary-rate analysis (RPhast phyloP)
  (add) phylopv4_linux_filtered.R          main pipeline: QC, neutral model, per-site phyloP
  (add) phylopv4_linux_filtered_byfeature.R  whole-locus (per-feature) phyloP variant;
                                             exploratory only — NOT used for the published results
  plot_figure2_representative_introns.py   builds Figure 2 (4 representative intron profiles)
                                             from phylopv4_linux_filtered.R's per-site output;
                                             also exports a submission-ready TIFF (see below)

03_species_tree/    Species-tree estimation
  (add) IQ-TREE commands/partition files (ModelFinder Plus, 1000 UFBoot)
  (add) sensitivity re-estimation after excluding fast-evolving loci

data/               Small derived tables only (e.g. Table S1/S2 sources).
                     Do NOT put raw reads or full alignments here — link to
                     SRA/BioProject instead (see manuscript Table 1).
```

## Status

This is a skeleton being populated. Scripts are currently split across two
machines; folders above marked "(add)" are placeholders for files not yet
copied in.

## Software versions used

Fill in as scripts are added (see manuscript Methods 2.2–2.3 for the full
list): fastp v0.20.1, HybPiper v2.1.3, MAFFT, trimAL v1.4.rev22, AliView 1.27,
MACSE v2, IQ-TREE 2.0.7 + ModelFinder Plus, PAML 4.9j (codeml), R + RPhast
(phyloP), Python 3.x.

## How this maps to the manuscript

| Manuscript section | Folder |
|---|---|
| 2.2 Data Filtering and Processing (reading-frame validation) | `01_exon_dNdS/` |
| 2.3 Signatures of selection — site/branch models, BH correction | `01_exon_dNdS/` |
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
