Table S2 -- phyloP evolutionary-rate scores for the 70 analyzed introns
========================================================================

FILES
-----
Table_S2_phyloP_results.csv.gz
    Per-site phyloP scores (LRT method, CONACC mode). One row per
    alignment position x clade x orthogroup.
    Columns:
      Orthogroup              Orthogroup ID (OG code, Romeiro-Brito et al. 2022)
      Clade                   One of A1, A2, B, C, D, E, Outgroup
      Alignment_Position_bp   1-based column index within the intron alignment
      PhyloP_Score            Signed -log10(P); positive = conservation,
                               negative = acceleration relative to the
                               genome-wide neutral model
      P_value                 Two-tailed P-value corresponding to PhyloP_Score

    515,739 rows. Score/P-value are rounded (3 decimals / 3 significant
    figures) from the raw pipeline output for file size; this does not
    change any threshold comparison used in the manuscript.

    |PhyloP_Score| > 1.3 corresponds to raw P < 0.05 (uncorrected at the
    site level; see manuscript Methods 2.3 and Results 3.2 for why a
    per-site FDR is not the criterion used; sustained blocks of
    accelerated sites are).

Table_S2b_summary_by_orthogroup.csv
    One row per orthogroup x clade (490 rows), summarizing the file above:
    number of sites, mean score, count and percentage of sites reaching
    the conservation/acceleration threshold. Useful for a quick overview
    without loading the full per-site table.

METHOD
------
GTR substitution model (HKY85, chosen automatically; see Methods 2.3),
fit once on a concatenated alignment of representative introns and reused
as the neutral background for every clade test. Each clade's ancestral
branch is tested independently against this fixed neutral model with a
likelihood-ratio test, per alignment column.

Generated from phylopv4_linux_filtered.R.
