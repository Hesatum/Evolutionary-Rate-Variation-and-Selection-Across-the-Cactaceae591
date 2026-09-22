Table S1 -- full codeml site- and branch-model results for the 94 analyzed exons
==================================================================================

FILE
----
Table_S1_codeml_results.xlsx
    One row per exon (94 rows in every comparison sheet).

SHEETS
------
Summary
    Genes analyzed, and counts/percentages significant at raw P < 0.05
    and P < 0.01, for each of the four LRT comparisons below. Matches
    the manuscript's reported raw significance counts (M1a vs M2a: 25/94;
    M7 vs M8: 31/94; M0 vs Branch: 2/94).

M0 vs M1a
    One-ratio model (M0) vs. nearly neutral model (M1a). Not used for a
    significance claim in the manuscript; included for completeness
    since it is part of the standard PAML site-model workflow.

M1a vs M2a
    Nearly neutral (M1a, null) vs. positive selection (M2a, alternative).
    LRT, df = 2. "Positive Sites (BEB)" lists Bayes Empirical Bayes sites
    with posterior probability > 0.95 (*) or > 0.99 (**), each with the
    site's posterior mean omega +/- standard error.

M7 vs M8
    Beta distribution (M7, null) vs. beta + positive-selection class (M8,
    alternative). Same LRT/BEB conventions as M1a vs M2a.

M0 vs Branch
    One-ratio model (M0, null) vs. free-ratio branch model (Branch,
    alternative) with every major Cereus clade (A1, A2, B, C, D, E,
    Outgroup) labeled separately. LRT, df = 6 (7 omega categories minus
    1). One exon (OG0070800) has a negative 2-delta-l, an impossible
    value for nested models that signals a codeml optimisation failure;
    it is excluded from Figure 3B's density plot (see the figure legend)
    but kept as a row in this table.

Every sheet also carries the functional annotation columns from the
Cactaceae591 reference panel (Romeiro-Brito et al. 2022): ortholog code,
Arabidopsis thaliana BLAST hit, PLAZA 5.0 annotation, and putative gene
function, consistent with Table 2's "Positively selected sites" summary.
The source data's "Angiosperm-353 loci ID" and "paralogs identification"
columns are not included here (out of scope for this manuscript).

SIGNIFICANCE COLUMNS
---------------------
"Sig. p<0.05" / "Sig. p<0.01" report the RAW likelihood-ratio p-value
against the chi-squared critical value at the stated df, before multiple-
testing correction. The manuscript's Benjamini-Hochberg-corrected counts
(21/94 under M1a vs M2a, 28/94 under M7 vs M8, 1/94 under M0 vs Branch;
Methods 2.3, Results 3.3) are not a column in this table; they are
reproducible from the p-value column here with
01_exon_dNdS/parse_codeml_results.py.

Generated from the codeml batch results in resultados_easypaml_2026/
(EasyPAML batch output, one *_results.txt per exon per model; see
01_exon_dNdS/run_codeml_batch_default.py for the batch/default-ctl logic
behind these runs) after reading-frame validation
(01_exon_dNdS/reading_frame_validation.py,
01_exon_dNdS/frame_validation_blastx.py).
