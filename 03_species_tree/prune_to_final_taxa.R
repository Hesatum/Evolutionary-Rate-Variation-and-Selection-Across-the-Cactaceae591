#!/usr/bin/env Rscript
# ==============================================================================
# prune_to_final_taxa.R
#
# Prunes a tree built on the full sample set (49 accessions in this study)
# down to the 18 final taxa used throughout the manuscript (Table 1),
# keeping whatever support values and branch lengths the full-sample
# analysis assigned to the retained splits.
#
# This is the actual method behind the manuscript's species tree (Figure 1)
# and its sensitivity re-estimations (Figure S2): every IQ-TREE run in this
# study, including the loci-removed sensitivity trees, was built on all 49
# available Cereus/Cipocereus accessions and pruned afterwards, not built
# from scratch on only the 18 final taxa. Running IQ-TREE directly on a
# fasta/partition file already reduced to 18 taxa (e.g. the output of
# build_reduced_supermatrix.py) produces a DIFFERENT tree: bootstrap support
# is re-estimated from just 18 taxa instead of inherited from the 49-taxon
# analysis, and the two branches that sit below 100% support in the full
# 49-taxon tree (92%, 98%) come out at 100% purely from the taxon-count
# reduction -- with or without removing any locus. Always prune; never
# rebuild directly on 18 taxa if the result needs to be compared to a
# figure or number from this manuscript.
#
# Usage:
#   Rscript prune_to_final_taxa.R <49taxon_tree.contree> <taxa_codes.txt> <output.nwk>
#
# <taxa_codes.txt> is the same file used by build_reduced_supermatrix.py:
# one sample code per line (e.g. S102A10); each code must be a unique
# suffix of exactly one tip label in the input tree.
# ==============================================================================

suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: prune_to_final_taxa.R <49taxon_tree> <taxa_codes.txt> <output.nwk>")
tree_file  <- args[1]
codes_file <- args[2]
out_file   <- args[3]

target_codes <- readLines(codes_file)
target_codes <- target_codes[nzchar(target_codes)]

tr <- read.tree(tree_file)
cat("Input tree: Ntip =", Ntip(tr), "\n")

code_to_tip <- sapply(target_codes, function(code) {
  hits <- tr$tip.label[endsWith(tr$tip.label, code)]
  if (length(hits) != 1) {
    cat("WARN: code", code, "matched", length(hits), "tip(s):",
        paste(hits, collapse = ", "), "\n")
    return(NA)
  }
  hits
})

if (anyNA(code_to_tip)) {
  stop("Not every code in ", codes_file, " matched exactly one tip. See warnings above.")
}

pruned <- keep.tip(tr, unname(code_to_tip))
cat("Pruned tree: Ntip =", Ntip(pruned), "\n")

write.tree(pruned, out_file)
cat("Wrote:", out_file, "\n")

if (!is.null(pruned$node.label)) {
  cat("\nSupport-value distribution:\n")
  print(table(pruned$node.label))
}
