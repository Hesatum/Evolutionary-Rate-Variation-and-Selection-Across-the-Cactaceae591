#!/usr/bin/env bash
# run_iqtree.sh
#
# Species-tree estimation used for the main tree and for every sensitivity
# tree in this study (manuscript Methods 2.3): ModelFinder Plus for
# per-partition substitution model selection, 1000 ultrafast bootstrap
# (UFBoot2) replicates, IQ-TREE 2.
#
# Usage:
#   ./run_iqtree.sh <supermatrix.fasta> <partitions.txt> <output_prefix> [seed]
#
# Example (exon-selection-removed tree):
#   ./run_iqtree.sh trees/exon_selection_removed/supermatrix_18tax.fasta \
#                    trees/exon_selection_removed/supermatrix_18tax_partitions.txt \
#                    trees/exon_selection_removed/run_exon_selection_removed \
#                    141309
#
# The seed used for the sensitivity re-runs behind Figure S2 was 141309;
# omit it to let IQ-TREE pick one at random.

set -euo pipefail

FASTA="$1"
PARTITIONS="$2"
PREFIX="$3"
SEED="${4:-}"

SEED_ARG=()
if [ -n "$SEED" ]; then
  SEED_ARG=(-seed "$SEED")
fi

iqtree2 -s "$FASTA" -p "$PARTITIONS" -m MFP -bb 1000 -nt AUTO -pre "$PREFIX" "${SEED_ARG[@]}"

echo "Tree: ${PREFIX}.treefile"
echo "Support-annotated consensus: ${PREFIX}.contree"
echo "Full report: ${PREFIX}.iqtree"
