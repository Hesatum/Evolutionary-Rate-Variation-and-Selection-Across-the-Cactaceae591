#!/usr/bin/env Rscript
# ==============================================================================
# make_figure1.R
#
# Reproducible plotting script for the Cereus/Cipocereus phylogeny (Figure 1
# of "Evolutionary Rate Variation and Selection Across the Cactaceae591
# Probe Set"), reused here for the sensitivity trees in Figure S2 (loci
# under selection removed).
#
# Input:  a rooted or unrooted Newick tree with 18 tips matching the sample
#         codes in Table 1 (branch lengths required; UFBoot node.labels
#         optional -- if absent, only branch lengths are drawn).
# Output: <prefix>.eps  (vector, journal-compliant), <prefix>.pdf (vector
#         preview), <prefix>.png (raster preview, 600 dpi).
#
# Method: R v4.4.1, ape v5.8.1. Tree rooted on the Cipocereus outgroup with
# ape::root(); tips reordered with ape::rotateConstr() to the clade order
# used in the original figure; plotted with ape::plot.phylo() using
# equal (non-proportional) node spacing so that node/branch labels never
# collide regardless of how short a given branch is -- the true branch
# length is still reported as a text label next to each edge.
#
# Usage:
#   Single tree:  Rscript make_figure1.R single <input.nwk> <output_prefix> "<panel title>"
#   Mirrored pair: Rscript make_figure1.R mirror <left.nwk> <right.nwk> <output_prefix> "<left title>" "<right title>"
# ==============================================================================

suppressPackageStartupMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]

# ---- fixed reference data (Table 1) -----------------------------------------
outgroup <- c("Cipocereuslaniflorus_S174A6", "Cipocereusminensis_S153A1")

order_species <- c(
  "Cereus_jamacaru_S180V2", "Cereushexagonus_S155A2", "Cereus_bicolor_S102A10",
  "Cereus_trigonodendron_S184A4",
  "Cereusfernambucensissericifer_S88F2", "Cereusfernambucensis_S80F9", "Cereussp.Nov_S169A2",
  "Cereusspegazzinii_S77A31", "Cereus_spegazzini_S180A14", "Cereus_vargasianus_S184A5",
  "Cereusrepandus_S162VA22", "Cereus_mortensenii_S184A3",
  "Cereussaddianus_S103D4", "Cereuskroenleinii_S149V10",
  "Cereusalbicaulis_S143F5", "Cereus_mirabella_S162A26",
  "Cipocereusminensis_S153A1", "Cipocereuslaniflorus_S174A6"
)

clades <- list(
  A1 = c("Cereus jamacaru S180V2","Cereus hexagonus S155A2","Cereus bicolor S102A10","Cereus trigonodendron S184A4"),
  A2 = c("Cereus sericifer S88F2","Cereus fernambucensis S80F9","Cereus ingens S169A2"),
  B  = c("Cereus spegazzinii S77A31","Cereus spegazzinii S180A14","Cereus vargasianus S184A5"),
  E  = c("Cereus repandus S162VA22","Cereus mortensenii S184A3"),
  C  = c("Cereus saddianus S103D4","Cereus phatnospermus S149V10"),
  D  = c("Cereus albicaulis S143F5","Cereus mirabella S162A26"),
  Outgroup = c("Cipocereus minensis S153A1","Cipocereus laniflorus S174A6")
)

pretty_label <- function(labs) {
  labs <- gsub("^Cereus_", "Cereus ", labs)
  labs <- gsub("^Cereus(?![_ ])", "Cereus ", labs, perl = TRUE)
  labs <- gsub("^Cipocereus", "Cipocereus ", labs)
  labs <- gsub("fernambucensissericifer", "sericifer", labs)
  labs <- gsub("spegazzini(?!i)", "spegazzinii", labs, perl = TRUE)
  labs <- gsub("kroenleinii", "phatnospermus", labs)
  labs <- gsub("sp\\.Nov", "ingens", labs)
  labs <- gsub("_", " ", labs)
  labs
}

prepare_tree <- function(tree_file) {
  tr <- read.tree(tree_file)
  stopifnot(Ntip(tr) == 18)
  tr <- root(tr, outgroup = outgroup, resolve.root = TRUE)
  tr <- rotateConstr(tr, constraint = rev(order_species))
  tr$tip.label <- pretty_label(tr$tip.label)
  tr
}

# Draw one tree panel. direction: "rightwards" or "leftwards".
# Node spacing is EQUAL (use.edge.length = FALSE) so short branches never
# cause label collisions; true branch lengths are printed as text instead.
draw_tree <- function(tr, direction, panel_title, label_side_pad) {
  has_support <- !is.null(tr$node.label) && any(nzchar(tr$node.label))
  node_lab <- if (has_support) tr$node.label else rep(NA, tr$Nnode)
  node_lab[node_lab %in% c("", "Root", "NA")] <- NA

  par(mar = c(1.2, label_side_pad[1], if (nzchar(panel_title)) 2 else 0.5, label_side_pad[2]),
      family = "Times", xpd = NA)

  # TRUE, proportional branch lengths -- layout and text both use the real
  # values (no transform). Collisions are avoided purely by giving the plot
  # a large absolute canvas (see fig_w_in/x.lim below), not by distorting
  # the tree geometry.
  bl <- tr$edge.length
  max_depth <- max(node.depth.edgelength(tr))

  plot(tr, direction = direction, show.tip.label = TRUE, font = 3,
       cex = 0.82, label.offset = max_depth * 0.035, edge.width = 1.1,
       y.lim = c(0.9, Ntip(tr) + 0.7), x.lim = c(0, max_depth * 1.95))
  if (nzchar(panel_title)) title(main = panel_title, cex.main = 1.05, font.main = 2)

  # bootstrap support: above the branch, close to the node
  if (has_support) {
    nodelabels(node_lab, frame = "none", adj = c(1.1, -0.7), cex = 0.62, font = 1)
  }

  # branch lengths: below the branch, at the edge midpoint -- blank for
  # zero-length (artificial root) edges so no "0.0000" is printed
  bl_txt <- ifelse(is.na(bl) | bl < 1e-6, "", sprintf("%.4f", bl))
  edgelabels(bl_txt, frame = "none", adj = c(0.5, 1.7), cex = 0.5, col = "grey35")

  # clade labels -- plain text, no box, placed just outside the node
  for (cn in names(clades)) {
    tips <- clades[[cn]][clades[[cn]] %in% tr$tip.label]
    if (length(tips) < 2) next
    mrca_node <- getMRCA(tr, tips)
    if (is.null(mrca_node)) next
    adj_x <- if (direction == "rightwards") -0.3 else 1.3
    nodelabels(cn, mrca_node, frame = "none", cex = 0.72, font = 2, adj = c(adj_x, 0.5))
  }
}

mm_to_in <- function(mm) mm / 25.4

if (mode == "single") {
  tree_file  <- args[2]
  out_prefix <- args[3]
  panel_title <- if (length(args) >= 4) args[4] else ""

  tr <- prepare_tree(tree_file)

  fig_w_in <- mm_to_in(220)   # extra lateral room so labels are never clipped
  fig_h_in <- mm_to_in(155)

  plot_it <- function() draw_tree(tr, "rightwards", panel_title, label_side_pad = c(0.3, 11))

  setEPS()
  postscript(paste0(out_prefix, ".eps"), width = fig_w_in, height = fig_h_in,
             family = "Times", horizontal = FALSE, paper = "special")
  plot_it(); dev.off()

  pdf(paste0(out_prefix, ".pdf"), width = fig_w_in, height = fig_h_in, family = "Times")
  plot_it(); dev.off()

  png(paste0(out_prefix, ".png"), width = fig_w_in * 600, height = fig_h_in * 600, res = 600)
  plot_it(); dev.off()

  cat("Wrote:", paste0(out_prefix, c(".eps", ".pdf", ".png")), "\n")
  cat(sprintf("Panel size: %.1f x %.1f mm\n", fig_w_in * 25.4, fig_h_in * 25.4))

} else if (mode == "mirror") {
  left_file   <- args[2]
  right_file  <- args[3]
  out_prefix  <- args[4]
  left_title  <- if (length(args) >= 5) args[5] else ""
  right_title <- if (length(args) >= 6) args[6] else ""

  t_left  <- prepare_tree(left_file)
  t_right <- prepare_tree(right_file)

  fig_w_in <- mm_to_in(240)
  fig_h_in <- mm_to_in(150)

  plot_it <- function() {
    par(mfrow = c(1, 2))
    draw_tree(t_left,  "rightwards", left_title,  label_side_pad = c(0.3, 8.5))
    draw_tree(t_right, "leftwards",  right_title, label_side_pad = c(8.5, 0.3))
  }

  setEPS()
  postscript(paste0(out_prefix, ".eps"), width = fig_w_in, height = fig_h_in,
             family = "Times", horizontal = FALSE, paper = "special")
  plot_it(); dev.off()

  pdf(paste0(out_prefix, ".pdf"), width = fig_w_in, height = fig_h_in, family = "Times")
  plot_it(); dev.off()

  png(paste0(out_prefix, ".png"), width = fig_w_in * 600, height = fig_h_in * 600, res = 600)
  plot_it(); dev.off()

  cat("Wrote:", paste0(out_prefix, c(".eps", ".pdf", ".png")), "\n")
  cat(sprintf("Panel size: %.1f x %.1f mm\n", fig_w_in * 25.4, fig_h_in * 25.4))

} else {
  stop("First argument must be 'single' or 'mirror'")
}
