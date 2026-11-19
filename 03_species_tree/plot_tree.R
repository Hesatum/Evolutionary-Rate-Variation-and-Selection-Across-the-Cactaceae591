#!/usr/bin/env Rscript
# ==============================================================================
# plot_tree.R
#
# Plots the Cereus/Cipocereus phylogeny in the style of Figure 1 (single
# tree) and Figure S3 (two trees, mirrored).
#
# Input:  a Newick tree with the 18 tips of Table 1 (branch lengths required;
#         UFBoot node labels optional).
# Output: <prefix>.eps (vector), <prefix>.pdf (vector), <prefix>.png (800 dpi).
#
# R v4.4.1, ape v5.8.1. The tree is rooted on the Cipocereus outgroup with
# ape::root() and tips are ordered with ape::rotateConstr(). The support of
# the root bipartition is shown on both descendants of the root. The
# zero-length outgroup branch at the root is drawn with a short display-only
# length (ROOT_STUB_FRACTION); its branch-length label is left blank.
#
# Usage:
#   Single tree:   Rscript plot_tree.R single <input.nwk> <output_prefix> "<panel title>"
#   Mirrored pair: Rscript plot_tree.R mirror <left.nwk> <right.nwk> <output_prefix> "<left title>" "<right title>"
# ==============================================================================

suppressPackageStartupMessages(library(ape))

if (.Platform$OS.type == "windows") windowsFonts(Times = windowsFont("Times New Roman"))

args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]

# Display-only length for the zero-length root branch, as a fraction of
# tree depth. 0 draws the true length.
ROOT_STUB_FRACTION <- 0.035

SINGLE_W_MM <- 168; SINGLE_H_MM <- 150; SINGLE_POINTSIZE <- 10
MIRROR_W_MM <- 168; MIRROR_H_MM <- 140; MIRROR_POINTSIZE <- 8
PNG_DPI <- 800

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
  tr <- root(tr, outgroup = outgroup, resolve.root = TRUE, edgelabel = TRUE)
  if (!is.null(tr$node.label)) {
    root_kids <- tr$edge[tr$edge[, 1] == Ntip(tr) + 1, 2] - Ntip(tr)
    root_kids <- root_kids[root_kids > 0]
    root_sup <- tr$node.label[root_kids]
    root_sup <- root_sup[nzchar(root_sup)]
    if (length(root_kids) == 2 && length(root_sup) == 1) tr$node.label[root_kids] <- root_sup
  }
  tr <- rotateConstr(tr, constraint = rev(order_species))
  tr$tip.label <- pretty_label(tr$tip.label)

  # True edge lengths are kept for the branch-length labels.
  tr$edge.length.true <- tr$edge.length
  if (ROOT_STUB_FRACTION > 0) {
    root_id <- Ntip(tr) + 1
    root_edges <- which(tr$edge[, 1] == root_id)
    zero_root_edges <- root_edges[tr$edge.length[root_edges] < 1e-8]
    if (length(zero_root_edges) > 0) {
      stub <- ROOT_STUB_FRACTION * max(node.depth.edgelength(tr))
      tr$edge.length[zero_root_edges] <- stub
    }
  }
  tr
}

# Draw one tree panel. direction: "rightwards" or "leftwards".
draw_tree <- function(tr, direction, panel_title, label_side_pad, x_lim_mult = 1.95, xaxs = "r") {
  has_support <- !is.null(tr$node.label) && any(nzchar(tr$node.label))
  node_lab <- if (has_support) tr$node.label else rep(NA, tr$Nnode)
  node_lab[node_lab %in% c("", "Root", "NA")] <- NA

  par(mar = c(1.2, label_side_pad[1], if (nzchar(panel_title)) 2 else 0.5, label_side_pad[2]),
      family = "Times", xpd = NA, xaxs = xaxs)

  bl <- if (!is.null(tr$edge.length.true)) tr$edge.length.true else tr$edge.length
  max_depth <- max(node.depth.edgelength(tr))

  tip_cex <- 0.82
  lab_off <- max_depth * 0.035
  plot(tr, direction = direction, show.tip.label = TRUE, font = 3,
       tip.color = "transparent",
       cex = tip_cex, label.offset = lab_off, edge.width = 1.1,
       y.lim = c(0.9, Ntip(tr) + 0.7),
       x.lim = if (is.null(x_lim_mult)) NULL else c(0, max_depth * x_lim_mult))
  if (nzchar(panel_title)) title(main = panel_title, cex.main = 1.05, font.main = 2)

  lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  right <- direction == "rightwards"
  tip_expr <- as.expression(lapply(tr$tip.label, function(l) {
    sp <- sub(" [^ ]+$", "", l); code <- sub("^.* ", "", l)
    bquote(italic(.(sp)) * " " * .(code))
  }))
  n <- Ntip(tr)
  text(lp$xx[1:n] + if (right) lab_off else -lab_off, lp$yy[1:n], tip_expr,
       adj = c(if (right) 0 else 1, 0.5), cex = tip_cex)

  # bootstrap support: above the branch, close to the node
  if (has_support) {
    nodelabels(node_lab, frame = "none", adj = c(1.1, -0.7), cex = 0.62, font = 1)
  }

  # branch lengths: below the branch, blank for zero-length edges
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

  fig_w_in <- mm_to_in(SINGLE_W_MM)
  fig_h_in <- mm_to_in(SINGLE_H_MM)
  ps <- SINGLE_POINTSIZE

  plot_it <- function() draw_tree(tr, "rightwards", panel_title, label_side_pad = c(0.3, 1.2),
                                  x_lim_mult = NULL, xaxs = "i")

  setEPS()
  postscript(paste0(out_prefix, ".eps"), width = fig_w_in, height = fig_h_in,
             family = "Times", horizontal = FALSE, paper = "special", pointsize = ps)
  plot_it(); dev.off()

  pdf(paste0(out_prefix, ".pdf"), width = fig_w_in, height = fig_h_in, family = "Times", pointsize = ps)
  plot_it(); dev.off()

  png(paste0(out_prefix, ".png"), width = fig_w_in * PNG_DPI, height = fig_h_in * PNG_DPI, res = PNG_DPI, pointsize = ps)
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

  fig_w_in <- mm_to_in(MIRROR_W_MM)
  fig_h_in <- mm_to_in(MIRROR_H_MM)
  ps <- MIRROR_POINTSIZE

  plot_it <- function() {
    par(mfrow = c(1, 2))
    draw_tree(t_left,  "rightwards", left_title,  label_side_pad = c(0.3, 8.5))
    draw_tree(t_right, "leftwards",  right_title, label_side_pad = c(8.5, 0.3))
  }

  setEPS()
  postscript(paste0(out_prefix, ".eps"), width = fig_w_in, height = fig_h_in,
             family = "Times", horizontal = FALSE, paper = "special", pointsize = ps)
  plot_it(); dev.off()

  pdf(paste0(out_prefix, ".pdf"), width = fig_w_in, height = fig_h_in, family = "Times", pointsize = ps)
  plot_it(); dev.off()

  png(paste0(out_prefix, ".png"), width = fig_w_in * PNG_DPI, height = fig_h_in * PNG_DPI, res = PNG_DPI, pointsize = ps)
  plot_it(); dev.off()

  cat("Wrote:", paste0(out_prefix, c(".eps", ".pdf", ".png")), "\n")
  cat(sprintf("Panel size: %.1f x %.1f mm\n", fig_w_in * 25.4, fig_h_in * 25.4))

} else {
  stop("First argument must be 'single' or 'mirror'")
}
